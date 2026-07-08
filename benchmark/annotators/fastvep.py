"""Parse fastVEP JSON output into the benchmark's standard formats.

fastVEP (github.com/Huang-lab/fastVEP) is run twice on the same VCF:
  pick.json  --pick  -> one chosen transcript per variant (Benchmark 2: pick)
  all.json   (none)  -> every overlapping transcript (Benchmark 1: completeness)

This script converts:
  pick.json -> data/tools/fastvep/annotations.parquet   (STD_FIELDS, one/variant)
  all.json  -> data/tools/fastvep/all_transcripts.parquet (long: var_key,transcript,...)

The VCF ID column carries the var_key, echoed back as JSON `id`, so output joins
straight to truth with no coordinate re-derivation.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (  # noqa: E402
    STD_FIELDS,
    hgvsp_to_short,
    strip_protein_prefix,
    strip_transcript_version,
)

ROOT = Path(__file__).resolve().parent.parent.parent
DIR = ROOT / "data" / "tools" / "fastvep"


def protein_position(tc: dict) -> str:
    s, e = tc.get("protein_start"), tc.get("protein_end")
    if s is None:
        return ""
    if e is not None and e != s:
        return f"{s}-{e}"
    return str(s)


def tc_to_std(var_key: str, tc: dict) -> dict:
    return {
        "var_key": var_key,
        "hugo_symbol": tc.get("gene_symbol", "") or "",
        "transcript_id": strip_transcript_version(tc.get("transcript_id", "") or ""),
        "consequence": ",".join(tc.get("consequence_terms", []) or []),
        "variant_classification": tc.get("variant_classification", "") or "",
        "hgvsc": strip_protein_prefix(tc.get("hgvsc", "") or ""),
        "hgvsp": strip_protein_prefix(tc.get("hgvsp", "") or ""),
        "hgvsp_short": hgvsp_to_short(tc.get("hgvsp", "") or ""),
        "protein_position": protein_position(tc),
        "codons": tc.get("codons", "") or "",
    }


def load_json(path: Path) -> list[dict]:
    txt = path.read_text().strip()
    if not txt:
        return []
    try:
        d = json.loads(txt)
        return d if isinstance(d, list) else [d]
    except json.JSONDecodeError:
        return [json.loads(ln) for ln in txt.splitlines() if ln.strip()]


def build_pick(pick_json: Path) -> pd.DataFrame:
    rows = []
    for r in load_json(pick_json):
        key = r.get("id", "")
        tcs = r.get("transcript_consequences") or []
        if tcs:
            rows.append(tc_to_std(key, tcs[0]))
        else:
            # No transcript overlap: record the variant as annotated-but-no-tx,
            # keeping the most-severe consequence for the consequence metric.
            rows.append(
                {
                    **{f: "" for f in STD_FIELDS},
                    "var_key": key,
                    "consequence": r.get("most_severe_consequence", "") or "",
                }
            )
    return pd.DataFrame(rows)[STD_FIELDS]


def build_all(all_json: Path) -> pd.DataFrame:
    rows = []
    for r in load_json(all_json):
        key = r.get("id", "")
        for tc in r.get("transcript_consequences") or []:
            rows.append(tc_to_std(key, tc))
    return pd.DataFrame(rows)[STD_FIELDS]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pick", default=str(DIR / "pick.json"))
    ap.add_argument("--all", default=str(DIR / "all.json"))
    ap.add_argument("--outdir", default=str(DIR))
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pick = build_pick(Path(args.pick))
    pick.to_parquet(outdir / "annotations.parquet", index=False)
    print(f"wrote {outdir/'annotations.parquet'}  ({len(pick):,} variants)")

    allt = build_all(Path(args.all))
    allt.to_parquet(outdir / "all_transcripts.parquet", index=False)
    n_var = allt["var_key"].nunique()
    print(
        f"wrote {outdir/'all_transcripts.parquet'}  "
        f"({len(allt):,} transcript rows across {n_var:,} variants)"
    )


if __name__ == "__main__":
    main()
