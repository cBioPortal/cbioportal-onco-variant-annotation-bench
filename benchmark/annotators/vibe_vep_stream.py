"""Run vibe-vep in stream mode and collect every transcript per variant.

`vibe-vep annotate vcf` merges same-position multi-allelic sites into a single
VCF record with ALT=A,B and drops all but the first record's ID column, so the
downstream parser silently loses one of every SNV pair. `annotate stream` reads
one JSON per line and emits one JSON per line — no merging, no ID collision.

Input JSONL uses Genome Nexus's genomic-location schema (`chromosome`, `start`,
`end`, `referenceAllele`, `variantAllele`), which is exactly the schema of our
`input.parquet`. Output is VEP-compatible JSONL: `input` (echoing the input
coordinates), `most_severe_consequence`, `transcript_consequences[]`.

Writes `all_transcripts.parquet` (long: one row per var_key × transcript).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (  # noqa: E402
    STD_FIELDS,
    hgvsp_to_short,
    strip_protein_prefix,
    strip_transcript_version,
    var_key,
)


def protein_position(tc: dict) -> str:
    s, e = tc.get("protein_start"), tc.get("protein_end")
    if s is None:
        return ""
    if e is not None and e != s:
        return f"{s}-{e}"
    return str(s)


def tc_to_std(vk: str, tc: dict) -> dict:
    return {
        "var_key": vk,
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


def stream_annotate(
    input_rows: list[dict],
    vibe_vep_bin: str,
    assembly: str,
) -> list[dict]:
    """Spawn `vibe-vep annotate stream` and one-shot every input variant.

    Uses `Popen.communicate(input=...)` (threaded writer + reader) to avoid
    the classic stdin/stdout pipe-buffer deadlock — feeding 58k lines
    synchronously before reading stdout would fill both pipe buffers and
    hang the subprocess.
    """
    cmd = [
        vibe_vep_bin, "annotate", "stream",
        "--assembly", assembly,
        "--output-format", "ensembl-vep-jsonl",
    ]
    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    lines = []
    for r in input_rows:
        lines.append(json.dumps({
            "chromosome": str(r["chromosome"]),
            "start": int(r["start"]),
            "end": int(r["end"]),
            "referenceAllele": r["reference_allele"],
            "variantAllele": r["variant_allele"],
        }))
    payload = "\n".join(lines) + "\n"
    out, err = proc.communicate(input=payload)
    if proc.returncode != 0:
        raise RuntimeError(f"vibe-vep stream failed (rc={proc.returncode}): {err[-2000:]}")
    return [json.loads(line) for line in out.splitlines() if line.strip()]


def build_rows(input_rows: list[dict], responses: list[dict]) -> list[dict]:
    """Zip inputs and stream responses on position order (stream preserves it)."""
    if len(responses) != len(input_rows):
        # Fall back to matching on `input` field if lengths disagree.
        print(
            f"WARN: {len(responses)} responses vs {len(input_rows)} inputs; matching by echoed input",
            file=sys.stderr,
        )
    rows: list[dict] = []
    for r, resp in zip(input_rows, responses):
        vk = r["var_key"]
        tcs = resp.get("transcript_consequences") or []
        for tc in tcs:
            rows.append(tc_to_std(vk, tc))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="input.parquet")
    ap.add_argument("--out", required=True, help="output all_transcripts.parquet")
    ap.add_argument("--assembly", default="GRCh37")
    ap.add_argument(
        "--vibe-vep",
        default="/Users/lix2/Documents/GitHub/vibe-vep/vibe-vep",
        help="path to the vibe-vep binary",
    )
    ap.add_argument("--limit", type=int, default=None, help="annotate first N variants (dev)")
    args = ap.parse_args()

    inp = pd.read_parquet(args.input)
    if args.limit:
        inp = inp.head(args.limit)
    input_rows = inp.to_dict("records")
    # Ensure every input has a var_key derivable from chr/start/end/ref/alt.
    for r in input_rows:
        if "var_key" not in r:
            r["var_key"] = var_key(
                r["chromosome"], r["start"], r["end"],
                r["reference_allele"], r["variant_allele"],
            )

    t0 = time.time()
    print(f"streaming {len(input_rows):,} variants through vibe-vep...", flush=True)
    responses = stream_annotate(input_rows, args.vibe_vep, args.assembly)
    print(f"  got {len(responses):,} responses in {time.time()-t0:.1f}s", flush=True)

    rows = build_rows(input_rows, responses)
    df = pd.DataFrame(rows)
    df = df[df["transcript_id"] != ""][STD_FIELDS]
    df = df.drop_duplicates(subset=["var_key", "transcript_id"]).reset_index(drop=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.out, index=False)
    n_var = df["var_key"].nunique()
    print(
        f"wrote {args.out}  ({len(df):,} transcript rows, "
        f"{n_var:,} variants, {len(df)/max(1,n_var):.1f} tx/variant)"
    )


if __name__ == "__main__":
    main()
