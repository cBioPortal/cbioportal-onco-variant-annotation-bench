"""Annotate variants with the Genome Nexus public REST API.

Genome Nexus is MSK's production pipeline and the source of the truth MAF, so it
serves as the baseline sanity check: concordance should be ~100%. It also
exercises the whole extract -> annotate -> compare framework end to end without
any local install.

Output: data/tools/genome_nexus/annotations.parquet in STD_FIELDS format.
The run is checkpointed per batch to a JSONL cache so it can resume.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, strip_transcript_version  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
IN_PARQUET = ROOT / "data" / "truth" / "input.parquet"
OUT_DIR = ROOT / "data" / "tools" / "genome_nexus"
API = "https://www.genomenexus.org/annotation/genomic"


def to_std(rec: dict) -> dict:
    """Map a GN annotation response to the normalized picked-transcript record."""
    var = rec.get("originalVariantQuery", "")
    key = ""
    if var:
        # originalVariantQuery is "chr,start,end,ref,alt"
        parts = var.split(",")
        if len(parts) == 5:
            key = f"{parts[0]}:{parts[1]}:{parts[2]}:{parts[3]}:{parts[4]}"

    summ = (rec.get("annotation_summary") or {}).get(
        "transcriptConsequenceSummary"
    ) or {}
    pp = summ.get("proteinPosition") or {}
    prot_start = pp.get("start")
    prot_end = pp.get("end")
    if prot_start is not None and prot_end is not None and prot_start != prot_end:
        prot_pos = f"{prot_start}-{prot_end}"
    elif prot_start is not None:
        prot_pos = str(prot_start)
    else:
        prot_pos = ""

    return {
        "var_key": key,
        "hugo_symbol": summ.get("hugoGeneSymbol", "") or "",
        "transcript_id": strip_transcript_version(summ.get("transcriptId", "") or ""),
        "consequence": summ.get("consequenceTerms", "") or "",
        "variant_classification": summ.get("variantClassification", "") or "",
        "hgvsc": summ.get("hgvsc", "") or "",
        "hgvsp": summ.get("hgvsp", "") or "",
        "hgvsp_short": summ.get("hgvspShort", "") or "",
        "protein_position": prot_pos,
        "codons": summ.get("codonChange", "") or "",
    }


def annotate_batch(rows: list[dict], session: requests.Session, retries: int = 4) -> list[dict]:
    payload = [
        {
            "chromosome": str(r["chromosome"]),
            "start": int(r["start"]),
            "end": int(r["end"]),
            "referenceAllele": r["reference_allele"],
            "variantAllele": r["variant_allele"],
        }
        for r in rows
    ]
    last_err = None
    for attempt in range(retries):
        try:
            resp = session.post(
                API,
                params={"fields": "annotation_summary"},
                json=payload,
                headers={"Accept": "application/json"},
                timeout=120,
            )
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"batch failed after {retries} retries: {last_err}")


def main() -> None:
    global API
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--limit", type=int, default=None, help="annotate first N variants (dev)")
    ap.add_argument("--api", default=API, help="GN base annotation endpoint")
    ap.add_argument("--input", default=str(IN_PARQUET), help="input.parquet")
    ap.add_argument("--outdir", default=str(OUT_DIR), help="output dir")
    args = ap.parse_args()
    API = args.api

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "cache.jsonl"

    inp = pd.read_parquet(args.input)
    if args.limit:
        inp = inp.head(args.limit)
    rows = inp.to_dict("records")

    done_keys: set[str] = set()
    if cache.exists():
        with cache.open() as fh:
            for line in fh:
                rec = json.loads(line)
                done_keys.add(rec["var_key"])
        print(f"resuming: {len(done_keys):,} variants already cached")

    todo = [r for r in rows if r["var_key"] not in done_keys]
    print(f"to annotate: {len(todo):,} of {len(rows):,}")

    session = requests.Session()
    t0 = time.time()
    with cache.open("a") as fh:
        for i in range(0, len(todo), args.batch_size):
            batch = todo[i : i + args.batch_size]
            resp = annotate_batch(batch, session)
            std = [to_std(rec) for rec in resp]
            # Any input the API silently dropped gets an empty (unannotated) row.
            got = {s["var_key"] for s in std if s["var_key"]}
            for r in batch:
                if r["var_key"] not in got:
                    std.append({**{f: "" for f in STD_FIELDS}, "var_key": r["var_key"]})
            for s in std:
                fh.write(json.dumps(s) + "\n")
            fh.flush()
            n = i + len(batch)
            rate = n / (time.time() - t0)
            print(f"  {n:,}/{len(todo):,}  ({rate:.0f} var/s)", flush=True)

    # Consolidate cache -> parquet (dedup, keep annotated row if present).
    records: dict[str, dict] = {}
    with cache.open() as fh:
        for line in fh:
            rec = json.loads(line)
            k = rec["var_key"]
            if k not in records or (not records[k]["hugo_symbol"] and rec["hugo_symbol"]):
                records[k] = rec
    out = pd.DataFrame(list(records.values()))[STD_FIELDS]
    out.to_parquet(out_dir / "annotations.parquet", index=False)
    print(f"wrote {out_dir/'annotations.parquet'}  ({len(out):,} variants)")


if __name__ == "__main__":
    main()
