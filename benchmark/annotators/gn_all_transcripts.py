"""Fetch Genome Nexus's per-transcript annotations into an all_transcripts parquet.

The pick annotator (genome_nexus.py) keeps only the canonical transcript summary.
The same REST response also carries annotation_summary.transcriptConsequenceSummaries
— one normalized record per overlapping transcript (transcriptId, hgvsp/hgvspShort,
consequenceTerms, variantClassification, …) — which is exactly the Stage 1 view.
This re-queries GN saving all of them, checkpointed per batch so it can resume.
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


def _prot_pos(pp: dict | None) -> str:
    pp = pp or {}
    s, e = pp.get("start"), pp.get("end")
    if s is not None and e is not None and s != e:
        return f"{s}-{e}"
    return "" if s is None else str(s)


def rows_from_record(rec: dict) -> list[dict]:
    var = rec.get("originalVariantQuery", "")
    parts = var.split(",")
    key = ":".join(parts) if len(parts) == 5 else ""
    if not key:
        return []
    summ = (rec.get("annotation_summary") or {}).get("transcriptConsequenceSummaries") or []
    out = []
    for t in summ:
        tx = strip_transcript_version(t.get("transcriptId") or "")
        if not tx:
            continue
        out.append(
            {
                "var_key": key,
                "hugo_symbol": t.get("hugoGeneSymbol") or "",
                "transcript_id": tx,
                "consequence": t.get("consequenceTerms") or "",
                "variant_classification": t.get("variantClassification") or "",
                "hgvsc": (t.get("hgvsc") or "").split(":", 1)[-1] if t.get("hgvsc") else "",
                "hgvsp": t.get("hgvsp") or "",
                "hgvsp_short": t.get("hgvspShort") or "",
                "protein_position": _prot_pos(t.get("proteinPosition")),
                "codons": t.get("codonChange") or "",
            }
        )
    return out


def annotate_batch(rows, session, api, retries=4):
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
    last = None
    for attempt in range(retries):
        try:
            resp = session.post(
                api, params={"fields": "annotation_summary"}, json=payload,
                headers={"Accept": "application/json"}, timeout=180,
            )
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"batch failed after {retries} retries: {last}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--api", default="https://www.genomenexus.org/annotation/genomic")
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--batch-size", type=int, default=500)
    args = ap.parse_args()

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "all_transcripts_cache.jsonl"

    rows = pd.read_parquet(args.input).to_dict("records")
    done: set[str] = set()
    if cache.exists():
        with cache.open() as fh:
            for line in fh:
                done.add(json.loads(line)["var_key"])
        print(f"resuming: {len(done):,} variants already cached", flush=True)

    todo = [r for r in rows if r["var_key"] not in done]
    print(f"to annotate: {len(todo):,} of {len(rows):,}", flush=True)

    session = requests.Session()
    t0 = time.time()
    with cache.open("a") as fh:
        for i in range(0, len(todo), args.batch_size):
            batch = todo[i : i + args.batch_size]
            resp = annotate_batch(batch, session, args.api)
            got: set[str] = set()
            for rec in resp:
                trs = rows_from_record(rec)
                for tr in trs:
                    got.add(tr["var_key"])
                    fh.write(json.dumps(tr) + "\n")
            # mark variants that returned nothing so resume skips them
            for r in batch:
                if r["var_key"] not in got:
                    fh.write(json.dumps({"var_key": r["var_key"], "_empty": True}) + "\n")
            fh.flush()
            n = i + len(batch)
            print(f"  {n:,}/{len(todo):,}  ({n/(time.time()-t0):.0f} var/s)", flush=True)

    records = []
    with cache.open() as fh:
        for line in fh:
            rec = json.loads(line)
            if not rec.get("_empty"):
                records.append(rec)
    out = pd.DataFrame(records)
    out = out[out["transcript_id"] != ""][STD_FIELDS]
    out = out.drop_duplicates(subset=["var_key", "transcript_id"]).reset_index(drop=True)
    out.to_parquet(out_dir / "all_transcripts.parquet", index=False)
    print(
        f"wrote {out_dir/'all_transcripts.parquet'}  ({len(out):,} rows, "
        f"{out['var_key'].nunique():,} variants)"
    )


if __name__ == "__main__":
    main()
