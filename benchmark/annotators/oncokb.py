"""Annotate truth variants with OncoKB oncogenicity (byGenomicChange).

Produces a per-variant flag used to filter the whole site to oncogenic variants.
Not a benchmarked tool — a variant-level label joined onto every view. Oncogenicity
is OncoKB's curated call; "oncogenic" here means the label is one of
Oncogenic / Likely Oncogenic / Resistance (Predicted Oncogenic for older data).

The API token is read from a local file (default data/oncokb/token.txt) or --token,
never stored in the output. Checkpointed per batch to a JSONL cache to resume.
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

API = "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
ONCOGENIC = {"Oncogenic", "Likely Oncogenic", "Resistance", "Predicted Oncogenic"}


def genomic_location(r: dict) -> str:
    chrom = str(r["chromosome"]).replace("chr", "")
    return f"{chrom},{int(r['start'])},{int(r['end'])},{r['reference_allele']},{r['variant_allele']}"


def annotate_batch(rows, session, token, ref_genome, retries=4):
    payload = [
        {"genomicLocation": genomic_location(r), "referenceGenome": ref_genome}
        for r in rows
    ]
    last = None
    for attempt in range(retries):
        try:
            resp = session.post(
                API,
                headers={"Authorization": f"Bearer {token}", "Content-Type": "application/json"},
                json=payload, timeout=180,
            )
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"batch failed after {retries} retries: {last}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ref-genome", required=True, choices=["GRCh37", "GRCh38"])
    ap.add_argument("--token-file", default="data/oncokb/token.txt")
    ap.add_argument("--token", default=None)
    ap.add_argument("--batch-size", type=int, default=800)
    args = ap.parse_args()

    token = args.token or Path(args.token_file).read_text().strip()
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "oncokb_cache.jsonl"

    rows = pd.read_parquet(args.input).to_dict("records")
    done: set[str] = set()
    if cache.exists():
        with cache.open() as fh:
            for line in fh:
                done.add(json.loads(line)["var_key"])
        print(f"resuming: {len(done):,} cached", flush=True)
    todo = [r for r in rows if r["var_key"] not in done]
    print(f"to annotate: {len(todo):,} of {len(rows):,}", flush=True)

    session = requests.Session()
    t0 = time.time()
    with cache.open("a") as fh:
        for i in range(0, len(todo), args.batch_size):
            batch = todo[i : i + args.batch_size]
            resp = annotate_batch(batch, session, token, args.ref_genome)
            for r, res in zip(batch, resp):
                label = (res or {}).get("oncogenic") or ""
                fh.write(json.dumps({
                    "var_key": r["var_key"],
                    "oncokb_oncogenic": label,
                    "oncokb_is_oncogenic": label in ONCOGENIC,
                    "oncokb_gene": ((res or {}).get("query") or {}).get("hugoSymbol") or "",
                    "oncokb_effect": ((res or {}).get("mutationEffect") or {}).get("knownEffect") or "",
                    "oncokb_highest_level": (res or {}).get("highestSensitiveLevel") or "",
                }) + "\n")
            fh.flush()
            n = i + len(batch)
            print(f"  {n:,}/{len(todo):,}  ({n/(time.time()-t0):.0f} var/s)", flush=True)

    records = [json.loads(l) for l in cache.open()]
    out = pd.DataFrame(records).drop_duplicates(subset="var_key", keep="last").reset_index(drop=True)
    out.to_parquet(out_dir / "oncokb.parquet", index=False)
    n_onc = int(out["oncokb_is_oncogenic"].sum())
    print(f"wrote {out_dir/'oncokb.parquet'}  ({len(out):,} variants, {n_onc:,} oncogenic = {n_onc/max(1,len(out)):.1%})")


if __name__ == "__main__":
    main()
