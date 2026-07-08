"""Extract the truth set and annotator input from an MSK-IMPACT MAF.

Produces two artifacts under data/truth/:
  - truth.parquet   one row per UNIQUE variant, normalized to STD_FIELDS
  - input.parquet   minimal genomic coordinates to feed annotators

Annotation is deterministic per genomic locus, so we dedup the per-sample MAF
calls down to unique variants. The first occurrence's truth annotation is kept
(they are identical across samples for the same locus in practice).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import MAF_TRUTH_COLS, load_maf, var_key

ROOT = Path(__file__).resolve().parent.parent
TRUTH_DIR = ROOT / "data" / "truth"


def build(maf_path: str, nrows: int | None) -> tuple[pd.DataFrame, pd.DataFrame]:
    maf = load_maf(maf_path, nrows=nrows)
    missing = [c for c in MAF_TRUTH_COLS if c not in maf.columns]
    if missing:
        raise SystemExit(f"MAF missing expected columns: {missing}")

    df = maf[MAF_TRUTH_COLS].copy()
    df["var_key"] = [
        var_key(c, s, e, r, a)
        for c, s, e, r, a in zip(
            df["Chromosome"],
            df["Start_Position"],
            df["End_Position"],
            df["Reference_Allele"],
            df["Tumor_Seq_Allele2"],
        )
    ]

    n_calls = len(df)
    df = df.drop_duplicates(subset="var_key", keep="first").reset_index(drop=True)
    n_unique = len(df)

    truth = pd.DataFrame(
        {
            "var_key": df["var_key"],
            "hugo_symbol": df["Hugo_Symbol"],
            "transcript_id": df["Transcript_ID"],
            "consequence": df["Consequence"],
            "variant_classification": df["Variant_Classification"],
            "hgvsc": df["HGVSc"],
            "hgvsp": df["HGVSp"],
            "hgvsp_short": df["HGVSp_Short"],
            "protein_position": df["Protein_position"],
            "codons": df["Codons"],
        }
    )

    inp = pd.DataFrame(
        {
            "var_key": df["var_key"],
            "chromosome": df["Chromosome"].str.replace("chr", "", regex=False),
            "start": df["Start_Position"],
            "end": df["End_Position"],
            "reference_allele": df["Reference_Allele"],
            "variant_allele": df["Tumor_Seq_Allele2"],
            "variant_type": df["Variant_Type"],
        }
    )

    print(f"MAF calls: {n_calls:,}  ->  unique variants: {n_unique:,}")
    return truth, inp


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--maf",
        default=str(
            Path.home()
            / "git/datahub/public/msk_impact_2017/data_mutations.txt"
        ),
    )
    ap.add_argument("--nrows", type=int, default=None, help="limit MAF rows (dev)")
    ap.add_argument("--outdir", default=str(TRUTH_DIR), help="output dir for truth/input parquet")
    ap.add_argument("--sample", type=int, default=None, help="subsample to N unique variants (seed=0)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    truth, inp = build(args.maf, args.nrows)
    if args.sample and len(truth) > args.sample:
        idx = truth.sample(n=args.sample, random_state=0).index
        truth = truth.loc[idx].reset_index(drop=True)
        inp = inp.loc[idx].reset_index(drop=True)
        print(f"subsampled to {len(truth):,} unique variants")
    truth.to_parquet(outdir / "truth.parquet", index=False)
    inp.to_parquet(outdir / "input.parquet", index=False)
    print(f"wrote {outdir/'truth.parquet'}")
    print(f"wrote {outdir/'input.parquet'}")


if __name__ == "__main__":
    main()
