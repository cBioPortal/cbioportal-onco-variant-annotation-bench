"""Write the benchmark's unique-variant input back out as a minimal MAF.

MAF-driven annotators (vibe-vep) take a MAF rather than a VCF, so this is the
MAF counterpart of `maf_to_vcf.py`. The alleles stay in the TCGA convention the
input carries (`-` for the absent side of an indel).

`Tumor_Sample_Barcode` is a single synthetic sample: the benchmark scores unique
variants, so sample identity carries no information here.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

# `vibe-vep annotate maf --replace` overwrites the core columns in place and
# does not create ones the input lacks, so they are emitted empty here.
REPLACED_COLUMNS = [
    "Consequence", "Transcript_ID", "HGVSc", "HGVSp", "HGVSp_Short",
    "Protein_position", "Codons",
]

COLUMNS = [
    "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Strand",
    "Variant_Classification", "Variant_Type", "Reference_Allele",
    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode",
] + REPLACED_COLUMNS


def build(input_parquet: str, sample: str) -> pd.DataFrame:
    inp = pd.read_parquet(input_parquet)
    return pd.DataFrame({
        "Hugo_Symbol": "",
        "Chromosome": inp["chromosome"].astype(str),
        "Start_Position": inp["start"].astype(str),
        "End_Position": inp["end"].astype(str),
        "Strand": "+",
        "Variant_Classification": "",
        "Variant_Type": inp["variant_type"].astype(str),
        "Reference_Allele": inp["reference_allele"].astype(str),
        "Tumor_Seq_Allele1": inp["reference_allele"].astype(str),
        "Tumor_Seq_Allele2": inp["variant_allele"].astype(str),
        "Tumor_Sample_Barcode": sample,
        **{c: "" for c in REPLACED_COLUMNS},
    })[COLUMNS]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="input.parquet")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="BENCHMARK-SAMPLE")
    args = ap.parse_args()
    df = build(args.input, args.sample)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"wrote {out}  ({len(df):,} variants)")


if __name__ == "__main__":
    main()
