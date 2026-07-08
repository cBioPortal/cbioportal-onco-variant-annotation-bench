"""Parse a vibe-vep annotated MAF into STD_FIELDS.

vibe-vep (github.com/inodb/vibe-vep) is a single-binary Go VEP/Genome-Nexus
reimplementation. Run with `annotate maf --replace --assembly GRCh37` to overwrite
the core columns (Hugo_Symbol, Consequence, Variant_Classification, Transcript_ID,
HGVSc, HGVSp, HGVSp_Short) with its own annotation of the picked transcript.

The MAF carries extra annotation columns (ClinVar etc.) with non-UTF-8 bytes, so
read latin-1; HGVS fields are URL-encoded per MAF convention, so URL-decode.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, strip_transcript_version, var_key  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
DIR = ROOT / "data" / "tools" / "vibe_vep"


def build(maf_path: str) -> pd.DataFrame:
    maf = pd.read_csv(
        maf_path, sep="\t", comment="#", dtype=str, na_filter=False,
        encoding="latin-1", low_memory=False,
    )
    for col in ("HGVSp", "HGVSp_Short", "HGVSc"):
        if col in maf.columns:
            maf[col] = maf[col].map(lambda v: unquote(v) if isinstance(v, str) else v)
    key = [
        var_key(c, s, e, r, a)
        for c, s, e, r, a in zip(
            maf["Chromosome"], maf["Start_Position"], maf["End_Position"],
            maf["Reference_Allele"], maf["Tumor_Seq_Allele2"],
        )
    ]
    out = pd.DataFrame(
        {
            "var_key": key,
            "hugo_symbol": maf.get("Hugo_Symbol", ""),
            "transcript_id": maf.get("Transcript_ID", "").map(strip_transcript_version),
            "consequence": maf.get("Consequence", ""),
            "variant_classification": maf.get("Variant_Classification", ""),
            "hgvsc": maf.get("HGVSc", ""),
            "hgvsp": maf.get("HGVSp", ""),
            "hgvsp_short": maf.get("HGVSp_Short", ""),
            "protein_position": maf.get("Protein_position", ""),
            "codons": maf.get("Codons", ""),
        }
    )
    return out.drop_duplicates(subset="var_key", keep="first").reset_index(drop=True)[
        STD_FIELDS
    ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--maf", default=str(DIR / "output.maf"))
    ap.add_argument("--out", default=None, help="output parquet (default beside the MAF)")
    args = ap.parse_args()
    out = build(args.maf)
    dest = Path(args.out) if args.out else Path(args.maf).parent / "annotations.parquet"
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dest, index=False)
    print(f"wrote {dest}  ({len(out):,} variants)")


if __name__ == "__main__":
    main()
