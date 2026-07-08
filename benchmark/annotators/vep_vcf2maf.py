"""Parse a VEP + vcf2maf MAF into STD_FIELDS.

Ensembl VEP (v112, GRCh37) driven by mskcc/vcf2maf.pl — the reference-standard
pipeline that fastVEP+mafsmith and vibe-vep reimplement. Fed the MSK isoform list
via vcf2maf's --custom-enst, so its transcript picks should match MSK closely.

vcf2maf writes a standard MAF with the core columns populated directly; HGVS
fields are URL-encoded per MAF spec and some annotation columns are latin-1.
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
DIR = ROOT / "data" / "tools" / "vep_vcf2maf"


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
    ap.add_argument("--out", default=str(DIR / "annotations.parquet"))
    args = ap.parse_args()
    out = build(args.maf)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_path, index=False)
    print(f"wrote {out_path}  ({len(out):,} variants)")


if __name__ == "__main__":
    main()
