"""Parse a mafsmith MAF (fastVEP + vcf2maf field mapping) into STD_FIELDS.

mafsmith (github.com/nf-osi/mafsmith) is a Rust drop-in for vcf2maf that
annotates with fastVEP and applies vcf2maf.pl's allele normalization, transcript
selection, and field mapping. So this "tool" is the full production flow
fastVEP+mafsmith, and its output is a standard MAF with the same columns as the
MSK truth — reduced here to one STD record per unique variant.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, load_maf, strip_transcript_version, var_key  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent.parent
DIR = ROOT / "data" / "tools" / "fastvep_mafsmith"


def picked_consequence(all_effects: str, symbol: str, hgvsp: str, vclass: str) -> str:
    """Recover the picked transcript's SO consequence from the all_effects field.

    all_effects is `SYMBOL,consequence,Variant_Classification,HGVSp;...` per
    transcript (no transcript id), so match on gene + full HGVSp, falling back to
    gene + Variant_Classification when HGVSp is empty (splice/intron/etc.).
    """
    if not all_effects:
        return ""
    fallback = ""
    for tok in all_effects.split(";"):
        parts = tok.split(",")
        if len(parts) < 4:
            continue
        sym, csq, cls, hp = parts[0], parts[1], parts[2], parts[3]
        if sym != symbol:
            continue
        if hgvsp and hp == hgvsp:
            return csq
        if not fallback and cls == vclass:
            fallback = csq
    return fallback


def url_decode(s: pd.Series) -> pd.Series:
    """MAF HGVS fields are URL-encoded per spec (e.g. p.S457= -> p.S457%3D)."""
    return s.map(lambda v: unquote(v) if isinstance(v, str) else v)


def build(maf_path: str) -> pd.DataFrame:
    maf = load_maf(maf_path)
    for col in ("HGVSp", "HGVSp_Short", "HGVSc"):
        if col in maf.columns:
            maf[col] = url_decode(maf[col])
    key = [
        var_key(c, s, e, r, a)
        for c, s, e, r, a in zip(
            maf["Chromosome"],
            maf["Start_Position"],
            maf["End_Position"],
            maf["Reference_Allele"],
            maf["Tumor_Seq_Allele2"],
        )
    ]
    out = pd.DataFrame(
        {
            "var_key": key,
            "hugo_symbol": maf.get("Hugo_Symbol", ""),
            "transcript_id": maf.get("Transcript_ID", "").map(strip_transcript_version),
            "consequence": [
                picked_consequence(ae, sym, hp, cls)
                for ae, sym, hp, cls in zip(
                    maf.get("all_effects", ""),
                    maf.get("Hugo_Symbol", ""),
                    maf.get("HGVSp", ""),
                    maf.get("Variant_Classification", ""),
                )
            ],
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
    ap.add_argument(
        "--out",
        default=None,
        help="output parquet (default: annotations.parquet beside the input MAF)",
    )
    args = ap.parse_args()
    out = build(args.maf)
    dest = Path(args.out) if args.out else Path(args.maf).parent / "annotations.parquet"
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dest, index=False)
    print(f"wrote {dest}  ({len(out):,} variants)")


if __name__ == "__main__":
    main()
