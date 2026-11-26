"""Expand vibe-vep's `all_effects` MAF column into an all_transcripts parquet.

`vibe-vep annotate maf` reports its picked transcript in the core MAF columns and
every overlapping transcript in `all_effects`, as `;`-separated entries of

    SYMBOL,<consequence terms>,HGVSp_short,ENST.version,HGVSc,IMPACT,CANONICAL

The consequence is itself a comma-separated list, so an entry has more than
seven fields whenever a transcript carries more than one SO term. The fields are
therefore located relative to the one holding the transcript id, not by a fixed
index -- reading `ENST` at a fixed position silently drops every multi-term row,
about 11% of them on MSK-IMPACT.

That column is what gives vibe-vep a Stage 1 (annotation) view alongside its
Stage 2 (pick) one, so it can be compared with VEP-family tools on the reference
transcript rather than only on the transcript it chose.

`all_effects` carries only the short protein form, so `hgvsp` is left empty and
Stage 1 protein concordance for vibe-vep is read from `hgvsp_short`. The MAF
class per transcript is derived from the consequence with the same vcf2maf
mapping the VEP-family parsers use.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (  # noqa: E402
    STD_FIELDS,
    STD_FIELDS_OPTIONAL,
    strip_transcript_version,
    var_key,
)
from so_to_maf import to_maf_class  # noqa: E402



def _vcf_shaped(ref: str, alt: str) -> tuple[str, str]:
    """TCGA-convention alleles -> left-anchored pair, for the MAF-class rules."""
    if ref == "-":
        return "N", "N" + alt
    if alt == "-":
        return "N" + ref, "N"
    return ref, alt


def build(maf_path: str) -> pd.DataFrame:
    maf = pd.read_csv(
        maf_path, sep="\t", comment="#", dtype=str, na_filter=False,
        encoding="latin-1", low_memory=False,
    )
    if "all_effects" not in maf.columns:
        raise SystemExit(f"{maf_path} has no all_effects column")

    rows: list[dict] = []
    for r in maf.itertuples(index=False):
        key = var_key(r.Chromosome, r.Start_Position, r.End_Position,
                      r.Reference_Allele, r.Tumor_Seq_Allele2)
        ref, alt = _vcf_shaped(r.Reference_Allele, r.Tumor_Seq_Allele2)
        for entry in str(r.all_effects).split(";"):
            if not entry:
                continue
            f = entry.split(",")
            # Anchor on the transcript id: everything between the symbol and the
            # protein change is the consequence term list.
            i = next((k for k, v in enumerate(f) if v.startswith("ENST")), None)
            if i is None or i < 2 or len(f) < i + 4:
                continue
            consequence = "&".join(t for t in f[1:i - 1] if t)
            hgvsp_short = unquote(f[i - 1])
            rows.append({
                "var_key": key,
                "hugo_symbol": f[0],
                "transcript_id": strip_transcript_version(f[i]),
                "transcript_id_versioned": f[i],
                "consequence": consequence,
                "variant_classification": to_maf_class(consequence, ref, alt),
                "hgvsc": unquote(f[i + 1]),
                "hgvsp": "",
                "hgvsp_short": hgvsp_short,
                "protein_position": "",
                "codons": "",
            })

    out = pd.DataFrame(rows)[STD_FIELDS + STD_FIELDS_OPTIONAL]
    return out.drop_duplicates(subset=["var_key", "transcript_id"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--maf", required=True, help="vibe-vep annotated MAF")
    ap.add_argument("--out", required=True, help="all_transcripts.parquet")
    args = ap.parse_args()
    out = build(args.maf)
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dest, index=False)
    n_var = out["var_key"].nunique()
    print(f"wrote {dest}  ({len(out):,} rows, {n_var:,} variants, "
          f"{len(out)/max(1,n_var):.1f} tx/variant)")


if __name__ == "__main__":
    main()
