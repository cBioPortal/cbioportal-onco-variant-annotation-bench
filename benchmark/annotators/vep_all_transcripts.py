"""Extract Ensembl VEP's per-transcript annotations into an all_transcripts parquet.

vcf2maf collapses VEP to a single picked transcript per variant (Stage 2), but the
VEP run itself annotates every overlapping transcript — those live in the CSQ INFO
field of the intermediate `*.vep.vcf` (VEP ran with --flag_pick_allele, so all
transcripts are kept, one flagged PICK). Parsing them gives the annotator-level
(Stage 1) view for VEP, the same as fastVEP's all_transcripts output.

The VCF ID column is the var_key (chr:start:end:ref:alt), preserved by VEP, so each
CSQ consequence group joins straight back to the input variant.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, hgvsp_to_short, strip_transcript_version  # noqa: E402

# 0-indexed CSQ subfield positions (from the ##INFO CSQ Format header).
F_CSQ, F_SYMBOL, F_FTYPE, F_FEATURE = 1, 3, 5, 6
F_HGVSC, F_HGVSP, F_PPOS, F_CODONS = 10, 11, 14, 16


def _hgvsp(raw: str) -> str:
    # VEP writes "ENSP00000123:p.Val600Glu"; keep the p. part.
    if not raw:
        return ""
    return raw.split(":", 1)[1] if ":" in raw else raw


def _hgvsc(raw: str) -> str:
    if not raw:
        return ""
    return raw.split(":", 1)[1] if ":" in raw else raw


def build(vep_vcf: str) -> pd.DataFrame:
    rows: list[dict] = []
    with open(vep_vcf, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            var_key = cols[2]
            info = cols[7]
            csq = ""
            for part in info.split(";"):
                if part.startswith("CSQ="):
                    csq = part[4:]
                    break
            if not csq:
                continue
            for grp in csq.split(","):
                f = grp.split("|")
                if len(f) <= F_CODONS or f[F_FTYPE] != "Transcript":
                    continue
                tx = f[F_FEATURE]
                if not tx.startswith("ENST"):
                    continue
                hgvsp = _hgvsp(f[F_HGVSP])
                rows.append(
                    {
                        "var_key": var_key,
                        "hugo_symbol": f[F_SYMBOL],
                        "transcript_id": strip_transcript_version(tx),
                        "consequence": f[F_CSQ],
                        "variant_classification": "",
                        "hgvsc": _hgvsc(f[F_HGVSC]),
                        "hgvsp": hgvsp,
                        "hgvsp_short": hgvsp_to_short(hgvsp) if hgvsp else "",
                        "protein_position": f[F_PPOS],
                        "codons": f[F_CODONS],
                    }
                )
    out = pd.DataFrame(rows)[STD_FIELDS]
    return out.drop_duplicates(subset=["var_key", "transcript_id"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep-vcf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = build(args.vep_vcf)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(
        f"wrote {args.out}  ({len(out):,} transcript rows, "
        f"{out['var_key'].nunique():,} variants, "
        f"{len(out)/max(1,out['var_key'].nunique()):.1f} tx/variant)"
    )


if __name__ == "__main__":
    main()
