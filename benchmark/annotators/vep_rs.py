"""Parse vep-rs VCF output into the benchmark's normalized schema.

vep-rs (natera-open-source/vep-rs) is a Rust reimplementation of Ensembl VEP.
Run it with `--vcf --symbol --hgvs --numbers --canonical --flag_pick_allele` so
every overlapping transcript is kept and one is flagged PICK, then this parser
emits both benchmark views from the single run:

  all_transcripts.parquet  one row per var_key x transcript  (Stage 1)
  annotations.parquet      the PICK row per variant          (Stage 2)

The CSQ subfield order is read from the `##INFO=<ID=CSQ...Format: ...>` header
rather than assumed, so a change in vep-rs's field set cannot silently shift the
columns. The VCF ID column carries the benchmark's var_key, so rows join back to
truth without re-deriving coordinates.

vep-rs emits SO terms; `so_to_maf` maps them to MAF Variant_Classification the
way vcf2maf does, which is what the MAF-derived truth carries.
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
    hgvsp_to_short,
    strip_transcript_version,
)
from so_to_maf import to_maf_class  # noqa: E402

CSQ_PREFIX = "##INFO=<ID=CSQ,"


def _csq_fields(path: str) -> list[str]:
    """Subfield names from the CSQ INFO header line."""
    with open(path, encoding="latin-1") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith(CSQ_PREFIX) and "Format:" in line:
                fmt = line.split("Format:", 1)[1].strip()
                return fmt.rstrip('">').strip().split("|")
    raise SystemExit(f"no CSQ Format header in {path}")


def _strip_prefix(raw: str) -> str:
    """'ENSP00000123.4:p.Val600Glu' -> 'p.Val600Glu', URL-decoded.

    VEP percent-encodes the characters VCF reserves inside an INFO field, so a
    synonymous change reaches this parser as `p.Pro617%3D`.
    """
    s = unquote(raw)
    return s.split(":", 1)[1] if ":" in s else s


def build(vep_vcf: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    names = _csq_fields(vep_vcf)
    idx = {n: i for i, n in enumerate(names)}
    need = ("Consequence", "SYMBOL", "Feature_type", "Feature", "HGVSc", "HGVSp",
            "Protein_position", "Codons", "BIOTYPE")
    missing = [n for n in need if n not in idx]
    if missing:
        raise SystemExit(f"CSQ header lacks {missing}; rerun vep with those flags")
    i_pick = idx.get("PICK")

    rows: list[dict] = []
    picks: list[dict] = []
    with open(vep_vcf, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            key, ref, alt, info = cols[2], cols[3], cols[4], cols[7]
            csq = next((p[4:] for p in info.split(";") if p.startswith("CSQ=")), "")
            if not csq:
                continue
            for grp in csq.split(","):
                f = grp.split("|")
                if len(f) < len(names) or f[idx["Feature_type"]] != "Transcript":
                    continue
                tx = f[idx["Feature"]]
                if not tx.startswith("ENST"):
                    continue
                hgvsp = _strip_prefix(f[idx["HGVSp"]])
                consequence = f[idx["Consequence"]]
                rec = {
                    "var_key": key,
                    "hugo_symbol": f[idx["SYMBOL"]],
                    "transcript_id": strip_transcript_version(tx),
                    "transcript_id_versioned": tx,
                    "consequence": consequence,
                    "variant_classification": to_maf_class(
                        consequence, ref, alt, f[idx["BIOTYPE"]]
                    ),
                    "hgvsc": _strip_prefix(f[idx["HGVSc"]]),
                    "hgvsp": hgvsp,
                    "hgvsp_short": hgvsp_to_short(hgvsp) if hgvsp else "",
                    "protein_position": f[idx["Protein_position"]],
                    "codons": f[idx["Codons"]],
                }
                rows.append(rec)
                if i_pick is not None and f[i_pick] == "1":
                    picks.append(rec)

    cols = STD_FIELDS + STD_FIELDS_OPTIONAL
    allt = pd.DataFrame(rows)[cols].drop_duplicates(
        subset=["var_key", "transcript_id"]
    ).reset_index(drop=True)
    pick = pd.DataFrame(picks)[cols].drop_duplicates(
        subset="var_key", keep="first"
    ).reset_index(drop=True)
    return allt, pick


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vep-vcf", required=True)
    ap.add_argument("--outdir", required=True, help="data/.../tools/vep_rs")
    args = ap.parse_args()
    allt, pick = build(args.vep_vcf)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    allt.to_parquet(out / "all_transcripts.parquet", index=False)
    pick.to_parquet(out / "annotations.parquet", index=False)
    n_var = allt["var_key"].nunique()
    print(f"wrote {out/'all_transcripts.parquet'}  ({len(allt):,} rows, "
          f"{n_var:,} variants, {len(allt)/max(1,n_var):.1f} tx/variant)")
    print(f"wrote {out/'annotations.parquet'}  ({len(pick):,} picked)")


if __name__ == "__main__":
    main()
