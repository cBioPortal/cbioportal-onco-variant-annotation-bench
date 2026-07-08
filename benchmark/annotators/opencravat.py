"""Parse OpenCRAVAT text-report output into STD_FIELDS.

OpenCRAVAT (KarchinLab) is a GENCODE-based annotator; run natively on GRCh38
(`oc run ... -l hg38 -m gencode -t text`). Its text report concatenates several
sections; we read only the variant-level section (one row per variant, primary
transcript), in input order — OpenCRAVAT renumbers UIDs but preserves order and
count, so we re-attach var_key from the input parquet by position.

Its protein change is three-letter (`p.Arg22Gly`) and its consequence vocabulary
is its own (`frameshift_truncation`, `splice_site_variant`), so hgvsp_short is
converted and consequence is left as OpenCRAVAT's SO term.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, hgvsp_to_short, strip_transcript_version  # noqa: E402

# 0-indexed columns in the variant-level section.
C_GENE, C_TX, C_SO, C_CDNA, C_PROT, C_ALLMAP, C_TAGS = 7, 8, 10, 12, 13, 14, 18


def row_keys(rows: list[list[str]], input_parquet: str) -> list[str]:
    """var_key per output row. Prefer the Tags column (survives liftover + dropped
    rows, e.g. hg19→hg38); fall back to input order when Tags is absent."""
    tags = [r[C_TAGS] if len(r) > C_TAGS else "" for r in rows]
    if any(":" in t for t in tags):
        return tags
    keys = pd.read_parquet(input_parquet)["var_key"].tolist()
    if len(rows) != len(keys):
        raise SystemExit(f"row mismatch: {len(rows)} OpenCRAVAT rows vs {len(keys)} inputs "
                         "(and no Tags column to join on)")
    return keys


def so_to_maf(terms: str) -> str:
    """Map a comma-joined SO string to a single MAF class (first mappable term)."""
    for t in terms.split(","):
        cls = SO_TO_MAF.get(t.strip())
        if cls:
            return cls
    return ""

# OpenCRAVAT's Sequence Ontology terms -> MAF Variant_Classification.
SO_TO_MAF = {
    "missense_variant": "Missense_Mutation",
    "synonymous_variant": "Silent",
    "stop_gained": "Nonsense_Mutation",
    "stop_lost": "Nonstop_Mutation",
    "start_lost": "Translation_Start_Site",
    "frameshift_truncation": "Frame_Shift_Del",
    "frameshift_elongation": "Frame_Shift_Ins",
    "inframe_deletion": "In_Frame_Del",
    "inframe_insertion": "In_Frame_Ins",
    "splice_site_variant": "Splice_Site",
    "complex_substitution": "Missense_Mutation",
    "2kb_upstream_variant": "5'Flank",
    "2kb_downstream_variant": "3'Flank",
    "5_prime_UTR_variant": "5'UTR",
    "3_prime_UTR_variant": "3'UTR",
    "intron_variant": "Intron",
    "unknown": "",
}


def read_variant_section(tsv: str) -> list[list[str]]:
    rows, in_var = [], False
    with open(tsv, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("UID\tChrom"):
                in_var = True
                continue
            if in_var:
                if line.startswith("#"):
                    break
                if line.strip():
                    rows.append(line.rstrip("\n").split("\t"))
    return rows


def build(tsv: str, input_parquet: str) -> pd.DataFrame:
    rows = read_variant_section(tsv)
    keys = row_keys(rows, input_parquet)
    out = pd.DataFrame(
        {
            "var_key": keys,
            "hugo_symbol": [r[C_GENE] if len(r) > C_GENE else "" for r in rows],
            "transcript_id": [strip_transcript_version(r[C_TX]) if len(r) > C_TX else "" for r in rows],
            "consequence": [r[C_SO] if len(r) > C_SO else "" for r in rows],
            "variant_classification": [
                SO_TO_MAF.get(r[C_SO], "") if len(r) > C_SO else "" for r in rows
            ],
            "hgvsc": [r[C_CDNA] if len(r) > C_CDNA else "" for r in rows],
            "hgvsp": "",
            "hgvsp_short": [hgvsp_to_short(r[C_PROT]) if len(r) > C_PROT else "" for r in rows],
            "protein_position": "",
            "codons": "",
        }
    )
    return out.drop_duplicates(subset="var_key", keep="first").reset_index(drop=True)[STD_FIELDS]


def build_all_transcripts(tsv: str, input_parquet: str) -> pd.DataFrame:
    """Expand the 'All Mappings' column into one row per var_key×transcript.

    Each mapping is 'ENST[.v]:GENE:UNIPROT:so,terms:p.Change:c.change', mappings
    separated by '; '. This is OpenCRAVAT's full per-transcript output — the
    Stage 1 (annotator) view, independent of which transcript it picks.
    """
    rows = read_variant_section(tsv)
    keys = row_keys(rows, input_parquet)
    out = []
    for key, r in zip(keys, rows):
        if len(r) <= C_ALLMAP:
            continue
        for mapping in r[C_ALLMAP].split(";"):
            mapping = mapping.strip()
            if not mapping:
                continue
            f = mapping.split(":")
            tx = strip_transcript_version(f[0]) if f else ""
            if not tx.startswith("ENST"):
                continue
            gene = f[1] if len(f) > 1 else ""
            so = f[3] if len(f) > 3 else ""
            prot = f[4] if len(f) > 4 else ""
            cdna = f[5] if len(f) > 5 else ""
            out.append(
                {
                    "var_key": key,
                    "hugo_symbol": gene,
                    "transcript_id": tx,
                    "consequence": so,
                    "variant_classification": so_to_maf(so),
                    "hgvsc": cdna,
                    "hgvsp": "",
                    "hgvsp_short": hgvsp_to_short(prot) if prot else "",
                    "protein_position": "",
                    "codons": "",
                }
            )
    df = pd.DataFrame(out)[STD_FIELDS]
    return df.drop_duplicates(subset=["var_key", "transcript_id"]).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--input", required=True, help="input.parquet used to run OpenCRAVAT (for order)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--all-out", default=None, help="also write all-transcripts parquet (Stage 1)")
    args = ap.parse_args()
    out = build(args.tsv, args.input)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(f"wrote {args.out}  ({len(out):,} variants)")
    if args.all_out:
        allt = build_all_transcripts(args.tsv, args.input)
        Path(args.all_out).parent.mkdir(parents=True, exist_ok=True)
        allt.to_parquet(args.all_out, index=False)
        print(f"wrote {args.all_out}  ({len(allt):,} transcript rows, "
              f"{allt['var_key'].nunique():,} variants, "
              f"{len(allt)/max(1,allt['var_key'].nunique()):.1f} tx/variant)")


if __name__ == "__main__":
    main()
