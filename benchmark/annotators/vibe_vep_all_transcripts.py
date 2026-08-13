"""Parse vibe-vep's `annotate vcf` output into an all_transcripts parquet.

vibe-vep writes VEP-style CSQ, but with two format quirks vs. Ensembl VEP:

  1. Consequence terms are joined with ',' — VEP uses '&'. And vibe-vep also uses
     ',' as the CSQ *record* separator (matching VEP), so a naive `csq.split(',')`
     misaligns fields whenever a variant has multi-term consequences (very common,
     e.g. `missense_variant,splice_region_variant`).
  2. Each CSQ record has exactly 20 pipe-separated fields; the field order is
     stable per the ##INFO header:
       Allele | Consequence | IMPACT | SYMBOL | Gene | Feature_type | Feature |
       BIOTYPE | EXON | INTRON | HGVSc | HGVSp | cDNA_position | CDS_position |
       Protein_position | Amino_acids | Codons |
       CANONICAL_MSK | CANONICAL_ENSEMBL | CANONICAL_MANE
     So we split by '|' and chunk into 20-field records; the "junction" pipe-token
     between records N and N+1 holds "<last_field_of_N>,<allele_of_N+1>", split
     at the last comma.

The VCF ID column already carries our var_key (chr:start:end:ref:alt), so records
join straight back to truth with no coordinate re-derivation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import STD_FIELDS, hgvsp_to_short, strip_transcript_version  # noqa: E402

# 0-indexed positions in vibe-vep's CSQ format.
FIELDS_PER_RECORD = 20
F_ALLELE, F_CSQ, F_IMPACT, F_SYMBOL, F_GENE = 0, 1, 2, 3, 4
F_FTYPE, F_FEATURE, F_BIOTYPE, F_EXON, F_INTRON = 5, 6, 7, 8, 9
F_HGVSC, F_HGVSP, F_CDNA_POS, F_CDS_POS, F_PPOS = 10, 11, 12, 13, 14
F_AMINO, F_CODONS, F_CAN_MSK, F_CAN_ENS, F_CAN_MANE = 15, 16, 17, 18, 19


def parse_csq(csq_str: str) -> list[list[str]]:
    """Split vibe-vep's CSQ into 20-field records, unambiguously.

    Character-level state machine that treats ',' as a record separator only
    when we're at the LAST field of a record (position 19 = CANONICAL_MANE) —
    and as content anywhere else (only the Consequence field, position 1,
    actually contains commas in practice). '|' always flushes the current
    field. Returns a list of 20-string lists.
    """
    records: list[list[str]] = []
    current: list[str] = []
    buf = ""
    for c in csq_str:
        if c == "|":
            current.append(buf)
            buf = ""
        elif c == ",":
            # len(current) counts fields already flushed; buf holds the field
            # currently being accumulated. A comma at position (FIELDS_PER_RECORD-1)
            # means we're building the last field of the record — this is the
            # record boundary. Anywhere else, treat as field content.
            if len(current) == FIELDS_PER_RECORD - 1:
                current.append(buf)
                buf = ""
                records.append(current)
                current = []
            else:
                buf += c
        else:
            buf += c
    # Flush the final record (which doesn't end with a comma).
    if buf or current:
        current.append(buf)
        if len(current) == FIELDS_PER_RECORD:
            records.append(current)
    return records


def _strip_prefix(raw: str) -> str:
    """VEP-style 'ENST...:c.1A>T' -> 'c.1A>T'; leave raw if no prefix."""
    if not raw:
        return ""
    return raw.split(":", 1)[1] if ":" in raw else raw


def build(vcf_path: str) -> pd.DataFrame:
    rows: list[dict] = []
    with open(vcf_path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            var_key = cols[2]
            info = cols[7]
            csq = ""
            for part in info.split(";"):
                if part.startswith("CSQ="):
                    csq = part[4:]
                    break
            if not csq:
                continue
            for f in parse_csq(csq):
                if f[F_FTYPE] != "Transcript":
                    continue
                tx = f[F_FEATURE]
                if not tx.startswith("ENST"):
                    continue
                hgvsp = _strip_prefix(f[F_HGVSP])
                rows.append(
                    {
                        "var_key": var_key,
                        "hugo_symbol": f[F_SYMBOL],
                        "transcript_id": strip_transcript_version(tx),
                        "consequence": f[F_CSQ],
                        "variant_classification": "",
                        "hgvsc": _strip_prefix(f[F_HGVSC]),
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
    ap.add_argument("--vcf", required=True, help="vibe-vep annotate vcf output")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = build(args.vcf)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out, index=False)
    print(
        f"wrote {args.out}  ({len(out):,} transcript rows, "
        f"{out['var_key'].nunique():,} variants, "
        f"{len(out)/max(1,out['var_key'].nunique()):.1f} tx/variant)"
    )


if __name__ == "__main__":
    main()
