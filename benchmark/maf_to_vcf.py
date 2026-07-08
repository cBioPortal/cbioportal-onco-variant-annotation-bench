"""Convert the normalized input.parquet coordinates to a VCF for VCF-based tools.

MSK-IMPACT MAF uses the TCGA indel convention (`-` alleles); VCF requires a
left-anchoring reference base. We fetch that base from the GRCh37 FASTA:

  SNP / MNV (ref,alt same length, no `-`):
      POS=start, REF=ref, ALT=alt
  Deletion (alt == '-'):
      POS=start-1, REF=<base@start-1>+ref, ALT=<base@start-1>
  Insertion (ref == '-'):
      POS=start,   REF=<base@start>,       ALT=<base@start>+alt

The VCF ID column carries the var_key so tool output joins straight back to
truth. Records are emitted sorted by (chrom, pos) as VCF requires.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam

ROOT = Path(__file__).resolve().parent.parent
IN_PARQUET = ROOT / "data" / "truth" / "input.parquet"
FASTA = ROOT / "data" / "refs" / "Homo_sapiens.GRCh37.dna.primary_assembly.fa"

# Ensembl GRCh37 primary-assembly contig order.
CHROM_ORDER = {str(c): i for i, c in enumerate(list(range(1, 23)) + ["X", "Y", "MT"])}


def base_at(fa: pysam.FastaFile, chrom: str, pos1: int) -> str:
    """Return the reference base at 1-based position `pos1` (uppercase)."""
    return fa.fetch(chrom, pos1 - 1, pos1).upper()


def to_vcf_record(fa, chrom, start, end, ref, alt) -> tuple[int, str, str] | None:
    """Return (pos, vcf_ref, vcf_alt) or None if the variant can't be encoded."""
    start = int(start)
    if ref == "-":  # insertion
        anchor = base_at(fa, chrom, start)
        return start, anchor, anchor + alt
    if alt == "-":  # deletion
        anchor = base_at(fa, chrom, start - 1)
        return start - 1, anchor + ref, anchor
    # substitution (SNP / DNP / MNV): equal-length, already VCF-shaped
    return start, ref, alt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "data" / "tools" / "fastvep" / "input.vcf"))
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument(
        "--sample",
        default=None,
        help="emit a FORMAT/GT column for this tumor sample (needed by vcf2maf tools)",
    )
    ap.add_argument("--input", default=str(IN_PARQUET), help="input.parquet")
    ap.add_argument("--fasta", default=str(FASTA), help="reference FASTA (faidx'd)")
    args = ap.parse_args()

    inp = pd.read_parquet(args.input)
    if args.limit:
        inp = inp.head(args.limit)

    fa = pysam.FastaFile(str(args.fasta))
    contigs = set(fa.references)

    recs = []
    skipped = 0
    for r in inp.itertuples(index=False):
        chrom = str(r.chromosome)
        if chrom not in contigs:
            skipped += 1
            continue
        try:
            pos, vref, valt = to_vcf_record(
                fa, chrom, r.start, r.end, r.reference_allele, r.variant_allele
            )
        except Exception:  # noqa: BLE001
            skipped += 1
            continue
        recs.append((chrom, pos, r.var_key, vref, valt))

    recs.sort(key=lambda x: (CHROM_ORDER.get(x[0], 99), x[1]))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##reference=GRCh37\n")
        if args.sample:
            fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        for c in list(range(1, 23)) + ["X", "Y", "MT"]:
            fh.write(f"##contig=<ID={c}>\n")
        header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        if args.sample:
            header += f"\tFORMAT\t{args.sample}"
        fh.write(header + "\n")
        for chrom, pos, key, vref, valt in recs:
            row = f"{chrom}\t{pos}\t{key}\t{vref}\t{valt}\t.\t.\t."
            if args.sample:
                row += "\tGT\t0/1"
            fh.write(row + "\n")

    print(f"wrote {out}  ({len(recs):,} records, {skipped} skipped)")


if __name__ == "__main__":
    main()
