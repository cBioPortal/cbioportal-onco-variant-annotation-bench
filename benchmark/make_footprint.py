"""Reduce a VCF and a MAF to one variant per 1 Mb window they touch.

The result carries the same genomic footprint as the full input with a tiny
fraction of the variants, so a tool that loads its transcript index lazily by
region still pays its whole index cost on it. Subtracting this point from the
full run leaves per-variant work and nothing else.
"""
import sys


def pick(lines, chrom_i, pos_i):
    seen, out = set(), []
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        try:
            w = (f[chrom_i], int(f[pos_i]) // 1_000_000)
        except (ValueError, IndexError):
            continue
        if w not in seen:
            seen.add(w)
            out.append(ln)
    return out


def main() -> None:
    vcf_in, vcf_out, maf_in, maf_out = sys.argv[1:5]

    src = open(vcf_in).readlines()
    body = [l for l in src if not l.startswith("#")]
    sub = pick(body, 0, 1)
    open(vcf_out, "w").writelines([l for l in src if l.startswith("#")] + sub)

    m = open(maf_in).readlines()
    msub = pick(m[1:], 1, 2)
    open(maf_out, "w").writelines([m[0]] + msub)

    print(f"  footprint: {len(sub)} of {len(body)} variants (vcf), "
          f"{len(msub)} of {len(m) - 1} (maf)")


if __name__ == "__main__":
    main()
