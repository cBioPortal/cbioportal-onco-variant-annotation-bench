"""Shared schema, MAF loading, and variant-key logic for the benchmark.

The benchmark is tool-agnostic: every annotator is converted to the same
`STD_FIELDS` record so `compare.py` never needs to know which tool produced it.
"""
from __future__ import annotations

import pandas as pd

# Columns pulled from the MSK-IMPACT MAF to form the truth set.
MAF_TRUTH_COLS = [
    "Hugo_Symbol",
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Strand",
    "Variant_Type",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "Consequence",
    "Variant_Classification",
    "HGVSc",
    "HGVSp",
    "HGVSp_Short",
    "Transcript_ID",
    "Protein_position",
    "Codons",
]

# Normalized "picked transcript" record every annotator (and the truth set) is
# reduced to. These are the fields the benchmark scores concordance on.
STD_FIELDS = [
    "var_key",
    "hugo_symbol",
    "transcript_id",
    "consequence",
    "variant_classification",
    "hgvsc",
    "hgvsp",
    "hgvsp_short",
    "protein_position",
    "codons",
]

# Columns a parser emits when its tool reports them. They are not part of the
# scored schema, so an older parser that omits one still loads; `with_optional`
# fills what is missing.
#
# `transcript_id_versioned` keeps the tool's transcript id as written, with its
# Ensembl subversion (ENST00000288602.6). Matching stays version-insensitive
# because the MAF truth sets carry no versions at all, but two tools built on
# different Ensembl releases can annotate the *same* ENST differently when its
# subversion differs, so the version is kept to tell that apart from a genuine
# disagreement between engines. `compare_versions.py` reports it.
STD_FIELDS_OPTIONAL = [
    "transcript_id_versioned",
]


def var_key(chrom, start, end, ref, alt) -> str:
    """Stable identifier for a genomic variant, shared across truth and tools."""
    chrom = str(chrom).replace("chr", "")
    return f"{chrom}:{start}:{end}:{ref}:{alt}"


def load_maf(path: str, nrows: int | None = None) -> pd.DataFrame:
    """Load a cBioPortal MAF, skipping the leading `#sequenced_samples` line."""
    return pd.read_csv(
        path,
        sep="\t",
        comment="#",
        dtype=str,
        nrows=nrows,
        na_filter=False,
        low_memory=False,
    )


def with_optional(df: "pd.DataFrame") -> "pd.DataFrame":
    """Add any missing optional column as an empty string column."""
    for col in STD_FIELDS_OPTIONAL:
        if col not in df.columns:
            df[col] = ""
    return df


def transcript_version(tid: str) -> str:
    """ENST00000288602.6 -> '6'; '' when the id carries no version."""
    if not tid or "." not in tid:
        return ""
    return tid.split(".", 1)[1]


def strip_transcript_version(tid: str) -> str:
    """ENST00000288602.6 -> ENST00000288602 for version-insensitive matching."""
    if not tid:
        return ""
    return tid.split(".")[0]


_AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U", "Pyl": "O", "Xaa": "X",
}


def strip_protein_prefix(hgvsp: str) -> str:
    """ENSP00000288602.6:p.Val600Glu -> p.Val600Glu (drop transcript/protein id)."""
    if not hgvsp:
        return ""
    return hgvsp.split(":")[-1]


def hgvsp_to_short(hgvsp: str) -> str:
    """Convert a 3-letter HGVS protein string to MSK short form.

    p.Val600Glu -> p.V600E, p.Arg273Ter -> p.R273*, p.Gly12= -> p.G12=.
    Longest-first token replacement handles fs/del/ins/dup/delins tails too.
    Best-effort: any 3-letter code found is swapped for its 1-letter code.
    """
    s = strip_protein_prefix(hgvsp)
    if not s:
        return ""
    for aa3, aa1 in _AA3TO1.items():
        s = s.replace(aa3, aa1)
    return s
