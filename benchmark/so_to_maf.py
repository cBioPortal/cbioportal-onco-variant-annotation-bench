"""Sequence Ontology consequence -> MAF Variant_Classification, vcf2maf-style.

VEP-family tools (Ensembl VEP, fastVEP, vep-rs) emit SO terms; cBioPortal MAFs
carry MAF variant classes. The mapping below follows vcf2maf's
`GetVariantClassification`, including its effect-priority ordering for the
`&`-joined consequence lists VEP writes.
"""
from __future__ import annotations

# vcf2maf's effect priority: lower number wins when a row lists several terms.
EFFECT_PRIORITY = {
    "transcript_ablation": 1,
    "exon_loss_variant": 1,
    "splice_donor_variant": 2,
    "splice_acceptor_variant": 2,
    "stop_gained": 3,
    "frameshift_variant": 3,
    "stop_lost": 3,
    "start_lost": 4,
    "initiator_codon_variant": 4,
    "missense_variant": 6,
    "protein_altering_variant": 6,
    "conservative_missense_variant": 6,
    "rare_amino_acid_variant": 6,
    "transcript_amplification": 7,
    "splice_region_variant": 8,
    "start_retained_variant": 9,
    "stop_retained_variant": 9,
    "synonymous_variant": 9,
    "incomplete_terminal_codon_variant": 10,
    "coding_sequence_variant": 11,
    "mature_miRNA_variant": 11,
    "exon_variant": 11,
    "5_prime_UTR_variant": 12,
    "5_prime_UTR_premature_start_codon_gain_variant": 12,
    "3_prime_UTR_variant": 12,
    "non_coding_exon_variant": 13,
    "non_coding_transcript_exon_variant": 13,
    "non_coding_transcript_variant": 14,
    "nc_transcript_variant": 14,
    "intron_variant": 14,
    "INTRAGENIC": 14,
    "intragenic_variant": 14,
    "NMD_transcript_variant": 15,
    "upstream_gene_variant": 16,
    "downstream_gene_variant": 16,
    "TFBS_ablation": 17,
    "TFBS_amplification": 17,
    "TF_binding_site_variant": 17,
    "regulatory_region_ablation": 17,
    "regulatory_region_amplification": 17,
    "regulatory_region_variant": 17,
    "regulatory_region": 17,
    "feature_elongation": 18,
    "feature_truncation": 18,
    "intergenic_variant": 19,
    "intergenic_region": 19,
}

_SPLICE_SITE = {
    "splice_acceptor_variant", "splice_donor_variant",
    "transcript_ablation", "exon_loss_variant",
}
_MISSENSE = {
    "missense_variant", "coding_sequence_variant",
    "conservative_missense_variant", "rare_amino_acid_variant",
}
_INTRON = {
    "transcript_amplification", "intron_variant", "INTRAGENIC", "intragenic_variant",
}
_SILENT = {
    "incomplete_terminal_codon_variant", "synonymous_variant",
    "stop_retained_variant", "NMD_transcript_variant", "start_retained_variant",
}
_RNA = {
    "mature_miRNA_variant", "exon_variant", "non_coding_exon_variant",
    "non_coding_transcript_exon_variant", "non_coding_transcript_variant",
    "nc_transcript_variant",
}
_IGR = {
    "TF_binding_site_variant", "regulatory_region_variant", "regulatory_region",
    "intergenic_variant", "intergenic_region",
}


def top_consequence(consequence: str) -> str:
    """Highest-priority SO term of an `&`-joined VEP consequence list.

    vcf2maf sorts the terms by priority with Perl's stable `sort`, so terms
    sharing a priority keep the order VEP wrote them in. `stop_gained` and
    `frameshift_variant` are both priority 3, and VEP emits
    `stop_gained&frameshift_variant`, which vcf2maf therefore scores as
    Nonsense_Mutation; ranking ties any other way flips those rows.
    """
    if not consequence:
        return ""
    terms = [t for t in consequence.split("&") if t]
    if not terms:
        return ""
    return min(
        enumerate(terms), key=lambda it: (EFFECT_PRIORITY.get(it[1], 20), it[0])
    )[1]


def variant_type(ref: str, alt: str) -> str:
    """MAF Variant_Type from VCF-shaped alleles (SNP/DNP/TNP/ONP/INS/DEL)."""
    lr, la = len(ref), len(alt)
    if lr == la:
        return {1: "SNP", 2: "DNP", 3: "TNP"}.get(lr, "ONP")
    return "DEL" if lr > la else "INS"


def to_maf_class(consequence: str, ref: str, alt: str, biotype: str = "") -> str:
    """MAF Variant_Classification for one transcript row (vcf2maf rules)."""
    effect = top_consequence(consequence)
    vt = variant_type(ref, alt)
    inframe = abs(len(ref) - len(alt)) % 3 == 0

    if effect in _SPLICE_SITE:
        return "Splice_Site"
    if effect == "stop_gained":
        return "Nonsense_Mutation"
    if (effect == "frameshift_variant"
            or (effect == "protein_altering_variant" and not inframe)):
        if vt == "DEL":
            return "Frame_Shift_Del"
        if vt == "INS":
            return "Frame_Shift_Ins"
    if effect == "stop_lost":
        return "Nonstop_Mutation"
    if effect in ("initiator_codon_variant", "start_lost"):
        return "Translation_Start_Site"
    if effect == "inframe_insertion" or (
            effect == "protein_altering_variant" and inframe and vt == "INS"):
        return "In_Frame_Ins"
    if effect == "inframe_deletion" or (
            effect == "protein_altering_variant" and inframe and vt == "DEL"):
        return "In_Frame_Del"
    if effect in _MISSENSE:
        return "Missense_Mutation"
    if effect in _INTRON:
        return "Intron"
    if effect == "splice_region_variant":
        return "Splice_Region"
    if effect in _SILENT:
        return "Silent"
    if effect in _RNA:
        return "RNA"
    if effect in ("5_prime_UTR_variant",
                  "5_prime_UTR_premature_start_codon_gain_variant"):
        return "5'UTR"
    if effect == "3_prime_UTR_variant":
        return "3'UTR"
    if effect in _IGR:
        return "IGR"
    if effect == "upstream_gene_variant":
        return "5'Flank"
    if effect == "downstream_gene_variant":
        return "3'Flank"
    # vcf2maf routes a non-coding biotype's leftovers to RNA before its fallback.
    if biotype and biotype not in ("protein_coding", "IG_V_gene", "TR_V_gene"):
        return "RNA"
    return "Targeted_Region"
