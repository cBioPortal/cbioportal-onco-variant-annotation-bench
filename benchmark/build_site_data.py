"""Precompute the JSON the HTML explorer loads.

Outputs into website/data/:
  summary.json      leaderboard (both stages), tool metadata, headline counts
  by_gene.json      per-gene per-tool concordance (pick / protein / class)
  by_class.json     per variant-classification per-tool concordance
  variants.json     discrepancy browser: every variant where >=1 tool disagrees
                    with MSK on transcript or protein, with per-tool picks and the
                    full per-transcript detail for drill-down

Concordance is version-insensitive on transcripts and exact on HGVSp_Short.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from common import STD_FIELDS, strip_transcript_version  # noqa: F401

ROOT = Path(__file__).resolve().parent.parent
OVERRIDE_URL = "https://github.com/genome-nexus/genome-nexus-importer/blob/master/data/common_input/isoform_overrides_at_mskcc_grch37.txt"

TRACKS = {
    "grch37": {
        "label": "GRCh37 · MSK-IMPACT",
        "genome": "GRCh37 / hg19",
        "truth_set": "MSK-IMPACT 2017",
        "truth": ROOT / "data" / "truth" / "truth.parquet",
        "tools_dir": ROOT / "data" / "tools",
        "results": ROOT / "results",
        "out": ROOT / "website" / "data" / "grch37",
        "allt": ROOT / "data" / "tools" / "fastvep" / "all_transcripts.parquet",
        "override_file": ROOT / "data" / "msk_overrides" / "isoform_overrides_at_mskcc_grch37.txt",
        "no_override_pick": 0.6823,
        "tools": [
            {"id": "genome_nexus", "label": "Genome Nexus (truth pipeline)", "light": "#2a78d6", "dark": "#3987e5",
             "ver": {"software": "Genome Nexus (MSK production, 2017)", "transcripts": "Ensembl GRCh37 + MSK isoform overrides"}},
            {"id": "gn_vep111", "label": "Genome Nexus (VEP 111 endpoint)", "light": "#00a3c4", "dark": "#33c4e6",
             "ver": {"software": "vep111.genomenexus.org (isoformOverrideSource=mskcc)", "transcripts": "Ensembl GRCh37 release 111 + MSK isoform overrides"}},
            {"id": "vep_vcf2maf", "label": "Ensembl VEP 111 + vcf2maf (reference)", "light": "#eb6834", "dark": "#d95926",
             "ver": {"software": "Ensembl VEP 111.0 + vcf2maf @589406f", "transcripts": "VEP indexed cache release 111 (GRCh37 gene set)"}},
            {"id": "vep_vcf2maf_116", "label": "Ensembl VEP 116 + vcf2maf (latest)", "light": "#c0399a", "dark": "#db5cb5",
             "ver": {"software": "Ensembl VEP 116.0 + vcf2maf @589406f", "transcripts": "VEP indexed cache release 116 (GRCh37 gene set)"}},
            {"id": "fastvep_mafsmith", "label": "fastVEP + mafsmith (+ MSK isoforms)", "light": "#008300", "dark": "#22a022",
             "ver": {"software": "fastVEP 0.2.0 @d0cd1f7 + mafsmith @e009e98", "transcripts": "Ensembl GRCh37 release 87 (GFF3)"}},
            {"id": "vibe_vep", "label": "vibe-vep (single Go binary)", "light": "#4a3aa7", "dark": "#9085e9",
             "ver": {"software": "vibe-vep @c3ca6d8", "transcripts": "GENCODE v19 (GRCh37)"}},
            {"id": "fastvep", "label": "fastVEP (--pick, no override option)", "light": "#1baf7a", "dark": "#199e70",
             "ver": {"software": "fastVEP 0.2.0 @d0cd1f7 (--pick)", "transcripts": "Ensembl GRCh37 release 87 (GFF3)"}},
            {"id": "opencravat", "label": "OpenCRAVAT (hg19→hg38 liftover)", "light": "#e34948", "dark": "#e66767",
             "ver": {"software": "OpenCRAVAT 3.1.1 (gencode mapper)", "transcripts": "GENCODE — hg19 input lifted to hg38 (7 variants unmapped)"}},
        ],
    },
    "grch38": {
        "label": "GRCh38 · TCGA-GDC",
        "genome": "GRCh38 / hg38",
        "truth_set": "TCGA-LUAD (GDC)",
        "truth": ROOT / "data" / "grch38" / "truth" / "truth.parquet",
        "tools_dir": ROOT / "data" / "grch38" / "tools",
        "results": ROOT / "results" / "grch38",
        "out": ROOT / "website" / "data" / "grch38",
        "allt": ROOT / "data" / "grch38" / "tools" / "fastvep" / "all_transcripts.parquet",
        "override_file": None,
        "no_override_pick": None,
        "tools": [
            {"id": "vep_vcf2maf", "label": "Ensembl VEP 111 + vcf2maf (reference)", "light": "#eb6834", "dark": "#d95926",
             "ver": {"software": "Ensembl VEP 111.0 + vcf2maf", "transcripts": "VEP indexed cache release 111 (GRCh38)"}},
            {"id": "vep_vcf2maf_116", "label": "Ensembl VEP 116 + vcf2maf (latest)", "light": "#c0399a", "dark": "#db5cb5",
             "ver": {"software": "Ensembl VEP 116.0 + vcf2maf", "transcripts": "VEP indexed cache release 116 (GRCh38)"}},
            {"id": "genome_nexus", "label": "Genome Nexus (GRCh38)", "light": "#2a78d6", "dark": "#3987e5",
             "ver": {"software": "Genome Nexus (grch38.genomenexus.org)", "transcripts": "Ensembl GRCh38 + GN overrides"}},
            {"id": "fastvep_mafsmith", "label": "fastVEP + mafsmith", "light": "#008300", "dark": "#22a022",
             "ver": {"software": "fastVEP 0.2.0 + mafsmith @e009e98", "transcripts": "Ensembl GRCh38 release 112 (GFF3)"}},
            {"id": "vibe_vep", "label": "vibe-vep (single Go binary)", "light": "#4a3aa7", "dark": "#9085e9",
             "ver": {"software": "vibe-vep @c3ca6d8", "transcripts": "GENCODE (GRCh38)"}},
            {"id": "opencravat", "label": "OpenCRAVAT", "light": "#e34948", "dark": "#e66767",
             "ver": {"software": "OpenCRAVAT 3.1.1 (gencode mapper)", "transcripts": "GENCODE (GRCh38, native)"}},
            {"id": "fastvep", "label": "fastVEP (--pick)", "light": "#1baf7a", "dark": "#199e70",
             "ver": {"software": "fastVEP 0.2.0 (--pick)", "transcripts": "Ensembl GRCh38 release 112 (GFF3)"}},
        ],
    },
    "grch37_50k": {
        "label": "GRCh37 · MSK-IMPACT 50K",
        "genome": "GRCh37 / hg19",
        "truth_set": "MSK-IMPACT 50K (2026)",
        "truth": ROOT / "data" / "50k" / "truth" / "truth.parquet",
        "tools_dir": ROOT / "data" / "50k" / "tools",
        "results": ROOT / "results" / "50k",
        "out": ROOT / "website" / "data" / "grch37_50k",
        "allt": ROOT / "data" / "50k" / "tools" / "fastvep" / "all_transcripts.parquet",
        "override_file": ROOT / "data" / "msk_overrides" / "isoform_overrides_at_mskcc_grch37.txt",
        "no_override_pick": None,
        "tools": [
            {"id": "genome_nexus", "label": "Genome Nexus (truth pipeline)", "light": "#2a78d6", "dark": "#3987e5",
             "ver": {"software": "Genome Nexus (MSK production, 2026)", "transcripts": "Ensembl GRCh37 + MSK isoform overrides"}},
            {"id": "gn_vep111", "label": "Genome Nexus (VEP 111 endpoint)", "light": "#00a3c4", "dark": "#33c4e6",
             "ver": {"software": "vep111.genomenexus.org (isoformOverrideSource=mskcc)", "transcripts": "Ensembl GRCh37 release 111 + MSK isoform overrides"}},
            {"id": "vep_vcf2maf", "label": "Ensembl VEP 116 + vcf2maf (reference)", "light": "#eb6834", "dark": "#d95926",
             "ver": {"software": "Ensembl VEP 116.0 + vcf2maf @589406f", "transcripts": "VEP indexed cache release 116 (GRCh37 gene set)"}},
            {"id": "fastvep_mafsmith", "label": "fastVEP + mafsmith (+ MSK isoforms)", "light": "#008300", "dark": "#22a022",
             "ver": {"software": "fastVEP 0.2.0 @d0cd1f7 + mafsmith", "transcripts": "Ensembl GRCh37 release 87 (GFF3)"}},
            {"id": "vibe_vep", "label": "vibe-vep (single Go binary)", "light": "#4a3aa7", "dark": "#9085e9",
             "ver": {"software": "vibe-vep @c3ca6d8", "transcripts": "GENCODE v19 (GRCh37)"}},
            {"id": "fastvep", "label": "fastVEP (--pick, no override option)", "light": "#1baf7a", "dark": "#199e70",
             "ver": {"software": "fastVEP 0.2.0 @d0cd1f7 (--pick)", "transcripts": "Ensembl GRCh37 release 87 (GFF3)"}},
        ],
    },
}

# These module globals are set per-track at the top of build().
TRUTH = TOOLS_DIR = RESULTS = OUT = ALLT = OVERRIDE_FILE = None
TOOLS = []
GENOME = TRUTH_SET = ""
NO_OVERRIDE_PICK = None


def norm(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def load_tool(tid: str) -> pd.DataFrame:
    a = pd.read_parquet(TOOLS_DIR / tid / "annotations.parquet")
    a["tx"] = a["transcript_id"].map(strip_transcript_version)
    return a


def build(track: str) -> None:
    global TRUTH, TOOLS_DIR, RESULTS, OUT, ALLT, OVERRIDE_FILE, TOOLS, GENOME, TRUTH_SET, NO_OVERRIDE_PICK
    cfg = TRACKS[track]
    TRUTH, TOOLS_DIR, RESULTS, OUT = cfg["truth"], cfg["tools_dir"], cfg["results"], cfg["out"]
    ALLT, OVERRIDE_FILE, TOOLS = cfg["allt"], cfg["override_file"], cfg["tools"]
    GENOME, TRUTH_SET, NO_OVERRIDE_PICK = cfg["genome"], cfg["truth_set"], cfg["no_override_pick"]

    OUT.mkdir(parents=True, exist_ok=True)
    truth = pd.read_parquet(TRUTH)
    truth["tx"] = truth["transcript_id"].map(strip_transcript_version)
    tools = [t for t in TOOLS if (TOOLS_DIR / t["id"] / "annotations.parquet").exists()]

    # ---- Per-tool all-transcripts tables (drives Stage 1 coverage + drill-down) ----
    # tid -> DataFrame(var_key, tx, hgvsp_short, variant_classification, consequence, hugo_symbol)
    tool_allt: dict[str, pd.DataFrame] = {}
    # tid -> {var_key: set(tx)} for O(1) coverage lookups
    tool_tx_by_var: dict[str, dict] = {}
    for t in tools:
        p = TOOLS_DIR / t["id"] / "all_transcripts.parquet"
        if not p.exists():
            continue
        a = pd.read_parquet(p)
        a["tx"] = a["transcript_id"].map(strip_transcript_version)
        tool_allt[t["id"]] = a
        tool_tx_by_var[t["id"]] = {vk: s for vk, s in a.groupby("var_key")["tx"].agg(set).items()}
    # Per-variant union transcript set across every tool.
    union_tx_by_var: dict[str, set] = {}
    for tid, m in tool_tx_by_var.items():
        for vk, s in m.items():
            union_tx_by_var.setdefault(vk, set()).update(s)

    # Unified per-variant frame: truth + each tool's pick + match flags.
    base = truth[
        ["var_key", "hugo_symbol", "tx", "hgvsp_short", "hgvsp", "variant_classification", "consequence"]
    ].rename(
        columns={
            "hugo_symbol": "gene",
            "tx": "truth_tx",
            "hgvsp_short": "truth_hgvsp_short",
            "hgvsp": "truth_hgvsp",
            "variant_classification": "truth_class",
            "consequence": "truth_csq",
        }
    )
    df = base.copy()
    for t in tools:
        tid = t["id"]
        a = load_tool(tid)[["var_key", "tx", "hgvsp_short", "variant_classification"]].rename(
            columns={
                "tx": f"{tid}__tx",
                "hgvsp_short": f"{tid}__hgvsp",
                "variant_classification": f"{tid}__class",
            }
        )
        df = df.merge(a, on="var_key", how="left")
        df[f"{tid}__tx_match"] = norm(df[f"{tid}__tx"]) == norm(df["truth_tx"])
        df[f"{tid}__hgvsp_match"] = norm(df[f"{tid}__hgvsp"]) == norm(df["truth_hgvsp_short"])
        df[f"{tid}__class_match"] = norm(df[f"{tid}__class"]) == norm(df["truth_class"])

    # ---- OncoKB filter presets (variant-level filters across every view) ----
    onco_path = ROOT / "data" / "oncokb" / track / "oncokb.parquet"
    if onco_path.exists():
        ok = pd.read_parquet(onco_path)[
            ["var_key", "oncokb_is_oncogenic", "oncokb_oncogenic", "oncokb_highest_level"]
        ]
        df = df.merge(ok, on="var_key", how="left")
        df["onco"] = df["oncokb_is_oncogenic"].fillna(False).astype(bool)
        df["onco_label"] = df["oncokb_oncogenic"].fillna("")
        df["onco_level"] = df["oncokb_highest_level"].fillna("")
    else:
        df["onco"] = False
        df["onco_label"] = ""
        df["onco_level"] = ""

    # Each preset: (id, human label, boolean mask). Precomputed as <base>_<id>.json.
    FILTERS = [
        ("onco", "Oncogenic (any)", df["onco"]),
        ("level1", "Level 1 (FDA-recognized)", df["onco_level"] == "LEVEL_1"),
        ("level12", "Level 1–2", df["onco_level"].isin(["LEVEL_1", "LEVEL_2"])),
        ("level123", "Level 1–3", df["onco_level"].isin(["LEVEL_1", "LEVEL_2", "LEVEL_3A"])),
        ("resistance", "Resistance", df["onco_label"] == "Resistance"),
    ]
    FILTERS = [(fid, lbl, m) for fid, lbl, m in FILTERS if int(m.sum()) > 0]
    filter_keys = {fid: set(df.loc[m, "var_key"]) for fid, _, m in FILTERS}
    onco_keys = set(df.loc[df["onco"], "var_key"])
    n_oncogenic = int(df["onco"].sum())

    _S2MAP = {
        "annotated_rate": "annotated",
        "gene_concordance": "match_hugo_symbol",
        "protein_change_concordance": "match_hgvsp_short",
        "protein_change_concordance_full": "match_hgvsp",
        "transcript_concordance": "match_transcript",
        "variant_class_concordance": "match_variant_classification",
        "consequence_concordance": "match_consequence",
        "hgvsc_concordance": "match_hgvsc",
    }

    def stage2_subset(tid: str, keys: set) -> dict | None:
        """Recompute the Stage-2 metrics on a variant subset from compare_<tid>.parquet."""
        cp = RESULTS / f"compare_{tid}.parquet"
        if not cp.exists() or not keys:
            return None
        c = pd.read_parquet(cp)
        c = c[c["var_key"].isin(keys)]
        if not len(c):
            return None
        return {k: float(c[col].mean()) for k, col in _S2MAP.items() if col in c.columns}

    # ---- summary.json ----
    stage2 = [json.loads(p.read_text()) for p in sorted(RESULTS.glob("metrics_*.json"))]
    stage2 = {m["tool"]: m for m in stage2}
    stage1 = {
        json.loads(p.read_text())["tool"]: json.loads(p.read_text())
        for p in sorted(RESULTS.glob("transcripts_*.json"))
    }
    # MSK isoform-override list (GRCh37 track only, from genome-nexus-importer).
    override_rows = []
    override_genes: set[str] = set()
    if OVERRIDE_FILE and OVERRIDE_FILE.exists():
        ov = pd.read_csv(OVERRIDE_FILE, sep="\t", dtype=str).fillna("")
        for r in ov.itertuples(index=False):
            override_rows.append({"gene": r.gene_name, "refseq": r.refseq_id, "enst": r.enst_id})
            override_genes.add(r.gene_name)
    with_pick = stage2.get("fastvep_mafsmith", {}).get("transcript_concordance")

    # Variants where at least one tool disagrees with MSK on transcript or protein.
    disc_mask = pd.Series(False, index=df.index)
    for t in tools:
        tid = t["id"]
        disc_mask |= ~df[f"{tid}__tx_match"] | ~df[f"{tid}__hgvsp_match"]
    n_discrepant = int(disc_mask.sum())

    n_discrepant_onco = int((disc_mask & df["onco"]).sum())
    filters_meta = [
        {"id": fid, "label": lbl, "n": int(m.sum()),
         "n_discrepant": int((disc_mask & m).sum()),
         "n_genes": int(df.loc[m, "gene"].replace("", pd.NA).nunique())}
        for fid, lbl, m in FILTERS
    ]
    summary = {
        "track": track,
        "n_variants": int(len(truth)),
        "n_discrepant": n_discrepant,
        "n_oncogenic": n_oncogenic,
        "n_discrepant_onco": n_discrepant_onco,
        "filters": filters_meta,
        "truth_set": TRUTH_SET,
        "genome": GENOME,
        "n_genes": int(truth["hugo_symbol"].nunique()),
        "override": {
            "source": "genome-nexus-importer · isoform_overrides_at_mskcc_grch37",
            "url": OVERRIDE_URL,
            "n_genes": len(override_rows),
            "no_override_pick": NO_OVERRIDE_PICK,
            "with_override_pick": with_pick,
        } if override_rows else None,
        "tools": [
            {
                **t,
                "stage2": {
                    k: stage2[t["id"]].get(k)
                    for k in [
                        "annotated_rate",
                        "gene_concordance",
                        "transcript_concordance",
                        "protein_change_concordance_full",
                        "protein_change_concordance",
                        "variant_class_concordance",
                        "consequence_concordance",
                        "hgvsc_concordance",
                    ]
                }
                if t["id"] in stage2
                else {},
                "stage2_onco": stage2_subset(t["id"], onco_keys),
                "stage2_filters": {fid: stage2_subset(t["id"], filter_keys[fid]) for fid, _, _ in FILTERS},
            }
            for t in tools
        ],
        "stage1": list(stage1.values()),
    }
    (OUT / "summary.json").write_text(json.dumps(summary))
    override_rows.sort(key=lambda r: r["gene"])
    (OUT / "overrides.json").write_text(json.dumps({"rows": override_rows}))

    # ---- by_class.json (+ oncogenic-only variant) ----
    def compute_by_class(dsub: pd.DataFrame) -> list:
        out = []
        for cls, g in dsub.groupby(norm(dsub["truth_class"])):
            if not cls:
                continue
            row = {"variant_class": cls, "n": int(len(g))}
            for t in tools:
                tid = t["id"]
                row[tid] = {
                    "tx": float(g[f"{tid}__tx_match"].mean()),
                    "hgvsp": float(g[f"{tid}__hgvsp_match"].mean()),
                    "class": float(g[f"{tid}__class_match"].mean()),
                }
            out.append(row)
        out.sort(key=lambda r: r["n"], reverse=True)
        return out

    (OUT / "by_class.json").write_text(json.dumps(compute_by_class(df)))
    (OUT / "by_class_onco.json").write_text(json.dumps(compute_by_class(df[df["onco"]])))
    for fid, _, m in FILTERS:
        (OUT / f"by_class_{fid}.json").write_text(json.dumps(compute_by_class(df[m])))

    # ---- by_gene.json (+ oncogenic-only variant) ----
    def compute_by_gene(dsub: pd.DataFrame) -> list:
        out = []
        for gene, g in dsub.groupby("gene"):
            if not gene:
                continue
            row = {"gene": gene, "n": int(len(g))}
            worst = 1.0
            for t in tools:
                tid = t["id"]
                txm = float(g[f"{tid}__tx_match"].mean())
                hpm = float(g[f"{tid}__hgvsp_match"].mean())
                clm = float(g[f"{tid}__class_match"].mean())
                row[tid] = {"tx": txm, "hgvsp": hpm, "class": clm}
                worst = min(worst, hpm)
            row["worst_hgvsp"] = worst  # for sorting "where it goes wrong"
            row["override"] = gene in override_genes  # MSK pins a transcript for this gene
            out.append(row)
        out.sort(key=lambda r: (r["worst_hgvsp"], -r["n"]))
        return out

    (OUT / "by_gene.json").write_text(json.dumps(compute_by_gene(df)))
    (OUT / "by_gene_onco.json").write_text(json.dumps(compute_by_gene(df[df["onco"]])))
    for fid, _, m in FILTERS:
        (OUT / f"by_gene_{fid}.json").write_text(json.dumps(compute_by_gene(df[m])))

    # ---- coverage.json — per-tool coverage vs. union-of-all-tools transcripts ----
    #
    # For each variant, the "expected" transcript set is the UNION of what every
    # tool emits (no external ground truth). A tool's coverage on that variant is
    # |tool_tx ∩ union| / |union|. Aggregate = mean over variants where union > 0.
    # Also aggregated per gene (using the truth's Hugo_Symbol per variant), and we
    # record the count of (var, tx) rows each tool emits vs. omits, plus how often
    # the tool covers the MSK reference transcript (from the truth MAF, which
    # already reflects mskcc isoform overrides).
    gene_by_vk = dict(zip(df["var_key"], df["gene"]))
    truth_tx_by_vk = dict(zip(df["var_key"], df["truth_tx"]))
    # Which genes have an explicit MSK isoform override on file — subset metric.
    override_genes: set[str] = set()
    if OVERRIDE_FILE and OVERRIDE_FILE.exists():
        _ov = pd.read_csv(OVERRIDE_FILE, sep="\t", dtype=str).fillna("")
        override_genes = {r.gene_name for r in _ov.itertuples(index=False) if r.gene_name}
    coverage_tools = list(tool_tx_by_var.keys())
    n_var_with_union = 0
    tool_cov_sum = {tid: 0.0 for tid in coverage_tools}
    tool_emitted = {tid: 0 for tid in coverage_tools}      # total (var, tx) rows the tool emits
    tool_omitted = {tid: 0 for tid in coverage_tools}      # union rows the tool doesn't have
    tool_union_sum = {tid: 0 for tid in coverage_tools}    # sum of |union| across variants
    # MSK reference-transcript presence — per tool, on every variant AND on the
    # subset where the gene has an explicit isoform override.
    tool_msk_ref_hit = {tid: 0 for tid in coverage_tools}
    tool_msk_ref_denom = {tid: 0 for tid in coverage_tools}
    tool_msk_ovr_hit = {tid: 0 for tid in coverage_tools}
    tool_msk_ovr_denom = {tid: 0 for tid in coverage_tools}
    # "All transcripts present" — variants where the tool covers 100% of the union.
    tool_all_present_hit = {tid: 0 for tid in coverage_tools}
    per_gene = {}
    for vk, u in union_tx_by_var.items():
        n_u = len(u)
        if not n_u:
            continue
        n_var_with_union += 1
        gene = gene_by_vk.get(vk, "") or ""
        truth_tx = truth_tx_by_vk.get(vk, "") or ""
        is_override_gene = gene in override_genes
        g = per_gene.setdefault(gene, {"n_variants": 0, **{tid: {"cov_sum": 0.0, "emitted": 0, "omitted": 0, "union_sum": 0} for tid in coverage_tools}})
        g["n_variants"] += 1
        for tid in coverage_tools:
            tx = tool_tx_by_var[tid].get(vk, set())
            k = len(tx & u)
            emitted = len(tx)
            omitted = n_u - k
            frac = k / n_u
            tool_cov_sum[tid] += frac
            tool_emitted[tid] += emitted
            tool_omitted[tid] += omitted
            tool_union_sum[tid] += n_u
            if k == n_u:
                tool_all_present_hit[tid] += 1
            if truth_tx:
                tool_msk_ref_denom[tid] += 1
                if truth_tx in tx:
                    tool_msk_ref_hit[tid] += 1
                if is_override_gene:
                    tool_msk_ovr_denom[tid] += 1
                    if truth_tx in tx:
                        tool_msk_ovr_hit[tid] += 1
            g[tid]["cov_sum"] += frac
            g[tid]["emitted"] += emitted
            g[tid]["omitted"] += omitted
            g[tid]["union_sum"] += n_u
    coverage = {
        "tools": coverage_tools,
        "n_variants_with_union": n_var_with_union,
        "n_override_genes": len(override_genes),
        "aggregate": {
            tid: {
                "mean_coverage": tool_cov_sum[tid] / n_var_with_union if n_var_with_union else 0.0,
                "micro_coverage": (tool_union_sum[tid] - tool_omitted[tid]) / tool_union_sum[tid] if tool_union_sum[tid] else 0.0,
                "n_emitted": tool_emitted[tid],
                "n_omitted": tool_omitted[tid],
                "n_union": tool_union_sum[tid],
                "msk_ref_tx_present": (tool_msk_ref_hit[tid] / tool_msk_ref_denom[tid]) if tool_msk_ref_denom[tid] else None,
                "msk_ref_tx_denom": tool_msk_ref_denom[tid],
                "msk_override_tx_present": (tool_msk_ovr_hit[tid] / tool_msk_ovr_denom[tid]) if tool_msk_ovr_denom[tid] else None,
                "msk_override_tx_denom": tool_msk_ovr_denom[tid],
                "all_transcripts_present": tool_all_present_hit[tid] / n_var_with_union if n_var_with_union else 0.0,
                "all_transcripts_present_hit": tool_all_present_hit[tid],
            }
            for tid in coverage_tools
        },
    }
    (OUT / "coverage.json").write_text(json.dumps(coverage))

    # ---- coverage_variants.json — per-variant coverage detail ----
    #
    # One row per "interesting" variant: max per-tool coverage < 1 OR at least
    # one tool doesn't return the MSK reference transcript. Fully-covered
    # variants are dropped because every column would read 100% + ✓.
    # Compact array shape to keep the file small:
    #   [vk, gene, chrom, pos, union_size, ref_tx,
    #    per_tool_covered_count[len(tools)],
    #    per_tool_ref_present_flag[len(tools)] (0/1)]
    cov_variants = []
    for vk, u in union_tx_by_var.items():
        n_u = len(u)
        if not n_u:
            continue
        truth_tx = truth_tx_by_vk.get(vk, "") or ""
        cov_counts = []
        ref_flags = []
        any_gap = False
        for tid in coverage_tools:
            tx = tool_tx_by_var[tid].get(vk, set())
            k = len(tx & u)
            cov_counts.append(k)
            in_ref = int(bool(truth_tx and truth_tx in tx))
            ref_flags.append(in_ref)
            if k < n_u or (truth_tx and not in_ref):
                any_gap = True
        if not any_gap:
            continue
        chrom, start, _end, _ref, _alt = vk.split(":")
        cov_variants.append([
            vk, gene_by_vk.get(vk, "") or "", chrom, int(start), n_u, truth_tx,
            cov_counts, ref_flags,
        ])
    cov_variants.sort(key=lambda v: (v[1], v[3]))
    (OUT / "coverage_variants.json").write_text(json.dumps({
        "tools": coverage_tools,
        "n_total": n_var_with_union,
        "n_shown": len(cov_variants),
        "variants": cov_variants,
    }))

    # ---- coverage_download.tsv — one row per interesting variant, wide layout ----
    # Every column that lets a user reproduce the per-tool coverage gap:
    # variant coords, union transcript set, MSK reference tx, and per tool:
    # covered count, whether MSK ref is present, emitted transcripts (⊆ union),
    # and missing transcripts (union − emitted).
    header = ["Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Variant_Allele",
              "Hugo_Symbol", "Union_Transcript_Count",
              "Union_Transcripts", "MSK_Reference_Transcript"]
    for tid in coverage_tools:
        header += [f"{tid}__Covered_Count",
                   f"{tid}__MSK_Reference_Present",
                   f"{tid}__Emitted_Transcripts",
                   f"{tid}__Missing_Transcripts"]
    lines = ["\t".join(header)]
    for v in cov_variants:
        vk, gene, chrom, pos, n_u, ref_tx, _covs, _refs = v
        _c, start, end, ref_a, alt_a = vk.split(":")
        union_set = union_tx_by_var[vk]
        row = [chrom, str(start), end, ref_a, alt_a, gene, str(n_u),
               ";".join(sorted(union_set)), ref_tx or ""]
        for tid in coverage_tools:
            tool_set = tool_tx_by_var[tid].get(vk, set()) & union_set
            missing = union_set - tool_set
            row += [
                str(len(tool_set)),
                "1" if (ref_tx and ref_tx in tool_set) else "0",
                ";".join(sorted(tool_set)),
                ";".join(sorted(missing)),
            ]
        lines.append("\t".join(row))
    (OUT / "coverage_download.tsv").write_text("\n".join(lines) + "\n")

    # ---- variants.json (discrepancies + transcript drill-down) ----
    tool_ids = [t["id"] for t in tools]
    disc_mask = pd.Series(False, index=df.index)
    for tid in tool_ids:
        disc_mask |= ~df[f"{tid}__tx_match"] | ~df[f"{tid}__hgvsp_match"]
    disc = df[disc_mask].copy()

    # Union all-transcripts view for the per-variant drill-down. Every tool's
    # all_transcripts.parquet contributes rows; we take the first tool's annotation
    # per (var_key, tx) and record `cov` = every tool that also emits that transcript
    # (so the front-end can highlight coverage gaps: cov ⊂ tool_ids => that tool
    # doesn't return this transcript for this variant).
    #
    # `interesting_keys` = Stage-2 discrepant OR Stage-1 coverage-interesting (union),
    # so the Transcripts (Stage 1) view's per-row download can lazy-load transcript
    # detail for any variant it displays.
    disc_keys = set(disc["var_key"])
    cov_interesting_keys = {v[0] for v in cov_variants}
    interesting_keys = disc_keys | cov_interesting_keys
    # (vk, tx) -> {"g": gene, "hp": ..., "cl": ..., "csq": ..., "cov": set()}
    tx_annot: dict[tuple[str, str], dict] = {}
    for tid, a in tool_allt.items():
        sub = a[a["var_key"].isin(interesting_keys)]
        for r in sub.itertuples(index=False):
            key = (r.var_key, r.tx)
            entry = tx_annot.get(key)
            if entry is None:
                tx_annot[key] = {
                    "g": r.hugo_symbol or "",
                    "hp": r.hgvsp_short or "",
                    "cl": r.variant_classification or "",
                    "csq": r.consequence or "",
                    "cov": [tid],
                }
            else:
                entry["cov"].append(tid)
                # Backfill missing annotation fields from later tools.
                if not entry["hp"] and r.hgvsp_short:
                    entry["hp"] = r.hgvsp_short
                if not entry["cl"] and r.variant_classification:
                    entry["cl"] = r.variant_classification
                if not entry["csq"] and r.consequence:
                    entry["csq"] = r.consequence
    # var_key -> list of (tx, annot dict) for the drill-down loop below.
    allt_by_var: dict = {}
    for (vk, tx), ent in tx_annot.items():
        allt_by_var.setdefault(vk, []).append((tx, ent))
    # Per-variant, per-tool picked transcript (to mark rows in the drill-down).
    pick_tx = {tid: dict(zip(df["var_key"], df[f"{tid}__tx"])) for tid in tool_ids}

    # Compact discrepancy table (variants.json) — no transcript detail, so it
    # stays small enough to load and filter client-side. Per-tool pick is a fixed
    # array [tx, hgvsp, class, tx_ok(0/1), hgvsp_ok(0/1)] in tool_ids order.
    variants = []
    # Per-gene transcript drill-down files (data/tx/<GENE>.json), loaded lazily.
    tx_dir = OUT / "tx"
    tx_dir.mkdir(exist_ok=True)
    tx_by_gene: dict[str, dict] = {}

    for r in disc.itertuples(index=False):
        vk = r.var_key
        chrom, start, end, ref, alt = vk.split(":")
        picks = [
            [
                getattr(r, f"{tid}__tx") or "",
                getattr(r, f"{tid}__hgvsp") or "",
                getattr(r, f"{tid}__class") or "",
                int(bool(getattr(r, f"{tid}__tx_match"))),
                int(bool(getattr(r, f"{tid}__hgvsp_match"))),
            ]
            for tid in tool_ids
        ]
        variants.append(
            {
                "vk": vk,
                "c": chrom,
                "p": int(start),
                "g": r.gene,
                "cl": r.truth_class,
                "tt": r.truth_tx,
                "th": r.truth_hgvsp_short,
                "pk": picks,
            }
        )
        # transcript detail — union of every tool's transcripts for this variant
        tx_rows = []
        gsub = allt_by_var.get(vk) or []
        seen = set()
        for tx, ent in gsub:
            tx_rows.append(
                {
                    "tx": tx,
                    "g": ent["g"],
                    "hp": ent["hp"],
                    "cl": ent["cl"],
                    "csq": ent["csq"],
                    "by": [tid for tid in tool_ids if pick_tx[tid].get(vk) == tx],
                    "cov": ent["cov"],
                    "t": int(tx == r.truth_tx),
                }
            )
            seen.add(tx)
        if r.truth_tx and r.truth_tx not in seen:
            tx_rows.append(
                {
                    "tx": r.truth_tx,
                    "g": r.gene,
                    "hp": r.truth_hgvsp_short,
                    "cl": r.truth_class,
                    "csq": r.truth_csq,
                    "by": [tid for tid in tool_ids if pick_tx[tid].get(vk) == r.truth_tx],
                    "cov": [],
                    "t": 1,
                }
            )
            seen.add(r.truth_tx)
        # Ensure each tool's picked transcript appears (some tools only emit their
        # pick, so they contribute nothing to the union above).
        for tid in tool_ids:
            ptx = pick_tx[tid].get(vk)
            if ptx and ptx not in seen:
                tx_rows.append(
                    {
                        "tx": ptx,
                        "g": r.gene,
                        "hp": getattr(r, f"{tid}__hgvsp") or "",
                        "cl": getattr(r, f"{tid}__class") or "",
                        "csq": "",
                        "by": [x for x in tool_ids if pick_tx[x].get(vk) == ptx],
                        "cov": [tid],
                        "t": 0,
                    }
                )
                seen.add(ptx)
        tx_by_gene.setdefault(r.gene or "_", {})[vk] = tx_rows

    # Also populate tx_by_gene for coverage-interesting variants that aren't in
    # the Stage-2 discrepancy set (so the Transcripts (Stage 1) per-row download
    # has data for every row it shows). No truth-vs-tool __hgvsp / __class delta
    # to record here, so we just build the union + per-tool pick markers.
    cov_only_keys = cov_interesting_keys - disc_keys
    if cov_only_keys:
        cov_only_df = df[df["var_key"].isin(cov_only_keys)][["var_key", "gene", "truth_tx", "truth_hgvsp_short", "truth_class", "truth_csq"]]
        for r in cov_only_df.itertuples(index=False):
            vk = r.var_key
            tx_rows = []
            gsub = allt_by_var.get(vk) or []
            seen = set()
            for tx, ent in gsub:
                tx_rows.append({
                    "tx": tx,
                    "g": ent["g"],
                    "hp": ent["hp"],
                    "cl": ent["cl"],
                    "csq": ent["csq"],
                    "by": [tid for tid in tool_ids if pick_tx[tid].get(vk) == tx],
                    "cov": ent["cov"],
                    "t": int(tx == r.truth_tx),
                })
                seen.add(tx)
            if r.truth_tx and r.truth_tx not in seen:
                tx_rows.append({
                    "tx": r.truth_tx,
                    "g": r.gene,
                    "hp": r.truth_hgvsp_short,
                    "cl": r.truth_class,
                    "csq": r.truth_csq,
                    "by": [tid for tid in tool_ids if pick_tx[tid].get(vk) == r.truth_tx],
                    "cov": [],
                    "t": 1,
                })
                seen.add(r.truth_tx)
            for tid in tool_ids:
                ptx = pick_tx[tid].get(vk)
                if ptx and ptx not in seen:
                    tx_rows.append({
                        "tx": ptx, "g": r.gene, "hp": "", "cl": "", "csq": "",
                        "by": [x for x in tool_ids if pick_tx[x].get(vk) == ptx],
                        "cov": [tid], "t": 0,
                    })
                    seen.add(ptx)
            tx_by_gene.setdefault(r.gene or "_", {})[vk] = tx_rows

    variants.sort(key=lambda v: (v["g"], v["p"]))
    payload = {
        "n_total": int(len(df)),
        "n_discrepant": int(len(variants)),
        "tool_ids": tool_ids,
        "variants": variants,
    }
    (OUT / "variants.json").write_text(json.dumps(payload))
    onco_variants = [v for v in variants if v["vk"] in onco_keys]
    (OUT / "variants_onco.json").write_text(json.dumps({
        "n_total": n_oncogenic,
        "n_discrepant": len(onco_variants),
        "tool_ids": tool_ids,
        "variants": onco_variants,
    }))
    for fid, _, m in FILTERS:
        fk = filter_keys[fid]
        fv = [v for v in variants if v["vk"] in fk]
        (OUT / f"variants_{fid}.json").write_text(json.dumps({
            "n_total": len(fk), "n_discrepant": len(fv), "tool_ids": tool_ids, "variants": fv,
        }))

    for gene, m in tx_by_gene.items():
        safe = gene.replace("/", "_")
        (tx_dir / f"{safe}.json").write_text(json.dumps(m))

    sizes = {p.name: f"{p.stat().st_size/1e6:.2f} MB" for p in OUT.glob("*.json")}
    tx_total = sum(p.stat().st_size for p in tx_dir.glob("*.json")) / 1e6
    print(f"[{track}] wrote {OUT}:")
    for k, v in sizes.items():
        print(f"  {k:18} {v}")
    print(f"  tx/ ({len(tx_by_gene)} gene files) {tx_total:.1f} MB total")
    print(f"discrepant variants: {len(variants):,} of {len(df):,}")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", default="grch37", choices=list(TRACKS))
    args = ap.parse_args()
    build(args.track)


if __name__ == "__main__":
    main()
