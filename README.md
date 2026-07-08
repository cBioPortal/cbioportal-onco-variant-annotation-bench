# cBioPortal Onco Variant Effect Benchmark

> **⚠️ Preliminary / work in progress.** Early, exploratory results — not a finished,
> validated, or peer-reviewed evaluation — and may contain errors.

**Why this exists.** cBioPortal harmonizes variant annotation across the many studies it
hosts using [Genome Nexus](https://www.genomenexus.org) — it re-annotates the *effect* of
each mutation (gene, protein change, transcript, consequence), but does **not** re-call
mutations (the variant calls come from each study's authors). As new annotation tools and
versions have appeared, we need a way to gauge how well they perform relative to the current
pipeline — that's what this benchmark is for.

A **preliminary, reproducible benchmark** of variant annotation / effect tools —
Genome Nexus, Ensembl VEP (+ vcf2maf), fastVEP (+ mafsmith), vibe-vep, and
OpenCRAVAT — scored against **published** cancer reference datasets.

**▶ Live interactive report: https://cbioportal.github.io/cbioportal-onco-variant-annotation-bench/**

There is no absolute ground truth for variant effect, so tools are scored against
the *published* annotation of each dataset (hence "Datahub" / "GDC", not "truth").
The report makes that explicit and lets you compare any two tools directly.

## Datasets (tracks)

| Track | Reference | Variants | Convention |
|-------|-----------|----------|------------|
| **GRCh37 · MSK-IMPACT (2017)** | cBioPortal Datahub (Genome Nexus) | 58,301 | MSK clinical isoform overrides |
| **GRCh37 · MSK-IMPACT 50K** | cBioPortal Datahub (newer/larger release) | 60,000 | MSK clinical isoform overrides |
| **GRCh38 · TCGA-GDC** | GDC VEP/GENCODE pipeline | 60,000 | GENCODE canonical (no overrides) |

Only **unique variants** are scored (recurrent variants counted once).

## Two-stage method

- **Stage 1 — Annotation:** for every transcript a tool reports, is the effect
  (HGVSp, class) right on the reference's transcript? Also measures completeness
  (transcripts / variant, and whether the reference transcript is covered). A
  property of the *annotator*.
- **Stage 2 — Selection:** which single transcript + effect does the tool pick,
  vs the reference's pick? A property of the *selection policy*.

Concordance is reported per metric (transcript pick, protein HGVSp, HGVSc, variant
class, consequence, gene) and broken down by variant class, by gene, and by
**OncoKB category / therapeutic level** (Oncogenic, Level 1–3, Resistance).

## What the report includes

- **Leaderboard** (Stage 1 + Stage 2), filterable by OncoKB preset
- **By variant class** and **by gene** drill-downs
- **Discrepancy browser** with per-variant, per-transcript detail
- **Compare tools** — direct two-tool agreement (no reference needed)
- **Isoforms** — per-gene canonical transcript from MANE / Ensembl / UniProt / MSKCC
- **2017 → 50K** — how the *same* published source re-annotated between releases
- **Download / LLM** — CSV exports, raw JSON links, a copy-paste LLM prompt, and
  in-page [WebMCP](https://github.com/webmachinelearning/webmcp) tools

## Reproduce

Python via [uv](https://github.com/astral-sh/uv). The pipeline is tool-agnostic —
each tool has a parser that maps its output to a normalized schema (`benchmark/common.py`).

```
benchmark/extract_input.py         # dataset MAF -> truth + input (unique variants)
benchmark/annotators/*.py          # run/parse each tool -> normalized annotations
benchmark/compare_transcripts.py   # Stage 1 (annotation)
benchmark/compare_pick.py          # Stage 2 (selection)
benchmark/build_site_data.py       # -> website/data/<track>/*.json (the report)
```

Adding a tool = one parser to the normalized schema, then re-run `build_site_data.py`.

## Contributing

Contributions and corrections are welcome via pull request — new tools, new
datasets, or fixes. The pipeline is tool-agnostic, so each addition is small and
self-contained.

- **Add a tool.** Write a parser in `benchmark/annotators/` that maps your tool's
  output to the normalized schema in `benchmark/common.py` (gene, transcript,
  HGVSp/HGVSc, variant class, consequence), add it to the tool list in
  `benchmark/build_site_data.py`, then run `compare_transcripts.py` (Stage 1),
  `compare_pick.py` (Stage 2), and `build_site_data.py`. It then appears in the
  leaderboard, the compare view, and the per-gene/per-class drill-downs.
- **Add a dataset.** Any MAF with the core columns works — `extract_input.py`
  builds the truth + input (unique variants), then add a track entry (reference
  label, genome build, isoform-override policy). Useful for testing a new panel or
  cohort.
- **Report a discrepancy.** If a "mismatch" is really a reference-data problem
  (e.g. a malformed value), open an issue — that feedback improves the underlying
  annotations too. (Example: [cBioPortal/datahub#2343](https://github.com/cBioPortal/datahub/issues/2343).)

Please keep raw variant data and any API keys out of commits (see the `.gitignore`
and `## Data & keys` below).

## Data & keys

This repo contains the **summarized report data** (`website/data/`) and the
pipeline code. It does **not** contain raw patient/variant data, tool caches, or
any API keys — the OncoKB token and all raw inputs live outside the repo. OncoKB
oncogenicity is used only to *group / filter* variants; see
[oncokb.org](https://www.oncokb.org) for its terms.
