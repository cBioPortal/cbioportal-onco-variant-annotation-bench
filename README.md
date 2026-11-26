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
Genome Nexus, Ensembl VEP (+ vcf2maf), fastVEP (+ mafsmith), vibe-vep, vep-rs, and
OpenCRAVAT — scored against **published** cancer reference datasets.

**▶ Live interactive report: https://cbioportal.github.io/cbioportal-onco-variant-annotation-bench/**

There is no absolute ground truth for variant effect, so tools are scored against
the *published* annotation of each dataset (hence "Datahub" / "GDC", not "truth").
The report makes that explicit and lets you compare any two tools directly.

## Which tools this run covers

The committed report is a **re-run from regenerated inputs**, not a continuation of
an earlier one. Every tool in it was run against the same freshly built truth sets,
on one pinned transcript model (GENCODE v19 for GRCh37, v45 for GRCh38), so the
Stage 1 figures compare engines rather than transcript-database vintage --
`benchmark/compare_versions.py` reports 0 subversion disagreements across every
pair of tools.

**Ensembl VEP (111/116) + vcf2maf, fastVEP + mafsmith, and OpenCRAVAT are not in
this run.** Their inputs were not regenerated (the VEP caches alone exceed the disk
available on the machine this was rebuilt on), so `build_site_data.py` drops them
rather than carry forward numbers scored against a different truth snapshot.
Earlier published figures for them are in the git history.

Two rows need reading with care:

- **fastVEP is 0.2.0**, from the `maf-variant-classification` branch of a fork,
  chosen because it emits MAF `Variant_Classification`. Upstream fastVEP is well
  ahead of that version, so this row is not a fair reading of current fastVEP.
- **vep-rs + MSK isoform overrides** is not a released tool. It is vep-rs's own
  per-transcript output re-picked by `benchmark/annotators/vep_rs_msk.py` against
  the MSK override list -- a demonstration that the selection policy, not the
  engine, accounts for vep-rs's Stage 2 gap.

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
- **Speed & memory** (GRCh37 track) — wall time, CPU time and peak RSS per tool,
  measured separately from accuracy; see below

## Reproduce

Python via [uv](https://github.com/astral-sh/uv). The pipeline is tool-agnostic —
each tool has a parser that maps its output to a normalized schema (`benchmark/common.py`).

```
benchmark/extract_input.py         # dataset MAF -> truth + input (unique variants)
benchmark/maf_to_vcf.py            # -> VCF for VCF-based tools (VEP, vep-rs)
benchmark/input_to_maf.py          # -> minimal MAF for MAF-based tools (vibe-vep)
benchmark/annotators/*.py          # run/parse each tool -> normalized annotations
benchmark/so_to_maf.py             # SO consequence -> MAF class (vcf2maf rules)
benchmark/compare_transcripts.py   # Stage 1 (annotation)
benchmark/compare_pick.py          # Stage 2 (selection)
benchmark/leaderboard.py           # both stages as one table (--results-dir per track)
benchmark/build_site_data.py       # -> website/data/<track>/*.json (the report)
```

### Speed & memory

Timings are a **separate measurement from the leaderboard** and do not feed into
any concordance number. `benchmark/run_perf.sh` re-runs each locally-installed
tool (point `BENCH_DATA_DIR` at the directory holding `refs/`, `bin/`,
`vep_rs_cache/` and the per-track inputs) and `benchmark/build_perf_data.py`
reduces the raw TSV to `website/data/<track>/perf.json`.

Method, and its limits:

- Each tool is timed twice per repetition — on a 10-variant input (**startup**:
  process launch plus transcript-index load) and on the full track input.
  Annotation cost is the difference, so a tool is credited for the work it does
  rather than for how its index happens to load.
- The two vibe-vep builds are given separate `HOME`s. They share one
  `~/.vibe-vep` transcript cache and each invalidates the other's, so whichever
  ran second was being charged a 12 s rebuild that has nothing to do with it.
  Each tool also gets a discarded warm-up run before its repetitions.
- All tools are pinned to the same two cores (`taskset`) at the same degree of
  parallelism (`--fork 2` / `GOMAXPROCS=2`), run one at a time, under a 4 GB
  `systemd-run` scope. Median of three runs.
- Every run is a **cold process against a warm OS page cache**. Dropping the page
  cache needs root, so true cold-disk numbers are not claimed.
- **Genome Nexus is excluded.** It is a remote production API; timing it would
  measure the network and someone else's server.
- **vep-rs and fastVEP are each run in both modes** — every overlapping
  transcript, and `--pick` — and the report groups the rows by mode, because a
  tool emitting one consequence per variant is not doing the same work as one
  emitting all ~7.25.
- **vep-rs loads its cache lazily, per chromosome**, so its startup column is
  small and most of its index cost falls in the annotate column. The other tools
  load a whole transcript index up front.
- vibe-vep is listed three times: `main` and PR #17 with extra annotations off,
  which is the configuration this benchmark scores, and `main` in its default
  configuration, which also attaches SIFT/PolyPhen/ClinVar/AlphaMissense/gnomAD/
  dbSNP/hotspots. Core annotations are identical between them, so the gap is the
  cost of the extras.

### Transcript subversions

The reference MAFs carry bare transcript ids (`ENST00000288602`, no `.6`), so
matching is version-insensitive and no tool can be *scored* on its subversion.
Parsers still keep the tool's own versioned id in `transcript_id_versioned`;
`benchmark/compare_versions.py` reports, per pair of tools, how often they
disagree on the subversion of a shared transcript. A non-zero rate means the two
are on different transcript models, so a Stage 1 difference between them is a
database-vintage effect rather than an engine one. The report shows the
subversion wherever a tool's picked transcript is displayed.

```
benchmark/compare_versions.py --tools-dir data/tools --results-dir results
```

### vep-rs

[vep-rs](https://github.com/natera-open-source/vep-rs) is a Rust reimplementation of
Ensembl VEP. It annotates from a JSON transcript cache, which `vep-cache-builder`
builds from a GTF and a genome FASTA — so it can be pinned to the *same* GENCODE
release vibe-vep uses, which makes the two a comparison of engines rather than of
transcript-database vintage:

```
vep-cache-builder --species homo_sapiens --assembly GRCh37 --release 116 \
  --gtf gencode.v19.gtf.gz --genome-fasta GRCh37.fa --output-dir <cache>
vep --json_cache <cache> --assembly GRCh37 --offline --fasta GRCh37.fa \
  --vcf --symbol --hgvs --numbers --canonical --flag_pick_allele \
  -i input.vcf -o output.vep.vcf
benchmark/annotators/vep_rs.py --vep-vcf output.vep.vcf --outdir <tool dir>
```

One run yields both stages: `--flag_pick_allele` keeps every overlapping transcript
(Stage 1) and flags one as the pick (Stage 2).

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
