"""Aggregate the two benchmark stages into one leaderboard.

Stage 1 (Annotation / transcript-level) from results/transcripts_<tool>.json —
written by compare_transcripts.py. A property of the annotator.

Stage 2 (Selection / pick) from results/metrics_<tool>.json — written by
compare_pick.py. A property of the transcript-selection policy.

Writes results/leaderboard.json and results/leaderboard.md.
"""
from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"


def pct(x) -> str:
    return "—" if x is None else f"{x*100:.1f}%"


def main() -> None:
    stage2 = [json.loads(p.read_text()) for p in sorted(RESULTS.glob("metrics_*.json"))]
    stage1 = {
        json.loads(p.read_text())["tool"]: json.loads(p.read_text())
        for p in sorted(RESULTS.glob("transcripts_*.json"))
    }

    n = stage2[0]["n_variants"] if stage2 else 0
    lines = [
        "# Annotation Benchmark Leaderboard",
        "",
        f"Truth set: MSK-IMPACT 2017 · {n:,} unique variants",
        "",
        "The benchmark is split into two stages that answer different questions.",
        "",
        "## Stage 1 — Annotation (transcript-level)",
        "",
        "For every transcript the tool reports, is the annotation right? Scored on "
        "MSK's transcript regardless of what the tool would *pick*. A property of "
        "the annotator, shared by every picker built on it.",
        "",
        "| Annotator | Transcripts / variant | MSK transcript present | Protein on MSK tx | Variant class on MSK tx |",
        "|-----------|-----------------------|------------------------|-------------------|--------------------------|",
    ]
    for tool, m in stage1.items():
        lines.append(
            f"| {tool} | {m['mean_transcripts_per_variant']:.1f} | "
            f"{pct(m['msk_transcript_present_rate'])} | "
            f"{pct(m['protein_full_concordance_on_msk_tx'])} | "
            f"{pct(m.get('variant_class_concordance_on_msk_tx'))} |"
        )

    lines += [
        "",
        "## Stage 2 — Selection (pick + effect)",
        "",
        "Which single transcript + effect does the tool report, vs MSK's pick? A "
        "property of the selection policy on top of the annotator.",
        "",
        "| Tool | Annotated | Gene | Transcript pick | Protein on pick (HGVSp) | Variant class | Consequence |",
        "|------|-----------|------|-----------------|-------------------------|---------------|-------------|",
    ]
    for m in stage2:
        lines.append(
            f"| {m['tool']} | {pct(m['annotated_rate'])} | {pct(m['gene_concordance'])} | "
            f"{pct(m['transcript_concordance'])} | "
            f"{pct(m.get('protein_change_concordance_full'))} | "
            f"{pct(m['variant_class_concordance'])} | {pct(m['consequence_concordance'])} |"
        )

    board = {"stage1": list(stage1.values()), "stage2": stage2}
    (RESULTS / "leaderboard.json").write_text(json.dumps(board, indent=2))
    (RESULTS / "leaderboard.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nwrote {RESULTS/'leaderboard.json'} and leaderboard.md")


if __name__ == "__main__":
    main()
