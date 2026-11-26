"""Report Ensembl transcript subversions across tools.

The MAF truth sets carry bare transcript ids (ENST00000288602, no `.6`), so the
benchmark necessarily matches version-insensitively and no tool can be *scored*
on its subversion. That hides a real effect: the same ENST can have a different
CDS between Ensembl releases, so two tools built on different releases may
disagree on a variant's protein change while appearing to annotate "the same"
transcript — a transcript-model difference charged to the engine.

This script makes that visible. For each track it reports, per tool, how many
picked transcripts carry a version, and for every pair of tools how often they
disagree on the subversion of a shared transcript. A non-zero disagreement means
the two are on different transcript models and their Stage 1 figures are not a
like-for-like engine comparison.

Output: results/<track>/transcript_versions.json and a printed table.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import pandas as pd

from common import with_optional

ROOT = Path(__file__).resolve().parent.parent


def _versions(path: Path) -> dict[str, str]:
    """Map bare transcript id -> subversion, for the rows that carry one."""
    df = with_optional(pd.read_parquet(path))
    out: dict[str, str] = {}
    for tid in df["transcript_id_versioned"].astype(str):
        if tid.startswith("ENST") and "." in tid:
            bare, _, ver = tid.partition(".")
            out[bare] = ver
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tools-dir", required=True, help="data/.../tools")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument(
        "--parquet",
        default="all_transcripts.parquet",
        help="which per-tool parquet to read (all_transcripts.parquet or annotations.parquet)",
    )
    args = ap.parse_args()

    tools_dir = Path(args.tools_dir)
    per_tool: dict[str, dict[str, str]] = {}
    for d in sorted(tools_dir.iterdir()):
        f = d / args.parquet
        if d.is_dir() and f.exists():
            v = _versions(f)
            if v:
                per_tool[d.name] = v

    if not per_tool:
        raise SystemExit(f"no tool under {tools_dir} reports a versioned transcript id")

    report = {
        "parquet": args.parquet,
        "tools": {t: {"transcripts_with_version": len(v)} for t, v in per_tool.items()},
        "pairs": [],
    }

    print(f"=== transcript subversions ({args.parquet}) ===")
    for tool, v in per_tool.items():
        print(f"  {tool:16} {len(v):,} distinct transcripts carry a version")

    print("  -- pairwise subversion agreement on shared transcripts --")
    for a, b in itertools.combinations(per_tool, 2):
        va, vb = per_tool[a], per_tool[b]
        shared = set(va) & set(vb)
        diff = sorted(t for t in shared if va[t] != vb[t])
        rate = len(diff) / len(shared) if shared else 0.0
        report["pairs"].append({
            "a": a, "b": b, "shared": len(shared), "disagree": len(diff),
            "disagree_rate": rate,
            "examples": [{"transcript": t, a: va[t], b: vb[t]} for t in diff[:10]],
        })
        flag = "" if not diff else "   <-- different transcript models"
        print(f"  {a} vs {b}: {len(shared):,} shared, "
              f"{len(diff):,} disagree ({rate:.3%}){flag}")
        for t in diff[:5]:
            print(f"      {t}: {a} .{va[t]} vs {b} .{vb[t]}")

    results = Path(args.results_dir)
    results.mkdir(parents=True, exist_ok=True)
    dest = results / "transcript_versions.json"
    dest.write_text(json.dumps(report, indent=2))
    print(f"  -> {dest}")


if __name__ == "__main__":
    main()
