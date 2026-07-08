"""Stage 2 of the benchmark: transcript selection + effect, vs MSK's pick.

Each tool reports ONE transcript per variant (its pick) and that transcript's
effect (HGVSp, Variant_Classification, consequence). This compares that pick to
MSK's, so it measures the SELECTION policy on top of the annotator. Stage 1
(compare_transcripts.py) isolates the annotator itself.

Tool-agnostic: takes any parquet in STD_FIELDS format and joins on var_key.
Emits per-variant match flags and aggregate concordance metrics, including
breakdowns by variant classification (Benchmark 5) and gene (Benchmark 3/6
inputs). Writes a per-variant results parquet and a JSON metrics summary.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import strip_transcript_version

ROOT = Path(__file__).resolve().parent.parent
TRUTH = ROOT / "data" / "truth" / "truth.parquet"
RESULTS_DIR = ROOT / "results"


def norm(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def norm_hgvsc(s: pd.Series) -> pd.Series:
    """Strip any 'ENST...:' transcript prefix so 'ENST..:c.1A>T' == 'c.1A>T'."""
    return norm(s).str.replace(r"^[^:]*:", "", regex=True)


def compare(tool_name: str, tool_parquet: Path, truth_path: Path = TRUTH,
            results_dir: Path = RESULTS_DIR) -> None:
    truth = pd.read_parquet(truth_path).add_suffix("_truth").rename(
        columns={"var_key_truth": "var_key"}
    )
    tool = pd.read_parquet(tool_parquet).add_suffix("_tool").rename(
        columns={"var_key_tool": "var_key"}
    )
    m = truth.merge(tool, on="var_key", how="left", indicator=True)

    m["annotated"] = m["_merge"] == "both"
    # Treat an empty picked transcript as "not annotated" even if a row exists.
    m.loc[norm(m["transcript_id_tool"]) == "", "annotated"] = False

    for col in ("hugo_symbol", "hgvsp", "hgvsp_short", "variant_classification", "consequence"):
        m[f"match_{col}"] = norm(m[f"{col}_truth"]) == norm(m[f"{col}_tool"])
    # HGVSc: prefix-insensitive, and only scored where the truth has a value.
    m["match_hgvsc"] = norm_hgvsc(m["hgvsc_truth"]) == norm_hgvsc(m["hgvsc_tool"])
    m.loc[norm(m["hgvsc_truth"]) == "", "match_hgvsc"] = False

    m["tx_truth_nov"] = m["transcript_id_truth"].map(strip_transcript_version)
    m["tx_tool_nov"] = m["transcript_id_tool"].fillna("").map(strip_transcript_version)
    m["match_transcript"] = m["tx_truth_nov"] == m["tx_tool_nov"]

    # Unannotated variants can't match anything.
    match_cols = [c for c in m.columns if c.startswith("match_")]
    for c in match_cols:
        m.loc[~m["annotated"], c] = False

    n = len(m)
    metrics = {
        "tool": tool_name,
        "n_variants": int(n),
        "annotated_rate": float(m["annotated"].mean()),
        "gene_concordance": float(m["match_hugo_symbol"].mean()),
        "protein_change_concordance": float(m["match_hgvsp_short"].mean()),
        "protein_change_concordance_full": float(m["match_hgvsp"].mean()),
        "transcript_concordance": float(m["match_transcript"].mean()),
        "variant_class_concordance": float(m["match_variant_classification"].mean()),
        "consequence_concordance": float(m["match_consequence"].mean()),
        "hgvsc_concordance": float(m["match_hgvsc"].mean()),
    }

    # Benchmark 5: breakdown by variant classification (truth label).
    by_class = (
        m.groupby(norm(m["variant_classification_truth"]))
        .agg(
            n=("var_key", "size"),
            protein_change_concordance=("match_hgvsp_short", "mean"),
            transcript_concordance=("match_transcript", "mean"),
            annotated_rate=("annotated", "mean"),
        )
        .sort_values("n", ascending=False)
    )
    metrics["by_variant_class"] = {
        idx: {
            "n": int(r["n"]),
            "protein_change_concordance": float(r["protein_change_concordance"]),
            "transcript_concordance": float(r["transcript_concordance"]),
            "annotated_rate": float(r["annotated_rate"]),
        }
        for idx, r in by_class.iterrows()
    }

    results_dir.mkdir(parents=True, exist_ok=True)
    keep = [
        "var_key",
        "hugo_symbol_truth",
        "variant_classification_truth",
        "hgvsp_short_truth",
        "hgvsp_short_tool",
        "transcript_id_truth",
        "transcript_id_tool",
        "annotated",
    ] + match_cols
    m[keep].to_parquet(results_dir / f"compare_{tool_name}.parquet", index=False)
    with (results_dir / f"metrics_{tool_name}.json").open("w") as fh:
        json.dump(metrics, fh, indent=2)

    print(f"=== {tool_name}  ({n:,} variants) ===")
    print(f"  annotated:              {metrics['annotated_rate']:.4%}")
    print(f"  gene concordance:       {metrics['gene_concordance']:.4%}")
    print(f"  protein-change (short): {metrics['protein_change_concordance']:.4%}")
    print(f"  protein-change (full):  {metrics['protein_change_concordance_full']:.4%}")
    print(f"  transcript pick:        {metrics['transcript_concordance']:.4%}")
    print(f"  variant class:          {metrics['variant_class_concordance']:.4%}")
    print(f"  consequence:            {metrics['consequence_concordance']:.4%}")
    print(f"  coding change (HGVSc):  {metrics['hgvsc_concordance']:.4%}")
    print(f"  -> results/metrics_{tool_name}.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, help="tool name, e.g. genome_nexus")
    ap.add_argument(
        "--parquet",
        default=None,
        help="path to tool annotations parquet (default data/tools/<tool>/annotations.parquet)",
    )
    ap.add_argument("--truth", default=str(TRUTH), help="truth parquet")
    ap.add_argument("--results-dir", default=str(RESULTS_DIR), help="output results dir")
    args = ap.parse_args()
    parquet = Path(args.parquet) if args.parquet else (
        ROOT / "data" / "tools" / args.tool / "annotations.parquet"
    )
    compare(args.tool, parquet, Path(args.truth), Path(args.results_dir))


if __name__ == "__main__":
    main()
