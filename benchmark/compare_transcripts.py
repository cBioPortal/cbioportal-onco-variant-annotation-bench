"""Stage 1 of the benchmark: transcript-level annotation, pick-independent.

The annotation and the transcript-selection questions are separated on purpose:

  Stage 1 (this script) — for every transcript a tool reports, is the annotation
  right? Measured against MSK's transcript regardless of what the tool would
  *pick*. Answers: does the tool cover MSK's transcript at all, and given that
  transcript does it produce the same HGVSp / consequence? This is a property of
  the ANNOTATOR (fastVEP, VEP, GN), shared by every picker built on top of it.

  Stage 2 (compare_pick.py) — which single transcript + effect does the tool
  report, vs MSK's pick? This is a property of the SELECTION policy (fastVEP
  --pick, vcf2maf canonical, vcf2maf + MSK isoform overrides).

Input: data/tools/<tool>/all_transcripts.parquet (long: one row per
var_key×transcript). Output: results/transcripts_<tool>.json + a per-transcript
match parquet.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import strip_transcript_version

ROOT = Path(__file__).resolve().parent.parent
TRUTH = ROOT / "data" / "truth" / "truth.parquet"
TOOLS = ROOT / "data" / "tools"
RESULTS = ROOT / "results"


def norm(s: pd.Series) -> pd.Series:
    return s.fillna("").astype(str).str.strip()


def compare(tool: str, parquet: Path, truth_path: Path = TRUTH,
            results_dir: Path = RESULTS) -> None:
    global RESULTS
    RESULTS = results_dir
    truth = pd.read_parquet(truth_path)
    truth["tx"] = truth["transcript_id"].map(strip_transcript_version)

    allt = pd.read_parquet(parquet)
    allt["tx"] = allt["transcript_id"].map(strip_transcript_version)

    # Coverage: is MSK's transcript present among the tool's transcripts?
    tx_sets = allt.groupby("var_key")["tx"].agg(set)
    cov = truth.merge(tx_sets.rename("tool_tx"), on="var_key", how="left")
    cov["tool_tx"] = cov["tool_tx"].apply(lambda x: x if isinstance(x, set) else set())
    cov["present"] = [tx in s for tx, s in zip(cov["tx"], cov["tool_tx"])]
    n_tx = allt.groupby("var_key")["tx"].nunique()

    # Per-transcript annotation concordance on MSK's transcript.
    j = allt.merge(
        truth[["var_key", "tx", "hgvsp", "hgvsp_short", "variant_classification", "consequence"]]
        .rename(
            columns={
                "hgvsp": "hgvsp_t",
                "hgvsp_short": "hgvsp_short_t",
                "variant_classification": "vc_t",
                "consequence": "csq_t",
            }
        ),
        on=["var_key", "tx"],
        how="inner",
    )
    prot = j[norm(j["hgvsp_t"]) != ""]

    metrics = {
        "tool": tool,
        "stage": "transcript",
        "n_truth_variants": int(len(truth)),
        "mean_transcripts_per_variant": float(n_tx.mean()),
        "msk_transcript_present_rate": float(cov["present"].mean()),
        "protein_full_concordance_on_msk_tx": float(
            (norm(prot["hgvsp"]) == norm(prot["hgvsp_t"])).mean()
        )
        if len(prot)
        else None,
        "protein_short_concordance_on_msk_tx": float(
            (norm(prot["hgvsp_short"]) == norm(prot["hgvsp_short_t"])).mean()
        )
        if len(prot)
        else None,
        "n_protein_variants": int(len(prot)),
    }
    if norm(j["variant_classification"]).ne("").any():
        cls = j[norm(j["vc_t"]) != ""]
        metrics["variant_class_concordance_on_msk_tx"] = float(
            (norm(cls["variant_classification"]) == norm(cls["vc_t"])).mean()
        )

    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"transcripts_{tool}.json").write_text(json.dumps(metrics, indent=2))

    print(f"=== {tool} — Stage 1 (transcript-level) ===")
    print(f"  mean transcripts / variant:      {metrics['mean_transcripts_per_variant']:.2f}")
    print(f"  MSK transcript present:          {metrics['msk_transcript_present_rate']:.4%}")
    if metrics["protein_full_concordance_on_msk_tx"] is not None:
        print(f"  protein (full) on MSK transcript:  {metrics['protein_full_concordance_on_msk_tx']:.4%}")
        print(f"  protein (short) on MSK transcript: {metrics['protein_short_concordance_on_msk_tx']:.4%}")
    if "variant_class_concordance_on_msk_tx" in metrics:
        print(f"  variant class on MSK transcript:   {metrics['variant_class_concordance_on_msk_tx']:.4%}")
    print(f"  -> results/transcripts_{tool}.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--parquet", default=None)
    ap.add_argument("--truth", default=str(TRUTH))
    ap.add_argument("--results-dir", default=str(RESULTS))
    args = ap.parse_args()
    parquet = Path(args.parquet) if args.parquet else (
        TOOLS / args.tool / "all_transcripts.parquet"
    )
    if not parquet.exists():
        raise SystemExit(
            f"no all_transcripts parquet for {args.tool} at {parquet} — "
            "Stage 1 needs the annotator's full transcript output"
        )
    compare(args.tool, parquet, Path(args.truth), Path(args.results_dir))


if __name__ == "__main__":
    main()
