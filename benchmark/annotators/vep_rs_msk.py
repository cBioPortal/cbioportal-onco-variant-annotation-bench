"""vep-rs annotation + the MSK isoform-override selection policy.

vep-rs annotates every overlapping transcript but picks with Ensembl VEP's
`--pick` order, which knows nothing about MSK's clinical isoforms. That is the
whole of its Stage 2 gap: the reference transcript is already present in its
output, it is simply not the one chosen.

This re-picks from vep-rs's existing all_transcripts output using the policy
Genome Nexus applies -- if the gene has an MSK override transcript and that
transcript is among the ones reported for the variant, take it; otherwise keep
vep-rs's own pick. No re-annotation: the effect values are vep-rs's own.

It exists to separate the two questions the benchmark keeps apart. Stage 1 asks
whether the engine is right; this shows what the same engine scores once a
cBioPortal-shaped selection layer sits on top of it.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import (  # noqa: E402
    STD_FIELDS,
    STD_FIELDS_OPTIONAL,
    strip_transcript_version,
    with_optional,
)


def load_overrides(path: str) -> dict[str, str]:
    """gene symbol -> bare override transcript id."""
    ov = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return {
        r.gene_name: strip_transcript_version(r.enst_id)
        for r in ov.itertuples(index=False)
        if r.gene_name and r.enst_id
    }


def build(all_tx: str, pick: str, overrides: dict[str, str]) -> pd.DataFrame:
    allt = with_optional(pd.read_parquet(all_tx))
    base = with_optional(pd.read_parquet(pick)).set_index("var_key")

    # Rows whose transcript is its gene's MSK override.
    allt["is_override"] = [
        bool(g) and overrides.get(g, "") == t
        for g, t in zip(allt["hugo_symbol"], allt["transcript_id"])
    ]
    cand = allt[allt["is_override"]]

    # Where a variant has several override transcripts (overlapping genes),
    # prefer the one on the gene vep-rs itself picked, else the first reported.
    picked_gene = base["hugo_symbol"].to_dict()
    cand = cand.assign(
        _same_gene=[
            int(g == picked_gene.get(vk, ""))
            for vk, g in zip(cand["var_key"], cand["hugo_symbol"])
        ]
    ).sort_values(["var_key", "_same_gene"], ascending=[True, False])
    chosen = cand.drop_duplicates(subset="var_key", keep="first").set_index("var_key")

    cols = STD_FIELDS + STD_FIELDS_OPTIONAL
    out = base.copy()
    out.loc[chosen.index, [c for c in cols if c != "var_key"]] = chosen[
        [c for c in cols if c != "var_key"]
    ]
    out = out.reset_index()[cols]
    print(f"re-picked {len(chosen):,} of {len(out):,} variants onto an MSK override transcript")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--all-transcripts", required=True)
    ap.add_argument("--pick", required=True, help="vep-rs annotations.parquet")
    ap.add_argument("--overrides", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = build(args.all_transcripts, args.pick, load_overrides(args.overrides))
    dest = Path(args.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dest, index=False)
    print(f"wrote {dest}  ({len(out):,} variants)")


if __name__ == "__main__":
    main()
