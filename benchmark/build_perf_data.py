"""Turn run_perf.sh's raw timings into website/data/<track>/perf.json.

Each tool is timed on a 10-variant input (startup: process launch + transcript
index load) and on the full track input. Annotation cost is the difference, so
a tool is credited for the work it does rather than for how its index loads.
Every figure is the median of the repetitions.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import statistics
import sys
from pathlib import Path

LABELS = {
    "vep_rs": "vep-rs 1.0.0",
    "vep_rs_pick": "vep-rs 1.0.0 (--pick)",
    "fastvep": "fastVEP 0.2.0 (--pick)",
    "fastvep_all": "fastVEP 0.2.0",
    "vibe_vep": "vibe-vep (main)",
    "vibe_vep_pr17": "vibe-vep (PR #17)",
    "vibe_vep_annotated": "vibe-vep (main) + SIFT/PolyPhen/ClinVar/AlphaMissense",
}

# What each tool was asked to produce, so the table is read with the right
# expectations: a tool that writes more fields is doing more work. `MODE` groups
# the rows that are genuinely comparable -- "all" computes and emits every
# overlapping transcript, "pick" reduces to one consequence per variant.
OUTPUT = {
    "vep_rs": "VCF, every transcript",
    "vep_rs_pick": "VCF, one pick per variant",
    "fastvep": "JSON, one pick per variant",
    "fastvep_all": "JSON, every transcript",
    "vibe_vep": "MAF, core columns replaced",
    "vibe_vep_pr17": "MAF, core columns replaced",
    "vibe_vep_annotated": "MAF + 18 annotation columns",
}

MODE = {
    "vep_rs": "all",
    "vep_rs_pick": "pick",
    "fastvep": "pick",
    "fastvep_all": "all",
    "vibe_vep": "all",
    "vibe_vep_pr17": "all",
    "vibe_vep_annotated": "all",
}


def cpu_model() -> str:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or "unknown"


def total_ram_gb() -> float:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemTotal"):
                return round(int(line.split()[1]) / 1048576, 1)
    except OSError:
        pass
    return 0.0


def med(vals: list[float]) -> float | None:
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", required=True, help="website/data/<track>/perf.json")
    args = ap.parse_args()

    by: dict[tuple[str, str], dict[str, list[float]]] = {}
    n_variants = 0
    reps: set[str] = set()
    with open(args.tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r.get("wall_s"):
                continue
            k = (r["tool"], r["phase"])
            d = by.setdefault(k, {"wall": [], "cpu": [], "rss": []})
            d["wall"].append(float(r["wall_s"]))
            d["cpu"].append(float(r["user_s"]) + float(r["sys_s"]))
            d["rss"].append(float(r["maxrss_kb"]) / 1024)
            n_variants = int(r["n_variants"])
            reps.add(r["rep"])

    tools = sorted({t for t, _ in by}, key=lambda t: list(LABELS).index(t) if t in LABELS else 99)
    rows = []
    for t in tools:
        s, f = by.get((t, "startup")), by.get((t, "full"))
        if not f:
            continue
        sw, sc = (med(s["wall"]), med(s["cpu"])) if s else (None, None)
        fw, fc = med(f["wall"]), med(f["cpu"])
        aw = fw - sw if sw is not None else None
        ac = fc - sc if sc is not None else None
        # A full run cheaper than a 10-variant run means the two were not doing
        # the same work -- a stale index rebuilt on one of them, say. Refuse to
        # publish a throughput derived from it.
        if (aw is not None and aw <= 0) or (ac is not None and ac <= 0):
            print(f"  !! {t}: full run cheaper than startup "
                  f"(wall {fw:.2f} vs {sw:.2f}) -- throughput suppressed", file=sys.stderr)
        rows.append({
            "tool": t,
            "label": LABELS.get(t, t),
            "output": OUTPUT.get(t, ""),
            "mode": MODE.get(t, ""),
            "startup_wall_s": sw,
            "startup_cpu_s": sc,
            "full_wall_s": fw,
            "full_cpu_s": fc,
            "annot_wall_s": aw,
            "annot_cpu_s": ac,
            "maxrss_mb": med(f["rss"]),
            "variants_per_s": (n_variants / aw) if aw and aw > 0 else None,
            "cpu_us_per_variant": (ac * 1e6 / n_variants) if ac and ac > 0 else None,
            "wall_runs": sorted(f["wall"]),
        })

    payload = {
        "n_variants": n_variants,
        "reps": len(reps),
        "host": {
            "cpu": cpu_model(),
            "cores_total": os.cpu_count(),
            "cores_used": 2,
            "ram_gb": total_ram_gb(),
        },
        "rows": rows,
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload))
    print(f"wrote {out} ({len(rows)} tools, {len(reps)} reps, {n_variants:,} variants)")


if __name__ == "__main__":
    main()
