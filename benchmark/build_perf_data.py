"""Turn run_perf.sh's raw timings into website/data/<track>/perf.json.

Each tool is timed at three input sizes:

  launch  ten variants in one or two places -- process start, nothing else
  index   one variant per 1 Mb window the full input touches -- everything the
          tool must load to cover the panel's footprint
  full    the whole track input -- index cost plus per-variant work

Subtracting `index` from `full` leaves per-variant work alone, which matters
because a tool that loads its transcript index lazily by region pays nothing
for regions it is not asked about. Timing such a tool on a handful of variants
understates its index cost by an order of magnitude, and that cost then lands
silently in its per-variant figure.

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

PHASES = ("launch", "index", "full")


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
    nv: dict[str, int] = {}
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
            nv[r["phase"]] = int(r["n_variants"])
            reps.add(r["rep"])

    n_full, n_index = nv.get("full", 0), nv.get("index", 0)
    span = n_full - n_index
    if span <= 0:
        sys.exit("full and index inputs are the same size; nothing to subtract")

    order = list(LABELS)
    tools = sorted({t for t, _ in by}, key=lambda t: order.index(t) if t in order else 99)
    rows = []
    for t in tools:
        if (t, "full") not in by:
            continue
        g = {p: by.get((t, p)) for p in PHASES}
        w = {p: (med(g[p]["wall"]) if g[p] else None) for p in PHASES}
        c = {p: (med(g[p]["cpu"]) if g[p] else None) for p in PHASES}

        per_wall = per_cpu = None
        if w["index"] is not None and w["full"] is not None:
            dw, dc = w["full"] - w["index"], c["full"] - c["index"]
            # A full run no dearer than the footprint run means the two were not
            # doing the same work -- a stale index rebuilt on one of them, say.
            # Refuse to publish a per-variant cost derived from it.
            if dw <= 0 or dc <= 0:
                print(f"  !! {t}: full run no dearer than the index run "
                      f"(wall {w['full']:.2f} vs {w['index']:.2f}) -- per-variant cost suppressed",
                      file=sys.stderr)
            else:
                per_wall, per_cpu = dw / span * 1e6, dc / span * 1e6

        rows.append({
            "tool": t,
            "label": LABELS.get(t, t),
            "output": OUTPUT.get(t, ""),
            "mode": MODE.get(t, ""),
            "launch_wall_s": w["launch"],
            "index_wall_s": w["index"],
            "index_cpu_s": c["index"],
            "full_wall_s": w["full"],
            "full_cpu_s": c["full"],
            "per_variant_wall_us": per_wall,
            "per_variant_cpu_us": per_cpu,
            "variants_per_s": (1e6 / per_wall) if per_wall else None,
            "maxrss_mb": med(g["full"]["rss"]),
            "wall_runs": sorted(g["full"]["wall"]),
        })

    payload = {
        "n_variants": n_full,
        "n_index_variants": n_index,
        "n_launch_variants": nv.get("launch", 0),
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
    print(f"wrote {out} ({len(rows)} tools, {len(reps)} reps, "
          f"index={n_index:,} full={n_full:,} variants)")


if __name__ == "__main__":
    main()
