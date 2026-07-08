"""Compare two published MSK-IMPACT Datahub releases (2017 vs 50K) on shared variants.

The benchmark scores tools against a *published* annotation, not absolute truth —
and the published annotation itself changes between releases. This quantifies that:
for variants present in both the 2017 and 50K Datahub MAFs, how often did the gene,
transcript, protein change (HGVSp_Short) or variant classification change, and which
genes had their canonical transcript re-assigned.

Output: website/data/grch37/datahub_diff.json (summary + changed variants + the
per-gene transcript switches that drive most of the protein-change differences).
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parent.parent
VER_CACHE = ROOT / "data" / "transcript_ref" / "enst_versions_grch37.json"
HUB = Path(os.path.expanduser("~/git/datahub/public"))
M2017 = HUB / "msk_impact_2017" / "data_mutations.txt"
M50K = HUB / "msk_impact_50k_2026" / "data_mutations.txt"
OUT = ROOT / "website" / "data" / "grch37" / "datahub_diff.json"

COLS = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele",
        "Tumor_Seq_Allele2", "Hugo_Symbol", "HGVSp_Short", "Transcript_ID",
        "Variant_Classification"]


def load(f: Path) -> pd.DataFrame:
    d = pd.read_csv(f, sep="\t", comment="#", dtype=str, low_memory=False,
                    usecols=lambda c: c in COLS).fillna("")
    d["vk"] = (d["Chromosome"] + ":" + d["Start_Position"] + ":" + d["End_Position"]
               + ":" + d["Reference_Allele"] + ":" + d["Tumor_Seq_Allele2"])
    d["tx"] = d["Transcript_ID"].str.split(".").str[0]
    return d.drop_duplicates("vk").set_index("vk")


def fetch_versions(ids: set[str]) -> dict:
    """Current Ensembl GRCh37 .version for each ENST (the MSK MAFs carry none).
    Cached to disk; only missing ids are fetched from grch37.rest.ensembl.org."""
    cache = json.loads(VER_CACHE.read_text()) if VER_CACHE.exists() else {}
    todo = sorted(i for i in ids if i and i not in cache)
    for i in range(0, len(todo), 900):
        chunk = todo[i : i + 900]
        for attempt in range(3):
            try:
                r = requests.post("https://grch37.rest.ensembl.org/lookup/id",
                                  headers={"Content-Type": "application/json", "Accept": "application/json"},
                                  json={"ids": chunk}, timeout=120)
                r.raise_for_status()
                for k, v in r.json().items():
                    cache[k] = str(v["version"]) if v and v.get("version") is not None else ""
                for k in chunk:
                    cache.setdefault(k, "")
                break
            except Exception:  # noqa: BLE001
                time.sleep(3)
        else:
            for k in chunk:
                cache.setdefault(k, "")
    VER_CACHE.write_text(json.dumps(cache))
    return cache


def main() -> None:
    a, b = load(M2017), load(M50K)
    shared = a.index.intersection(b.index)
    A, B = a.loc[shared], b.loc[shared]

    def rate(col: str) -> float:
        return float((A[col] == B[col]).mean())

    cats = {
        "gene": ("Hugo_Symbol", rate("Hugo_Symbol")),
        "transcript": ("tx", rate("tx")),
        "protein": ("HGVSp_Short", rate("HGVSp_Short")),
        "variant_class": ("Variant_Classification", rate("Variant_Classification")),
    }

    # variants where transcript or protein changed
    ch = A[(A["tx"] != B["tx"]) | (A["HGVSp_Short"] != B["HGVSp_Short"])]
    changed = []
    for vk in ch.index:
        r2, r5 = A.loc[vk], B.loc[vk]
        changed.append({
            "vk": vk, "gene": r5["Hugo_Symbol"],
            "tx_2017": r2["tx"], "tx_50k": r5["tx"],
            "hgvsp_2017": r2["HGVSp_Short"], "hgvsp_50k": r5["HGVSp_Short"],
            "class_2017": r2["Variant_Classification"], "class_50k": r5["Variant_Classification"],
            "tx_changed": int(r2["tx"] != r5["tx"]),
        })
    changed.sort(key=lambda x: x["gene"])

    # per-gene transcript switches (the systematic re-assignments)
    tx_switch: dict[str, dict] = {}
    for c in changed:
        if c["tx_changed"] and c["tx_2017"] and c["tx_50k"]:
            key = (c["gene"], c["tx_2017"], c["tx_50k"])
            k = "\t".join(key)
            tx_switch.setdefault(k, {"gene": c["gene"], "tx_2017": c["tx_2017"],
                                     "tx_50k": c["tx_50k"], "n": 0})["n"] += 1
    switches = sorted(tx_switch.values(), key=lambda x: -x["n"])

    # attach current Ensembl GRCh37 .version to every displayed transcript (MAFs carry none)
    allids: set[str] = set()
    for c in changed:
        allids.update([c["tx_2017"], c["tx_50k"]])
    for s in switches:
        allids.update([s["tx_2017"], s["tx_50k"]])
    ver = fetch_versions(allids)

    def wv(tx: str) -> str:
        v = ver.get(tx, "")
        return f"{tx}.{v}" if tx and v else tx

    for c in changed:
        c["tx_2017"], c["tx_50k"] = wv(c["tx_2017"]), wv(c["tx_50k"])
    for s in switches:
        s["tx_2017"], s["tx_50k"] = wv(s["tx_2017"]), wv(s["tx_50k"])

    payload = {
        "n_2017": int(len(a)), "n_50k": int(len(b)), "n_shared": int(len(shared)),
        "n_new_in_50k": int(len(b) - len(shared)), "n_dropped_from_2017": int(len(a) - len(shared)),
        "unchanged": {k: v[1] for k, v in cats.items()},
        "n_changed": {k: int(round((1 - v[1]) * len(shared))) for k, v in cats.items()},
        "n_tx_switch_genes": len({s["gene"] for s in switches}),
        "transcript_switches": switches,
        "changed": changed,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload))
    print(f"wrote {OUT}")
    print(f"  shared {len(shared):,} | protein unchanged {cats['protein'][1]:.2%} "
          f"| transcript unchanged {cats['transcript'][1]:.2%}")
    print(f"  {len(changed):,} changed variants | {len(switches)} gene transcript-switches")


if __name__ == "__main__":
    main()
