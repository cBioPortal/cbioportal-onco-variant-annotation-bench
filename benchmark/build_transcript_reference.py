"""Build a multi-source canonical-transcript reference, per track.

Source: genome-nexus-importer's per-HGNC canonical export, which lists, per gene,
the canonical transcript each source picks (Ensembl, Genome Nexus, UniProt, MSKCC)
plus MANE Select — the "different canonical transcripts going around". The GRCh37
export (ensembl111) carries a .version per transcript; the GRCh38 export
(ensembl95) does not, so versions are simply blank there.

  GRCh37 (MSK-IMPACT truth): scoped to the benchmark panel + MSK-override genes;
    MSKCC column IS the truth.
  GRCh38 (TCGA-GDC truth): scoped to the SAME gene set as GRCh37 so the two tracks
    are directly comparable; the truth is the GDC's VEP/GENCODE pick, so MANE (not
    MSKCC) is the closest reference here.

Output: website/data/<track>/transcript_reference.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REFDIR = ROOT / "data" / "transcript_ref"
IMPORTER = "https://github.com/genome-nexus/genome-nexus-importer/blob/master/data"

TRACKS = {
    "grch37": {
        "canon": REFDIR / "canonical_per_hgnc_grch37_ensembl111.txt",
        "truth": ROOT / "data" / "truth" / "truth.parquet",
        "msk_override": ROOT / "data" / "msk_overrides" / "isoform_overrides_at_mskcc_grch37.txt",
        "scope_from": None,
        "out": ROOT / "website" / "data" / "grch37" / "transcript_reference.json",
        "source": "genome-nexus-importer · grch37_ensembl111 (Ensembl 111, GRCh37)",
        "source_url": f"{IMPORTER}/grch37_ensembl111/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt",
        "has_version": True,
        "truth_source": "mskcc",
    },
    "grch38": {
        "canon": REFDIR / "canonical_per_hgnc_grch38_ensembl95.txt",
        "truth": ROOT / "data" / "grch38" / "truth" / "truth.parquet",
        "msk_override": None,
        "scope_from": ROOT / "website" / "data" / "grch37" / "transcript_reference.json",
        "out": ROOT / "website" / "data" / "grch38" / "transcript_reference.json",
        "source": "genome-nexus-importer · grch38_ensembl95 (Ensembl 95, GRCh38)",
        "source_url": f"{IMPORTER}/grch38_ensembl95/export/ensembl_biomart_canonical_transcripts_per_hgnc.txt",
        "has_version": False,
        # the ensembl95 export has no transcript .version; fill from Ensembl GRCh38 REST.
        "version_map": REFDIR / "enst_versions_grch38.json",
        "truth_source": "mane",
    },
}


def tv(tx: str, ver: str) -> str:
    tx, ver = (tx or "").strip(), (ver or "").strip()
    return f"{tx}.{ver}" if tx and ver else tx


def base(x: str) -> str:
    return x.split(".")[0]


def build(track: str) -> None:
    cfg = TRACKS[track]
    df = pd.read_csv(cfg["canon"], sep="\t", dtype=str, na_filter=False, low_memory=False)
    hasv = cfg["has_version"]
    # GRCh38 has no inline version columns; fill from an Ensembl REST version map.
    vmap = {}
    vpath = cfg.get("version_map")
    if vpath and Path(vpath).exists():
        vmap = json.loads(Path(vpath).read_text())

    def col(r, name):
        return getattr(r, name, "") or ""

    def ver(r, name):
        if hasv:
            return getattr(r, f"{name}_version", "")
        tx = col(r, name)
        return vmap.get(tx.split(".")[0], "") if tx else ""

    # gene scope
    if cfg["scope_from"]:
        want = {row["gene"] for row in json.loads(Path(cfg["scope_from"]).read_text())["rows"]}
        panel = set(pd.read_parquet(cfg["truth"])["hugo_symbol"].dropna())
    else:
        panel = set(pd.read_parquet(cfg["truth"])["hugo_symbol"].dropna())
        override_genes = set(
            pd.read_csv(cfg["msk_override"], sep="\t", dtype=str).fillna("")["gene_name"]
        )
        want = {g for g in (panel | override_genes) if g}

    rows, counts = [], {"msk_eq_ensembl": 0, "msk_override": 0, "msk_eq_mane": 0}
    for r in df.itertuples(index=False):
        gene = r.hgnc_symbol
        if gene not in want:
            continue
        ens = tv(col(r, "ensembl_canonical_transcript"), ver(r, "ensembl_canonical_transcript"))
        gn = tv(col(r, "genome_nexus_canonical_transcript"), ver(r, "genome_nexus_canonical_transcript"))
        uni = tv(col(r, "uniprot_canonical_transcript"), ver(r, "uniprot_canonical_transcript"))
        msk_tx = tv(col(r, "mskcc_canonical_transcript"), ver(r, "mskcc_canonical_transcript"))
        mane_raw = col(r, "mane_select")
        mane_enst = mane_raw.split(",")[0].strip()
        mane_refseq = mane_raw.split(",")[1].strip() if "," in mane_raw else ""

        rec = {
            "gene": gene,
            "in_panel": gene in panel,
            "mane": mane_enst,
            "mane_refseq": mane_refseq,
            "ensembl": ens,
            "ensembl_why": col(r, "ensembl_canonical_transcript_explanation"),
            "gn": gn,
            "gn_why": col(r, "genome_nexus_canonical_transcript_explanation"),
            "uniprot": uni,
            "uniprot_why": col(r, "uniprot_canonical_transcript_explanation"),
            "mskcc": msk_tx,
            "mskcc_why": col(r, "mskcc_canonical_transcript_explanation"),
        }
        rec["type"] = "msk_eq_ensembl" if (msk_tx and base(msk_tx) == base(ens)) else "msk_override"
        rec["msk_eq_mane"] = bool(mane_enst and msk_tx and base(msk_tx) == base(mane_enst))
        counts[rec["type"]] += 1
        if rec["msk_eq_mane"]:
            counts["msk_eq_mane"] += 1
        rows.append(rec)

    rows.sort(key=lambda x: x["gene"])
    payload = {
        "n": len(rows),
        "counts": counts,
        "has_version": hasv or bool(vmap),
        "version_source": None if hasv else ("Ensembl GRCh38 REST (current)" if vmap else None),
        "truth_source": cfg["truth_source"],
        "source": cfg["source"],
        "source_url": cfg["source_url"],
        "rows": rows,
    }
    cfg["out"].parent.mkdir(parents=True, exist_ok=True)
    cfg["out"].write_text(json.dumps(payload))
    print(f"[{track}] wrote {cfg['out']}  ({len(rows)} genes)  counts: {counts}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--track", default="grch37", choices=list(TRACKS))
    args = ap.parse_args()
    build(args.track)


if __name__ == "__main__":
    main()
