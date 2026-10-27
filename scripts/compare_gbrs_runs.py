#!/usr/bin/env python3
"""Compare two GBRS results for the same libraries, sample by sample.

Each side is a path template for a sample's diploid quantification file, with {sample} where the sample name goes:

    python3 scripts/compare_gbrs_runs.py --sheet WORKDIR/samples.tsv \\
        --a 'WORKDIR/samples/{sample}/either/{sample}.diploid.genes.expected_read_counts' \\
        --b '/old/run/{sample}/{sample}_gbrs.quantified.diploid.genes.expected_read_counts' \\
        --out WORKDIR/comparison/either_vs_old.tsv

For each sample with both files, over genes with at least 10 expected reads on either side: total reads on each
side, the correlation of log2(total + 1), the correlation of the eight founder counts and the share of genes with
the same diplotype. Identical inputs give correlations of 1 and a diplotype share of 1. Samples missing on either
side are listed in the summary. Writes OUT and OUT with .json for the summary.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

FOUNDERS = list("ABCDEFGH")


def read(path):
    with open(path, newline="") as stream:
        return {r["locus"]: r for r in csv.DictReader(stream, delimiter="\t")}


def compare(path_a, path_b):
    a, b = read(path_a), read(path_b)
    genes = sorted(set(a) & set(b))
    ta = np.array([float(a[g]["total"]) for g in genes])
    tb = np.array([float(b[g]["total"]) for g in genes])
    used = (ta >= 10) | (tb >= 10)
    fa = np.array([[float(a[g][f]) for f in FOUNDERS] for g, u in zip(genes, used) if u]).ravel()
    fb = np.array([[float(b[g][f]) for f in FOUNDERS] for g, u in zip(genes, used) if u]).ravel()
    same = np.mean([a[g]["notes"] == b[g]["notes"] for g, u in zip(genes, used) if u])
    return {
        "genes_shared": len(genes), "genes_used": int(used.sum()),
        "total_a": round(float(ta.sum())), "total_b": round(float(tb.sum())),
        "log_total_r": round(float(np.corrcoef(np.log2(ta[used] + 1), np.log2(tb[used] + 1))[0, 1]), 5),
        "founder_counts_r": round(float(np.corrcoef(fa, fb)[0, 1]), 5),
        "same_diplotype_share": round(float(same), 4),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True, type=Path, help="a table with a sample column (samples.tsv)")
    ap.add_argument("--a", required=True, help="path template with {sample}")
    ap.add_argument("--b", required=True, help="path template with {sample}")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)
    if "{sample}" not in args.a or "{sample}" not in args.b:
        raise SystemExit("both templates need {sample}")
    with args.sheet.open(newline="") as stream:
        samples = [r["sample"] for r in csv.DictReader(stream, delimiter="\t")]
    rows, missing = [], {"a": [], "b": []}
    for s in samples:
        pa, pb = Path(args.a.format(sample=s)), Path(args.b.format(sample=s))
        if not pa.is_file():
            missing["a"].append(s)
        if not pb.is_file():
            missing["b"].append(s)
        if pa.is_file() and pb.is_file():
            rows.append({"sample": s, **compare(pa, pb)})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "genes_shared", "genes_used", "total_a", "total_b", "log_total_r", "founder_counts_r",
              "same_diplotype_share"]
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {"a": args.a, "b": args.b, "samples_compared": len(rows), "missing": missing}
    for key in ("log_total_r", "founder_counts_r", "same_diplotype_share"):
        values = [r[key] for r in rows]
        if values:
            summary[key] = {"min": min(values), "median": float(np.median(values)), "max": max(values)}
    summary["identical"] = [r["sample"] for r in rows
                            if r["total_a"] == r["total_b"] and r["same_diplotype_share"] == 1.0]
    args.out.with_suffix(".json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "identical"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
