#!/usr/bin/env python3
"""Check one mate rule of a GBRS run and write a summary.

    python3 scripts/validate_gbrs.py --sheet WORKDIR/samples.tsv --samples WORKDIR/samples --rule either \\
        --grid REF_DIR/ref.genome_grid.GRCm39.tsv --out WORKDIR/validation/either \\
        [--sample-list SAMPLE_LIST] [--dna DNA_ALLELEPROBS.tsv.gz] [--star STAR_COUNTS.tsv] [--same-mouse 0.9]

Checks (docs/validation.md explains each one and how to read the numbers):
  finished    every sample in the sheet has done.json and this rule's outputs
  sex         Xist against four Y-chromosome genes agrees with the sex used for the transition probabilities;
              strong signal from both marks a possible mixture of two mice
  duplicates  every pair of samples is compared on their reconstructed founder genotypes: samples from different
              genotype groups should not look like one mouse, and samples sharing a group should
  dna         with --dna: each sample's genotypes best match its own mouse's DNA founder probabilities
  star        with --star: each sample's GBRS gene counts best match its own library in a STAR count matrix

Optional sample-list columns: mouse_id (the mouse's ID in the DNA file; default the sample name),
genotype_group (samples expected to share a genotype, such as mice of one inbred CC strain or one mouse sequenced
twice; default each sample alone) and star_library (default the sample name).

Writes OUT/summary.json and one TSV per check, and exits 1 when a check fails.
"""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import numpy as np

FOUNDERS = list("ABCDEFGH")
AUTOSOMES = [str(c) for c in range(1, 20)]
XIST = "ENSMUSG00000086503"
Y_GENES = ["ENSMUSG00000069045", "ENSMUSG00000056673", "ENSMUSG00000069049", "ENSMUSG00000068457"]  # Ddx3y Kdm5d Eif2s3y Uty
MIXED_XIST_CPM = 100.0  # a male library has under 10 Xist counts per million; a female one has over 1,000
MIXED_Y_CPM = 20.0  # a female library has no Y-gene counts; a male one has about 200 per million


def read_tsv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def overlap(one, many):
    """Founder-probability overlap of one sample (points x 8) with each of many (n x points x 8).

    At each point, the sum over founders of the smaller probability; averaged over points. 1 means identical.
    """
    return np.minimum(one[None], many).sum(axis=2).mean(axis=1)


def load_grid(path, every):
    """Row numbers and positions (Mb) of the autosomal grid points used, keeping every Nth point."""
    rows = read_tsv(path)
    keep, chrom, pos = [], [], []
    for i, row in enumerate(rows):
        if row["chr"] in AUTOSOMES:
            keep.append(i)
            chrom.append(row["chr"])
            pos.append(float(row.get("bp") or row["pos"]) / 1e6)
    keep, chrom, pos = np.array(keep)[::every], np.array(chrom)[::every], np.array(pos)[::every]
    return len(rows), keep, chrom, pos


def load_samples(sheet, samples_dir, rule, grid_rows, keep):
    """Per sample: done record, diploid gene totals and founder probabilities at the kept grid points."""
    finished, unfinished, genes = {}, {}, None
    for row in sheet:
        s = row["sample"]
        folder = samples_dir / s
        done_path = folder / "done.json"
        if not done_path.is_file():
            unfinished[s] = "no done.json"
            continue
        done = json.loads(done_path.read_text())
        if rule not in done.get("mate_rules", []):
            unfinished[s] = f"no results for rule {rule}"
            continue
        counts = folder / rule / f"{s}.diploid.genes.expected_read_counts"
        probs = folder / rule / f"{s}.gbrs.interpolated.genoprobs.tsv"
        missing = [p.name for p in (counts, probs) if not p.is_file() or p.stat().st_size == 0]
        if missing:
            unfinished[s] = "missing " + ", ".join(missing)
            continue
        table = read_tsv(counts)
        g = [r["locus"] for r in table]
        if genes is None:
            genes = g
        elif g != genes:
            raise SystemExit(f"{s}: genes or their order differ from the first sample")
        p = np.loadtxt(probs, comments="#", dtype=np.float32)
        if p.shape != (grid_rows, 8):
            raise SystemExit(f"{s}: genotype probabilities have shape {p.shape}, expected ({grid_rows}, 8)")
        finished[s] = {"done": done, "totals": np.array([float(r["total"]) for r in table]), "probs": p[keep]}
    return finished, unfinished, genes or []


def check_sex(sheet, finished, genes):
    index = {g: i for i, g in enumerate(genes)}
    rows, failed = [], []
    for row in sheet:
        s = row["sample"]
        if s not in finished:
            continue
        t = finished[s]["totals"]
        per_million = t.sum() / 1e6
        x = t[index[XIST]] / per_million if XIST in index else 0.0
        y = sum(t[index[g]] for g in Y_GENES if g in index) / per_million
        if x >= MIXED_XIST_CPM and y >= MIXED_Y_CPM:
            call = "mixed"
        else:
            call = "F" if x > y else "M"
        ok = call == row["sex"]
        rows.append([s, row["sex"], f"{x:.1f}", f"{y:.1f}", call, ok])
        if not ok:
            failed.append({"sample": s, "sheet_sex": row["sex"], "expression_sex": call})
    return rows, failed


def check_duplicates(names, probs, group, threshold):
    n = len(names)
    scores = np.full((n, n), np.nan, dtype=np.float32)
    for i in range(n):
        scores[i] = overlap(probs[i], probs)
        scores[i, i] = np.nan
    rows, pairs = [], []
    for i, s in enumerate(names):
        if n < 2:
            rows.append([s, group[s], "", "", "", "single sample"])
            continue
        j = int(np.nanargmax(scores[i]))
        same = group[names[j]] == group[s]
        rows.append([s, group[s], names[j], f"{scores[i, j]:.4f}", same, ""])
        for k in range(i + 1, n):
            same_group = group[names[k]] == group[s]
            if not same_group and scores[i, k] >= threshold:
                pairs.append({"a": s, "b": names[k], "score": round(float(scores[i, k]), 4),
                              "reason": "different genotype groups but the same genotype"})
            elif same_group and scores[i, k] < threshold:
                pairs.append({"a": s, "b": names[k], "score": round(float(scores[i, k]), 4),
                              "reason": "same genotype group but different genotypes"})
    flagged = {p["a"] for p in pairs} | {p["b"] for p in pairs}
    for r in rows:
        r[5] = r[5] or ("flagged" if r[0] in flagged else "ok")
    other = np.array([[group[a] != group[b] for b in names] for a in names], dtype=bool).reshape(n, n)
    stats = {}
    if other.any():
        closest = [float(scores[i][other[i]].max()) for i in range(n) if other[i].any()]
        stats["other_groups_median"] = round(float(np.median(scores[other])), 4)
        stats["closest_other_group_median"] = round(float(np.median(closest)), 4)
        stats["closest_other_group_max"] = round(float(np.max(closest)), 4)
    return rows, pairs, stats


def load_dna(path, chrom, pos):
    """DNA founder probabilities at the nearest marker to each kept grid point: (mice, points, 8)."""
    by_chr = {}
    for row in read_tsv(path):
        c = row["chr"]
        if c not in AUTOSOMES:
            continue
        by_chr.setdefault(c, []).append((row["mouse_id"], float(row["pos_mb"]), [float(row[f]) for f in FOUNDERS]))
    mice = sorted({m for rows in by_chr.values() for m, _, _ in rows})
    mouse_index = {m: i for i, m in enumerate(mice)}
    out = np.zeros((len(mice), len(pos), 8), dtype=np.float32)
    for c in sorted(set(chrom)):
        if c not in by_chr:
            raise SystemExit(f"DNA file has no markers on chromosome {c}")
        markers = sorted({p for _, p, _ in by_chr[c]})
        marker_index = {p: k for k, p in enumerate(markers)}
        table = np.full((len(mice), len(markers), 8), np.nan, dtype=np.float32)
        for m, p, values in by_chr[c]:
            table[mouse_index[m], marker_index[p]] = values
        if np.isnan(table).any():
            raise SystemExit(f"DNA file: not every mouse has every marker on chromosome {c}")
        points = np.where(chrom == c)[0]
        mpos = np.array(markers)
        nearest = np.abs(mpos[None, :] - pos[points][:, None]).argmin(axis=1)
        out[:, points] = table[:, nearest]
    return mice, out


def best_two(scores, labels):
    order = np.argsort(-scores)
    second = order[1] if len(order) > 1 else order[0]
    return labels[order[0]], float(scores[order[0]]), labels[second], float(scores[second])


def check_dna(names, probs, mouse_of, mice, dna):
    index = {m: i for i, m in enumerate(mice)}
    rows, failed, unchecked = [], [], []
    for i, s in enumerate(names):
        own = mouse_of[s]
        scores = overlap(probs[i], dna)
        best, b, second, c = best_two(scores, mice)
        own_score = float(scores[index[own]]) if own in index else float("nan")
        if own not in index:
            status = "mouse not in DNA file"
            unchecked.append(s)
        elif best == own:
            status = "ok"
        else:
            status = "best match is another mouse"
            failed.append({"sample": s, "mouse_id": own, "best_mouse": best, "score_own": round(own_score, 4),
                           "score_best": round(b, 4)})
        rows.append([s, own, f"{own_score:.4f}", best, f"{b:.4f}", second, f"{c:.4f}", status])
    return rows, failed, unchecked


def check_star(names, finished, genes, star_path, library_of, top=2000):
    opener = gzip.open if str(star_path).endswith(".gz") else open
    gbrs_index = {g: i for i, g in enumerate(genes)}
    star_genes, values = [], []
    with opener(star_path, "rt", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        libraries = next(reader)[1:]
        for row in reader:
            gene = row[0].split(".")[0]
            if gene in gbrs_index:
                star_genes.append(gene)
                values.append([float(x) for x in row[1:]])
    if len(star_genes) < top:
        raise SystemExit(f"STAR matrix shares only {len(star_genes)} genes with GBRS")
    star = np.array(values)
    gb = np.array([finished[s]["totals"][[gbrs_index[g] for g in star_genes]] for s in names]).T
    star_cpm = star / star.sum(axis=0) * 1e6
    expressed = (star_cpm > 1).mean(axis=1) > 0.5
    star_log = np.log2(star_cpm[expressed] + 1)
    gbrs_log = np.log2(gb[expressed] / gb[expressed].sum(axis=0) * 1e6 + 1)
    star_c = star_log - star_log.mean(axis=1, keepdims=True)
    # Centre each side on its own samples when there are enough, which removes the genes' overall
    # difference between the two pipelines; with few GBRS samples, centre them on the STAR means.
    centre = gbrs_log.mean(axis=1, keepdims=True) if len(names) >= 10 else star_log.mean(axis=1, keepdims=True)
    gbrs_c = gbrs_log - centre
    variable = np.argsort(-star_c.var(axis=1))[:top]
    a = star_c[variable]
    a = (a - a.mean(axis=0)) / a.std(axis=0)
    rows, failed, unchecked = [], [], []
    lib_index = {l: i for i, l in enumerate(libraries)}
    for i, s in enumerate(names):
        v = gbrs_c[variable, i]
        v = (v - v.mean()) / v.std()
        r = a.T @ v / len(v)
        best, b, second, c = best_two(r, libraries)
        own = library_of[s]
        own_r = float(r[lib_index[own]]) if own in lib_index else float("nan")
        if own not in lib_index:
            status = "library not in STAR matrix"
            unchecked.append(s)
        elif best == own:
            status = "ok"
        else:
            status = "best match is another library"
            failed.append({"sample": s, "star_library": own, "best_library": best, "r_own": round(own_r, 4),
                           "r_best": round(b, 4)})
        rows.append([s, own, f"{own_r:.4f}", best, f"{b:.4f}", second, f"{c:.4f}", status])
    method = "each side centred on its own samples" if len(names) >= 10 else "centred on the STAR means"
    return rows, failed, unchecked, method


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True, type=Path, help="WORKDIR/samples.tsv written by run.sh sheet")
    ap.add_argument("--samples", required=True, type=Path, help="WORKDIR/samples")
    ap.add_argument("--rule", required=True, choices=["either", "both"])
    ap.add_argument("--grid", required=True, type=Path, help="REF_DIR/ref.genome_grid.GRCm39.tsv")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sample-list", type=Path, help="SAMPLE_LIST, for mouse_id, genotype_group, star_library")
    ap.add_argument("--dna", type=Path, help="DNA founder probabilities: mouse_id, chr, pos_mb, A..H")
    ap.add_argument("--star", type=Path, help="STAR gene-by-library count matrix")
    ap.add_argument("--same-mouse", type=float, default=0.9, help="genotype overlap counted as one mouse")
    ap.add_argument("--every", type=int, default=10, help="use every Nth genome grid point")
    args = ap.parse_args(argv)

    sheet = read_tsv(args.sheet)
    extra = {r["sample"]: r for r in read_tsv(args.sample_list)} if args.sample_list else {}
    def column(s, name):
        value = (extra.get(s, {}).get(name) or "").strip()
        return value or s
    grid_rows, keep, chrom, pos = load_grid(args.grid, args.every)
    finished, unfinished, genes = load_samples(sheet, args.samples, args.rule, grid_rows, keep)
    names = [r["sample"] for r in sheet if r["sample"] in finished]
    probs = np.stack([finished[s]["probs"] for s in names]) if names else np.zeros((0, len(keep), 8))
    args.out.mkdir(parents=True, exist_ok=True)
    checks = {"finished": {"status": "fail" if unfinished else "pass", "samples": len(sheet),
                           "finished": len(finished), "unfinished": unfinished}}

    sex_rows, sex_failed = check_sex(sheet, finished, genes)
    write_tsv(args.out / "sex.tsv", ["sample", "sheet_sex", "xist_cpm", "y_genes_cpm", "expression_sex", "agrees"],
              sex_rows)
    checks["sex"] = {"status": "fail" if sex_failed else "pass", "failed": sex_failed, "file": "sex.tsv"}

    group = {s: column(s, "genotype_group") for s in names}
    dup_rows, pairs, stats = check_duplicates(names, probs, group, args.same_mouse)
    write_tsv(args.out / "duplicates.tsv",
              ["sample", "genotype_group", "closest_sample", "score", "closest_in_same_group", "status"], dup_rows)
    checks["duplicates"] = {"status": "fail" if pairs else "pass", "same_mouse_score": args.same_mouse,
                            "flagged_pairs": pairs, **stats, "file": "duplicates.tsv"}

    if args.dna:
        mice, dna = load_dna(args.dna, chrom, pos)
        mouse_of = {s: column(s, "mouse_id") for s in names}
        rows, failed, unchecked = check_dna(names, probs, mouse_of, mice, dna)
        write_tsv(args.out / "dna_identity.tsv", ["sample", "mouse_id", "score_own", "best_mouse", "score_best",
                                                   "second_mouse", "score_second", "status"], rows)
        checks["dna"] = {"status": "fail" if failed else "pass", "failed": failed, "unchecked": unchecked,
                         "mice_in_file": len(mice), "file": "dna_identity.tsv"}
    else:
        checks["dna"] = {"status": "skipped"}

    if args.star:
        library_of = {s: column(s, "star_library") for s in names}
        rows, failed, unchecked, method = check_star(names, finished, genes, args.star, library_of)
        write_tsv(args.out / "star_identity.tsv", ["sample", "star_library", "r_own", "best_library", "r_best",
                                                    "second_library", "r_second", "status"], rows)
        checks["star"] = {"status": "fail" if failed else "pass", "failed": failed, "unchecked": unchecked,
                          "method": method, "file": "star_identity.tsv"}
    else:
        checks["star"] = {"status": "skipped"}

    revisions = sorted({f["done"].get("pipeline_revision", "not recorded") for f in finished.values()})
    summary = {"rule": args.rule, "samples": len(sheet), "grid_points_used": int(len(keep)),
               "pipeline_revisions": revisions, "checks": checks,
               "passed": all(c["status"] != "fail" for c in checks.values())}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    for name, c in checks.items():
        print(f"{name:11s} {c['status']}")
    print(f"summary: {args.out / 'summary.json'}")
    return 0 if summary["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
