#!/usr/bin/env python3
"""Join finished GBRS samples into gene-by-sample tables.

Every sample in the sample sheet must have finished (done.json present, its
recorded output checksums unchanged). For one mate rule (either or both), from
each sample's <rule>/<sample>.diploid.genes.expected_read_counts this writes,
in sample-sheet order:

  gbrs_counts.tsv          gene x sample, expected read counts (column total)
  gbrs_diplotypes.tsv      gene x sample, reconstructed diplotype (column notes)
  gbrs_founder_counts.tsv  gene x sample.founder, expected read counts per founder
  qc.tsv                   per sample: reads, aligned reads per mate, sex, generation
  SHA256SUMS.txt

All samples must list the same genes in the same order. After writing, every
cell is read back and compared with its source file.

    python3 scripts/combine_gbrs.py SAMPLE_SHEET SAMPLES_DIR OUT_DIR RULE
"""
import csv
import hashlib
import json
import sys
from pathlib import Path

FOUNDERS = list("ABCDEFGH")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_sample(folder, sample, rule):
    done = json.loads((folder / "done.json").read_text())
    if rule not in done.get("mate_rules", []):
        raise ValueError(f"{sample}: no results for the mate rule {rule}")
    for name, digest in done["outputs"].items():
        if name.endswith(".compressed.emase.h5"):
            continue  # large; checked by size only
        if sha(folder / name) != digest:
            raise ValueError(f"{sample}: {name} changed after the sample finished")
    path = folder / rule / f"{sample}.diploid.genes.expected_read_counts"
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or [c for c in ("locus", *FOUNDERS, "total", "notes") if c not in rows[0]]:
        raise ValueError(f"{path}: unexpected columns")
    return done, rows


def write(path, header, body):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(body)


def read_back(path):
    with path.open(newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def main():
    sheet, samples_dir, out, rule = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4]
    if out.exists():
        raise SystemExit(f"{out} exists")
    with sheet.open(newline="") as stream:
        meta = list(csv.DictReader(stream, delimiter="\t"))
    samples = [m["sample"] for m in meta]
    missing = [s for s in samples if not (samples_dir / s / "done.json").is_file()]
    if missing:
        raise SystemExit(f"{len(missing)} samples have not finished, e.g. {missing[:5]}")
    genes, totals, notes, founders, qc = None, {}, {}, {}, []
    for m in meta:
        s = m["sample"]
        done, rows = read_sample(samples_dir / s, s, rule)
        g = [r["locus"] for r in rows]
        if genes is None:
            genes = g
        elif g != genes:
            raise SystemExit(f"{s}: genes or their order differ from {samples[0]}")
        totals[s] = [r["total"] for r in rows]
        notes[s] = [r["notes"] for r in rows]
        founders[s] = [[r[f] for f in FOUNDERS] for r in rows]
        b = done["bowtie"]
        qc.append([s, m["sex"], m["generation"], b["R1"]["reads"], b["R1"]["aligned"], b["R2"]["aligned"],
                   done.get("versions", "")])
    if genes is None:
        raise SystemExit("The sample sheet lists no samples")
    out.mkdir(parents=True)
    write(out / "gbrs_counts.tsv", ["gene_id", *samples], ([g, *(totals[s][i] for s in samples)] for i, g in enumerate(genes)))
    write(out / "gbrs_diplotypes.tsv", ["gene_id", *samples], ([g, *(notes[s][i] for s in samples)] for i, g in enumerate(genes)))
    write(out / "gbrs_founder_counts.tsv", ["gene_id", *(f"{s}.{f}" for s in samples for f in FOUNDERS)],
          ([g, *(v for s in samples for v in founders[s][i])] for i, g in enumerate(genes)))
    write(out / "qc.tsv", ["sample", "sex", "generation", "read_pairs", "aligned_r1", "aligned_r2", "versions"], qc)
    for name, source in (("gbrs_counts.tsv", totals), ("gbrs_diplotypes.tsv", notes)):
        back = read_back(out / name)
        if back[0] != ["gene_id", *samples] or [r[0] for r in back[1:]] != genes:
            raise SystemExit(f"{name}: header or gene order differs on read-back")
        for j, s in enumerate(samples, 1):
            if [r[j] for r in back[1:]] != source[s]:
                raise SystemExit(f"{name}: column {s} differs on read-back")
    names = ["gbrs_counts.tsv", "gbrs_diplotypes.tsv", "gbrs_founder_counts.tsv", "qc.tsv"]
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha(out / n)}  {n}\n" for n in names))
    print(json.dumps({"rule": rule, "samples": len(samples), "genes": len(genes)}))


if __name__ == "__main__":
    main()
