#!/usr/bin/env python3
"""Check a GBRS sample list against the preprocessed reads and write the run's sample sheet.

The input is a tab-separated file with columns sample, sex (M or F) and
generation (the DO or CC generation used to pick the transition probabilities).
By default each sample's reads must be exactly READS_DIR/<sample>/<sample>_R1.fastq.gz
and _R2.fastq.gz, the layout written by the RNAseq_STAR HTStream stage. Names are
never matched by prefix: a missing file, an empty file, a duplicate sample or
a reads directory with samples that are not listed stops the script. Sample
folders may be links (to join pools); the sheet records each file's real path.

A sample list with r1 and r2 columns names each sample's two read files directly
(absolute, or relative to the sample list's folder); READS_DIR is then not used.
Other columns (mouse_id, genotype_group, star_library) are read by the validation
step and ignored here.

    python3 scripts/sample_sheet.py SAMPLE_LIST READS_DIR OUT_TSV [--allow-unlisted]
"""
import argparse
import csv
import re
from pathlib import Path

NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


def explicit_reads(sample, row, base):
    reads = []
    for mate in ("r1", "r2"):
        value = (row.get(mate) or "").strip()
        if not value:
            raise ValueError(f"{sample}: no {mate} path")
        path = Path(value) if Path(value).is_absolute() else base / value
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{sample}: missing or empty {path}")
        reads.append(path.resolve())
    if reads[0] == reads[1]:
        raise ValueError(f"{sample}: r1 and r2 are the same file")
    return reads


def build(sample_list, reads_dir, allow_unlisted=False):
    sample_list = Path(sample_list)
    with sample_list.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or not {"sample", "sex", "generation"} <= set(rows[0]):
        raise ValueError("Sample list needs the columns sample, sex and generation")
    explicit = {"r1", "r2"} <= set(rows[0])
    if not explicit:
        reads_dir = Path(reads_dir).resolve()
    out, seen, used = [], set(), {}
    for row in rows:
        sample, sex, gen = row["sample"].strip(), row["sex"].strip().upper(), row["generation"].strip()
        if not NAME.fullmatch(sample):
            raise ValueError(f"Unsupported sample name: {sample!r}")
        if sample in seen:
            raise ValueError(f"Sample listed twice: {sample}")
        seen.add(sample)
        if sex not in ("M", "F"):
            raise ValueError(f"{sample}: sex must be M or F, not {row['sex']!r}")
        if not gen.isdigit():
            raise ValueError(f"{sample}: generation must be a whole number")
        if explicit:
            reads = explicit_reads(sample, row, sample_list.resolve().parent)
        else:
            reads = []
            for mate in ("R1", "R2"):
                path = reads_dir / sample / f"{sample}_{mate}.fastq.gz"
                if not path.is_file() or path.is_symlink() or path.stat().st_size == 0:
                    raise ValueError(f"{sample}: missing or empty {path}")
                reads.append(path.resolve())
        for r in reads:
            if r in used:
                raise ValueError(f"{sample}: {r} is also listed for {used[r]}")
            used[r] = sample
        out.append([sample, str(reads[0]), str(reads[1]), reads[0].stat().st_size,
                    reads[1].stat().st_size, sex, int(gen)])
    if explicit:
        return out
    present = {p.name for p in reads_dir.iterdir() if p.is_dir()}
    unlisted = sorted(present - seen)
    if unlisted and not allow_unlisted:
        raise ValueError(f"{len(unlisted)} read folders are not in the sample list, e.g. {unlisted[:5]}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sample_list")
    parser.add_argument("reads_dir")
    parser.add_argument("out_tsv", type=Path)
    parser.add_argument("--allow-unlisted", action="store_true",
                        help="run a subset of the samples in READS_DIR (pilot runs)")
    args = parser.parse_args()
    rows = build(args.sample_list, args.reads_dir, args.allow_unlisted)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "r1", "r2", "r1_bytes", "r2_bytes", "sex", "generation"])
        writer.writerows(rows)
    print(f"{len(rows)} samples; reads found by exact name")


if __name__ == "__main__":
    main()
