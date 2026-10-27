#!/usr/bin/env python3
"""Write a finished sample's done.json: inputs, settings, tool versions, alignment counts and output checksums.

Called by 1_GBRS_Pipeline/gbrs_sample.slurm after every mate rule's outputs exist. Stops without writing when a
bowtie log lacks its summary, the two mates' read counts differ, or an expected output is missing or empty.

    python3 scripts/record_sample.py DIR SAMPLE R1 R2 SEX GENERATION TPROB IMAGE VERSIONS "RULES" REVISION
"""
import hashlib
import json
import re
import sys
from pathlib import Path

# Outputs every mate rule must have, as <rule>/<sample>.<name>.
NEEDED = ["multiway.genes.tpm", "genotypes.tsv", "genoprobs.npz", "diploid.genes.expected_read_counts",
          "interpolated.genoprobs.npz", "gbrs.interpolated.genoprobs.tsv", "compressed.emase.h5"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def bowtie_counts(text):
    """Reads processed and reads aligned from a bowtie 1 log.

    Bowtie 1.3 writes '# reads with at least one alignment: N (x%)'; some versions write 'reported alignment'.
    """
    processed = re.search(r"reads processed: (\d+)", text)
    aligned = re.search(r"reads with at least one (?:reported )?alignment: (\d+)", text)
    if not processed or not aligned:
        return None
    return {"reads": int(processed.group(1)), "aligned": int(aligned.group(1))}


def record(folder, sample, r1, r2, sex, generation, tprob, image, versions, rules, revision):
    folder = Path(folder)
    stats = {}
    for mate in ("R1", "R2"):
        counts = bowtie_counts((folder / f"{sample}.bowtie_{mate}.log").read_text())
        if counts is None:
            raise SystemExit(f"bowtie log for {mate} lacks its summary")
        stats[mate] = counts
    if stats["R1"]["reads"] != stats["R2"]["reads"]:
        raise SystemExit("R1 and R2 read counts differ")
    for rule in rules.split():
        for name in NEEDED:
            p = folder / rule / f"{sample}.{name}"
            if not p.is_file() or p.stat().st_size == 0:
                raise SystemExit(f"missing output {rule}/{p.name}")
    data = {"sample": sample, "reads": {"r1": r1, "r2": r2}, "sex": sex, "generation": int(generation),
            "transition_probabilities": tprob, "container": image, "versions": versions, "bowtie": stats,
            "pipeline_revision": revision, "mate_rules": rules.split(),
            "outputs": {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob("*")) if p.is_file()}}
    (folder / "done.json").write_text(json.dumps(data, indent=1) + "\n")
    return data


if __name__ == "__main__":
    if len(sys.argv) != 12:
        raise SystemExit(__doc__)
    record(*sys.argv[1:])
