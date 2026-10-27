#!/usr/bin/env python3
"""Combine one sample's two mate alignment files (EMASE format) into one. Runs inside the GBRS container.

Rules (which transcript and founder a read pair counts toward):
  either  a pair is compatible wherever either mate aligned; the rule of the lab's 2025 GBRS results
  both    only where both mates aligned; the same result as `emase get-common-alignments`, the step in the
          GBRS authors' paired-end pipeline (the pipeline itself calls that command for this rule)

Both files must hold the same transcripts and the same reads in the same order. bam2emase lists every read of the
BAM in sorted-name order, so the two mates of one sample always do; anything else stops the script.

    python3 merge_mates.py either|both R1.emase.h5 R2.emase.h5 OUT.emase.h5
"""
import sys

import numpy as np
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM

rule, r1, r2, out = sys.argv[1:5]
if rule not in ("either", "both"):
    raise SystemExit("rule must be either or both")
# Positional arguments: GBRS 1.1.0 renamed the file keyword (h5file to h5_file).
a, b = APM(None, r1), APM(None, r2)
if tuple(a.shape) != tuple(b.shape):
    raise SystemExit(f"shapes differ: {a.shape} and {b.shape}")
if list(a.lname) != list(b.lname):
    raise SystemExit("transcript lists differ")
if not np.array_equal(np.asarray(a.rname), np.asarray(b.rname)):
    raise SystemExit("read lists differ; the two files are not mates of one sample")
for h in range(a.num_haplotypes):
    x, y = a.data[h].tocsc(), b.data[h].tocsc()
    merged = (x + y) if rule == "either" else x.multiply(y)
    merged = merged.tocsc()
    merged.eliminate_zeros()
    merged.data[:] = 1  # alignment incidence, as bam2emase writes it
    a.data[h] = merged
a.save(out, complib="zlib")
nnz = sum(a.data[h].nnz for h in range(a.num_haplotypes))
print(f"{rule}: {a.num_reads} reads, {a.num_loci} transcripts, {nnz} read-transcript-founder alignments")
