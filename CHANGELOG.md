# Changelog

## 2026-10-26

The 14 numbered scripts in `1_GBRS_Pipeline/` are replaced by one Slurm array job, `gbrs_sample.slurm`, run through
`run.sh`. The earlier scripts remain in the repository history. The GBRS steps and their settings are unchanged:
Bowtie 1 with `-a --best --strata -v 3`, `bam2emase` with haplotypes A-H, quantification with model 4, and
reconstruction with sex-specific transition probabilities.

- Settings live in `config.sh` and a per-run `run_settings.sh` (`examples/ceres.sh`). Study paths, array sizes,
  partition and the mail address are no longer written into the scripts.
- Samples come from a sample list with sex and generation, checked by `scripts/sample_sheet.py`. Each sample's reads
  are found by their exact name; a missing, empty, duplicate or unlisted sample stops the run. Sex no longer comes
  from looking up spelling variants of the mouse ID in a spreadsheet.
- GBRS 1.0.0, EMASE and Bowtie 1 run from the JAX container `quay.io/jaxcompsci/gbrs_py3:v1.0.1`, pinned by digest.
  `scripts/fetch_references.py` downloads the reference files from Zenodo record 15091310, checks the published
  MD5 of each and records SHA-256 for every file.
- The two mates' alignments are combined by a named rule, and each rule in `MATE_RULES` gets its own results from
  one alignment: `either` (a read pair counts toward a transcript and founder where either mate aligned) gives the
  same alignments as the earlier `5_common_align_merge.slurm.py`; `both` (only where both mates aligned) is
  `emase get-common-alignments`, the step in the GBRS authors' paired-end pipeline. `scripts/merge_mates.py`
  implements `either` with GBRS's own file reader and writer and checks that both files hold the same reads in the
  same order. BAMs are no longer coordinate-sorted before `bam2emase`, which lists reads in sorted-name order anyway.
- Each sample writes into `<sample>.partial/` and becomes `<sample>/` only when every output exists; `done.json`
  records inputs, sex, generation, transition-probability file, tool versions, alignment counts and output checksums.
  A finished sample is never overwritten; a failed attempt is kept.
- `run.sh combine` (`scripts/combine_gbrs.py`) joins finished samples into gene-by-sample counts, diplotypes and
  founder counts, refuses unfinished or changed samples, and reads every cell back.
- README: the pipeline uses Bowtie 1 (the scripts always did); setup and outputs are documented.
- Tests in `tests/` cover exact-name matching and the joining step.
