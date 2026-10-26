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
- The two mates are combined with `emase get-common-alignments`, as in the JAX GBRS pipeline, instead of the
  separate `5_common_align_merge.slurm.py`. BAMs are no longer coordinate-sorted before `bam2emase`.
- Each sample writes into `<sample>.partial/` and becomes `<sample>/` only when every output exists; `done.json`
  records inputs, sex, generation, transition-probability file, tool versions, alignment counts and output checksums.
  A finished sample is never overwritten; a failed attempt is kept.
- `run.sh combine` (`scripts/combine_gbrs.py`) joins finished samples into gene-by-sample counts, diplotypes and
  founder counts, refuses unfinished or changed samples, and reads every cell back.
- README: the pipeline uses Bowtie 1 (the scripts always did); setup and outputs are documented.
- Tests in `tests/` cover exact-name matching and the joining step.
