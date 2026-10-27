# Changelog

## 2026-10-27

Every step that reads data or computes now runs as a Slurm job, and a run is checked before its results are used.

- `./run.sh all` submits setup (when the container or reference is missing), the per-sample array, `combine` and
  `validate`, each waiting for the one before, and records every job with its pipeline revision in
  `WORKDIR/jobs.tsv`. `./run.sh setup` (`0_Setup/gbrs_setup.slurm`) replaces `containers` and `references`;
  Apptainer does not run on Ceres login nodes. `combine` and `validate` are Slurm jobs.
- `scripts/validate_gbrs.py` (`./run.sh validate`) checks that every sample finished, that sex from Xist and Y-gene
  expression agrees with the sample list, that no two samples from different genotype groups share a genotype,
  and, when `DNA_ALLELEPROBS` or `STAR_COUNTS` is set, that each sample best matches its own mouse's DNA genotype
  and its own STAR library. `scripts/export_dna_alleleprobs.R` writes the DNA file from R/qtl2.
- `./run.sh compare` (`scripts/compare_gbrs_runs.py`) compares a rule's results with another result for the same
  libraries, sample by sample, for reruns and version changes.
- The sample list may name read files in `r1` and `r2` columns instead of using `READS_DIR`, and may carry
  `mouse_id`, `genotype_group` and `star_library` for the checks. A read file listed for two samples stops the run.
- The pipeline runs under GBRS 1.0.0 and 1.1.0: alignments are converted with `emase bam2emase` (in 1.0.0 the same
  function as `gbrs bam2emase`, which 1.1.0 removed), and `scripts/merge_mates.py` passes file names by position
  (1.1.0 renamed the keyword).
- The record step is `scripts/record_sample.py` (was inline in the job script) with tests. It reads bowtie 1.3's
  summary line ('reads with at least one alignment'); the inline version expected older wording and stopped
  every sample at its last step.
- `done.json` records the pipeline revision (`git describe`, or the `REVISION` file of a copy without `.git`).
- Defaults from measured runs: 64 GB and 12 hours per sample. Documentation: `docs/method.md` (GBRS and STAR,
  inputs, mate rules, versions), `docs/validation.md` and `docs/runtime.md`; `examples/other_slurm.sh`.

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
- The steps after combining the mates are in `1_GBRS_Pipeline/gbrs_downstream.sh`, which sets `GBRS_DATA` to the
  reference directory and checks each step's output files before the next step: `gbrs reconstruct` needs
  `GBRS_DATA/ref.fa.fai` (linked by `fetch_references.py`) and, without it, prints an error but exits 0. The same
  script can finish a sample from an existing compressed alignment file.
- `run.sh combine` (`scripts/combine_gbrs.py`) joins finished samples into gene-by-sample counts, diplotypes and
  founder counts, refuses unfinished or changed samples, and reads every cell back.
- README: the pipeline uses Bowtie 1 (the scripts always did); setup and outputs are documented.
- Tests in `tests/` cover exact-name matching and the joining step.
