# Run time, memory and storage

Measured on Ceres in October 2026 with DO_Chln liver libraries of about 20 million read pairs, 16 CPUs per sample:

| Step | Time per sample | Peak memory |
| --- | --- | --- |
| Align R1, convert to EMASE | 45 to 75 minutes | |
| Align R2, convert to EMASE | 45 to 70 minutes | |
| Combine the reads and compress, for each mate rule | 3 to 6 minutes | 24 GB |
| Quantify, reconstruct, interpolate, export, plot, for each mate rule | 2 minutes | under 1 GB |
| Whole sample, one mate rule | 1.5 to 2.6 hours | 25 to 32 GB |
| Whole sample, both mate rules | 2.2 to 3.1 hours | 19 to 31 GB |

GBRS 1.1.0 and 1.0.0 took the same time and memory.

The whole DO_Chln cohort (388 libraries, both mate rules, 60 at a time with 16 CPUs and 64 GB each) took about 20
hours from submission to the last sample on Ceres (October 26 and 27, 2026).

The defaults (`MEMORY=64G`, `TIME_LIMIT=12:00:00`) leave room for larger libraries. Alignment time grows with the
number of reads. `THROTTLE` sets how many samples run at once.

Storage per sample is about 200 MB for each mate rule (counts, genotype files, plot and the compressed alignment
file). The reference takes about 5.7 GB once unpacked and the container about 0.5 GB; several runs can share both through `REF_DIR`
and `GBRS_SIF`. Temporary BAM and EMASE files go to node-local scratch (`TMPDIR`) and are removed.

## Ceres

- Only `run.sh` runs on the login node; every step that reads data or computes is a Slurm job. Apptainer does not
  run on Ceres login nodes, which is why `./run.sh setup` is a job.
- `/90daydata` removes files not read for 90 days, and its cleanup also removed empty folders while a run was going
  (October 27, 2026). The pipeline creates each output folder together with its first file. Copy accepted results
  to `/project` (or Box) before 90 days pass.
- Many short ssh sessions in a row can make the login nodes refuse new connections for a while; reuse one
  connection (OpenSSH `ControlMaster`) when scripting.

## Other Slurm clusters

Copy `examples/other_slurm.sh`: set the partition and account, and either the module names or `SOFTWARE_MODE=path`
with `apptainer` (or `gbrs`, `emase`, `bowtie`), `samtools` and `python3` on `PATH`. The container brings GBRS,
EMASE, Bowtie 1 and the Python packages. Clusters without Slurm are not supported.
