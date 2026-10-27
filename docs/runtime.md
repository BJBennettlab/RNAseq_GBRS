# Run time, memory and storage

Measured on Ceres in October 2026 with DO_Chln liver libraries of about 20 million read pairs, 16 CPUs per sample:

| Step | Time per sample | Peak memory |
| --- | --- | --- |
| Align R1, convert to EMASE | 45 to 75 minutes | |
| Align R2, convert to EMASE | 45 to 70 minutes | |
| Combine the reads and compress, for each mate rule | 3 to 6 minutes | 24 GB |
| Quantify, reconstruct, interpolate, export, plot, for each mate rule | 2 minutes | under 1 GB |
| Whole sample | 1.5 to 2.6 hours | 25 to 32 GB |

The defaults (`MEMORY=64G`, `TIME_LIMIT=12:00:00`) leave room for larger libraries. Alignment time grows with the
number of reads. `THROTTLE` sets how many samples run at once.

Storage per sample is about 200 MB for each mate rule (counts, genotype files, plot and the compressed alignment
file). The reference takes about 5.7 GB once unpacked and the container about 0.5 GB; several runs can share both through `REF_DIR`
and `GBRS_SIF`. Temporary BAM and EMASE files go to node-local scratch (`TMPDIR`) and are removed.

## Ceres

- Only `run.sh` runs on the login node; every step that reads data or computes is a Slurm job. Apptainer does not
  run on Ceres login nodes, which is why `./run.sh setup` is a job.
- `/90daydata` removes files not read for 90 days. Copy accepted results to `/project` (or Box) before then.
- Many short ssh sessions in a row can make the login nodes refuse new connections for a while; reuse one
  connection (OpenSSH `ControlMaster`) when scripting.

## Other Slurm clusters

Copy `examples/other_slurm.sh`: set the partition and account, and either the module names or `SOFTWARE_MODE=path`
with `apptainer` (or `gbrs`, `emase`, `bowtie`), `samtools` and `python3` on `PATH`. The container brings GBRS,
EMASE, Bowtie 1 and the Python packages. Clusters without Slurm are not supported.
