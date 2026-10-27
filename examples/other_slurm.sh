# Another Slurm cluster. Copy to run_settings.sh and replace the paths and names.
WORKDIR="/scratch/YOUR_LAB/gbrs_run"
SAMPLE_LIST="/scratch/YOUR_LAB/gbrs_samples.tsv"   # sample, sex, generation, r1, r2 (paths to the read files)
SLURM_ACCOUNT="YOUR_ACCOUNT"
SLURM_PARTITION="YOUR_PARTITION"
# Either the cluster's module names (check with: module avail apptainer samtools python) ...
SOFTWARE_MODE="modules"
MODULE_APPTAINER="apptainer"
MODULE_SAMTOOLS="samtools"
MODULE_PYTHON="python"
# ... or SOFTWARE_MODE="path" with apptainer, samtools and python3 already on PATH.
# Without Apptainer, set GBRS_SIF="" and install gbrs 1.1.0 (with emase), bowtie 1 and samtools on PATH.
CPUS=16
MEMORY="64G"
TIME_LIMIT="12:00:00"
THROTTLE=20
