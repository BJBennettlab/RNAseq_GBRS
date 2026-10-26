# Ceres (SCINet) example. Copy to run_settings.sh and replace the paths.
WORKDIR="/90daydata/YOUR_PROJECT/gbrs_run"
READS_DIR="/90daydata/YOUR_PROJECT/rnaseq_run/01-HTS_Preproc"
SAMPLE_LIST="/90daydata/YOUR_PROJECT/gbrs_samples.tsv"   # columns: sample, sex (M/F), generation
# Share one verified reference between runs:
# REF_DIR="/90daydata/YOUR_PROJECT/references/gbrs_zenodo_15091310"
SLURM_ACCOUNT="YOUR_PROJECT"
SLURM_PARTITION="ceres"
CPUS=16
MEMORY="96G"
TIME_LIMIT="2-00:00:00"
THROTTLE=40
