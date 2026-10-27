# Ceres (SCINet) example. Copy to run_settings.sh and replace the paths.
WORKDIR="/90daydata/YOUR_PROJECT/gbrs_run"
READS_DIR="/90daydata/YOUR_PROJECT/rnaseq_run/01-HTS_Preproc"   # not needed when the sample list has r1 and r2
SAMPLE_LIST="/90daydata/YOUR_PROJECT/gbrs_samples.tsv"           # columns: sample, sex (M/F), generation
# Share one verified reference and container between runs:
# REF_DIR="/90daydata/YOUR_PROJECT/references/gbrs_zenodo_15091310"
# GBRS_SIF="/90daydata/YOUR_PROJECT/containers/gbrs_py3_v1.1.0-338c782.sif"
# Checks after the run (docs/validation.md), both optional:
# DNA_ALLELEPROBS="/90daydata/YOUR_PROJECT/dna_alleleprobs.tsv.gz"   # scripts/export_dna_alleleprobs.R
# STAR_COUNTS="/90daydata/YOUR_PROJECT/rnaseq_run/03-Counts/counts.txt"
# CC mice: TPROB_TEMPLATE="transition_probabilities/tranprob.CC.G20.{sex}.npz"
SLURM_ACCOUNT="YOUR_PROJECT"
SLURM_PARTITION="ceres"
CPUS=16
MEMORY="64G"
TIME_LIMIT="12:00:00"
THROTTLE=40
