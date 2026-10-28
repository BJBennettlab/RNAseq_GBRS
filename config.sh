# Defaults for the GBRS workflow. Copy examples/ceres.sh to run_settings.sh and
# edit that file for each run; its values override these defaults.
# Then: ./run.sh sheet, ./run.sh check, ./run.sh all.


# ======================= EDIT FROM HERE =======================

# Where output goes. Each sample keeps about 200 MB per mate rule (measured on
# 20-million-pair libraries). /90daydata is purged at 90 days by access time.
WORKDIR="/90daydata/CHANGEME/gbrs_run"

# Preprocessed reads from RNAseq_STAR's HTStream stage:
#   READS_DIR/<sample>/<sample>_R1.fastq.gz and _R2.fastq.gz
# Not used when the sample list has r1 and r2 columns.
READS_DIR="/path/to/run/01-HTS_Preproc"

# Tab-separated list. Required columns: sample, sex (M/F), generation.
# Optional: r1 and r2 (read files, instead of READS_DIR), mouse_id (ID in the DNA
# file), genotype_group (samples expected to share a genotype), star_library.
SAMPLE_LIST="/path/to/gbrs_samples.tsv"

# GBRS reference files (Zenodo record 15091310, GRCm39 / Ensembl 105).
# Point several runs at the same directory to download it once (about 2.2 GB to
# download, 5.7 GB unpacked).
REF_DIR="${WORKDIR}/references"

# Transition probabilities. {sex} becomes M or F and {generation} the sample's
# generation. DO mice use the files in REF_DIR; for CC mice point this at
# 3_CC_Implementation/tranprob.CC.G20.{sex}.npz.
TPROB_TEMPLATE="transition_probabilities/tranprob.DO.G{generation}.{sex}.npz"

# Founder letters A-H: AJ, B6, 129, NOD, NZO, CAST, PWK, WSB.
HAPLOTYPES="A,B,C,D,E,F,G,H"
EMASE_MODEL=4

# How the two mates of a read pair are combined; each rule gets its own results from one alignment.
#   either  a pair counts toward a transcript and founder where either mate aligned (the lab's 2025 GBRS rule)
#   both    only where both mates aligned (emase get-common-alignments, the GBRS authors' paired-end step)
MATE_RULES="either both"

# Checks after the run (./run.sh validate; docs/validation.md). Both files are
# optional. DNA_ALLELEPROBS: founder probabilities from the mice's DNA genotypes
# (scripts/export_dna_alleleprobs.R writes it from R/qtl2). STAR_COUNTS: a
# gene-by-library count matrix from RNAseq_STAR for the same libraries.
DNA_ALLELEPROBS=""
STAR_COUNTS=""
SAME_MOUSE_SCORE=0.9   # genotype overlap at which two samples count as one mouse

# Software. GBRS and EMASE run from a container pinned by digest (Apptainer),
# or from PATH when GBRS_SIF is empty. The default is GBRS 1.1.0 (tag
# v1.1.0-338c782), the version JAX's pipeline uses. It gave byte-identical
# results to GBRS 1.0.0 (tag v1.0.1, the lab's 2025 version) on four DO_Chln
# libraries; a study already run with 1.0.0 keeps it (docs/method.md, Versions).
# The container also brings Bowtie 1. SAMtools and Python come from site modules
# (SOFTWARE_MODE=modules) or PATH (SOFTWARE_MODE=path).
# The lab keeps copies with the same digests at ghcr.io/bjbennettlab/gbrs (1.1.0 and 1.0.1); to use one, replace
# quay.io/jaxcompsci/gbrs_py3 below with ghcr.io/bjbennettlab/gbrs and keep the digest. On Ceres the pulled images
# are also in /project/do2_projects/containers (gbrs_1.1.0.sif, gbrs_1.0.1.sif).
GBRS_IMAGE="docker://quay.io/jaxcompsci/gbrs_py3@sha256:95ac892ffb209ee7a39846dfd62203b8738d440d436c4202f64caa9d8ca329b0"
GBRS_SIF="${WORKDIR}/containers/gbrs_py3_v1.1.0-338c782.sif"
SOFTWARE_MODE="modules"
MODULE_APPTAINER="apptainer/1.4.1"
MODULE_BOWTIE="bowtie/1.3.1"
MODULE_SAMTOOLS="samtools/1.17"
MODULE_PYTHON="python_3/3.11.1"

# SLURM.
SLURM_ACCOUNT="CHANGEME"
SLURM_PARTITION="ceres"
CPUS=16
MEMORY="64G"         # measured peak 32 GB per sample
TIME_LIMIT="12:00:00" # measured 1.5 to 2.6 hours for 20-million-pair libraries
THROTTLE=40          # samples running at once

# ======================= TO HERE ==============================

settings="${PIPELINE_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/run_settings.sh}"
if [[ -f "$settings" ]]; then
    default_ref_dir="$REF_DIR"; default_sif="$GBRS_SIF"
    source "$settings" || return $?
    if [[ "$REF_DIR" == "$default_ref_dir" ]]; then REF_DIR="${WORKDIR}/references"; fi
    if [[ "$GBRS_SIF" == "$default_sif" ]]; then GBRS_SIF="${WORKDIR}/containers/gbrs_py3_v1.1.0-338c782.sif"; fi
fi

SAMPLE_SHEET="${WORKDIR}/samples.tsv"
GBRS_OUT="${WORKDIR}/samples"
BOWTIE_INDEX="${REF_DIR}/bowtie/bowtie.transcripts"

load_software() {
    if [[ "$SOFTWARE_MODE" == modules ]]; then
        module load "$@" || return $?
    fi
}
