# Defaults for the GBRS workflow. Copy examples/ceres.sh to run_settings.sh and
# edit that file for each run; its values override these defaults.
# Then: ./run.sh check, ./run.sh references, ./run.sh sheet, ./run.sh submit.


# ======================= EDIT FROM HERE =======================

# Where output goes. Each sample keeps about 2 GB of results plus a compressed
# alignment file of 1 to 3 GB. /90daydata is purged at 90 days by access time.
WORKDIR="/90daydata/CHANGEME/gbrs_run"

# Preprocessed reads from RNAseq_STAR's HTStream stage:
#   READS_DIR/<sample>/<sample>_R1.fastq.gz and _R2.fastq.gz
READS_DIR="/path/to/run/01-HTS_Preproc"

# Tab-separated list with the columns sample, sex (M/F) and generation.
SAMPLE_LIST="/path/to/gbrs_samples.tsv"

# GBRS reference files (Zenodo record 15091310, GRCm39 / Ensembl 105).
# Point several runs at the same directory to download it once (about 2.2 GB).
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

# Software. GBRS and EMASE run from a container pinned by digest (Apptainer),
# or from PATH when GBRS_SIF is empty. Bowtie 1 and SAMtools come from site
# modules (SOFTWARE_MODE=modules) or PATH (SOFTWARE_MODE=path).
GBRS_IMAGE="docker://quay.io/jaxcompsci/gbrs_py3@sha256:916de82f093a30d6132cb351412861a6c539cde950372c0f50108a9f0257af0a"
GBRS_SIF="${WORKDIR}/containers/gbrs_py3_v1.0.1.sif"
SOFTWARE_MODE="modules"
MODULE_APPTAINER="apptainer/1.4.1"
MODULE_BOWTIE="bowtie/1.3.1"
MODULE_SAMTOOLS="samtools/1.17"
MODULE_PYTHON="python_3/3.11.1"

# SLURM.
SLURM_ACCOUNT="CHANGEME"
SLURM_PARTITION="ceres"
CPUS=16
MEMORY="64G"
TIME_LIMIT="24:00:00"
THROTTLE=40          # samples running at once

# ======================= TO HERE ==============================

settings="${PIPELINE_CONFIG:-$(dirname "${BASH_SOURCE[0]}")/run_settings.sh}"
if [[ -f "$settings" ]]; then
    default_ref_dir="$REF_DIR"; default_sif="$GBRS_SIF"
    source "$settings" || return $?
    if [[ "$REF_DIR" == "$default_ref_dir" ]]; then REF_DIR="${WORKDIR}/references"; fi
    if [[ "$GBRS_SIF" == "$default_sif" ]]; then GBRS_SIF="${WORKDIR}/containers/gbrs_py3_v1.0.1.sif"; fi
fi

SAMPLE_SHEET="${WORKDIR}/samples.tsv"
GBRS_OUT="${WORKDIR}/samples"
BOWTIE_INDEX="${REF_DIR}/bowtie/bowtie.transcripts"

load_software() {
    if [[ "$SOFTWARE_MODE" == modules ]]; then
        module load "$@" || return $?
    fi
}

# gbrs/emase through the container when one is set, else from PATH.
gbrs_tool() {
    if [[ -n "$GBRS_SIF" ]]; then
        apptainer exec --cleanenv --bind "$WORKDIR,$REF_DIR,$READS_DIR,${TMPDIR:-/tmp}" "$GBRS_SIF" "$@"
    else
        "$@"
    fi
}
