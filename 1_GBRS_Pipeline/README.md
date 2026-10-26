# 1_GBRS_Pipeline

`gbrs_sample.slurm` runs GBRS for one sample per Slurm array task. Start it with `./run.sh submit` from the
repository root after `./run.sh sheet` has written `WORKDIR/samples.tsv` (see the main README).

## Steps in each task

| Step | Command | Output kept in `WORKDIR/samples/<sample>/` |
| --- | --- | --- |
| Align R1 and R2 separately to the pooled-founder transcriptome | `bowtie -q -a --best --strata --sam -v 3` | `<sample>.bowtie_R1.log`, `_R2.log` |
| Convert each BAM to an EMASE alignment file | `gbrs bam2emase -h A,B,C,D,E,F,G,H` | (temporary) |
| Keep alignments common to both mates | `emase get-common-alignments` | (temporary) |
| Compress | `gbrs compress` | `<sample>.compressed.emase.h5` |
| Quantify (multiway) | `gbrs quantify -M 4 --report-alignment-counts` | `<sample>.multiway.*` |
| Reconstruct the genome | `gbrs reconstruct` with `tranprob.DO.G<generation>.<sex>.npz` | `<sample>.genotypes.tsv`, `.genoprobs.npz`, `.genotypes.npz` |
| Quantify with the reconstructed genotypes | `gbrs quantify -G` | `<sample>.diploid.*` (allele-specific expression) |
| Interpolate, export and plot | `gbrs interpolate`, `export`, `plot` | `<sample>.interpolated.genoprobs.npz`, `.gbrs.interpolated.genoprobs.tsv`, `.plotted.genome.pdf` |

BAMs and uncompressed EMASE files are written to node-local scratch and removed. A task writes into
`<sample>.partial/` and renames it to `<sample>/` only after every output exists and `done.json` is written; a
failed attempt is kept as `<sample>.failed_<time>/`. A finished sample is never overwritten.

## Settings that change results

`EMASE_MODEL` (4), `HAPLOTYPES`, `TPROB_TEMPLATE` with each sample's sex and generation, the reference directory
(recorded with checksums in `reference_manifest.json`) and the container image (pinned by digest). `done.json`
records all of them for each sample.

## Resources

On Ceres, 16 CPUs and 96 GB per sample, typically several hours per sample; alignment of the two mates takes most
of the time. `THROTTLE` limits how many samples run at once.
