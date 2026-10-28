# Checking a GBRS run

`./run.sh validate` (queued by `./run.sh all` after `combine`) runs `scripts/validate_gbrs.py` for each mate rule.
It writes `WORKDIR/validation/<rule>/summary.json` and one table per check, and the job fails when any check fails.
`./run.sh status` prints whether each rule passed. Read `summary.json` first: it lists every failure.

## Checks

| Check | Compares | Passes when | Table |
| --- | --- | --- | --- |
| finished | the sample sheet with the sample folders | every sample has `done.json` and this rule's outputs | (in `summary.json`) |
| sex | Xist with four Y-chromosome genes (Ddx3y, Kdm5d, Eif2s3y, Uty), in counts per million | the call matches the sex used for the transition probabilities | `sex.tsv` |
| duplicates | every pair of samples, on their reconstructed genotypes | no two samples from different genotype groups look like one mouse, and samples sharing a group do | `duplicates.tsv` |
| dna | each sample with every mouse's DNA founder probabilities (`DNA_ALLELEPROBS`) | the best match is the sample's own mouse (`mouse_id`) | `dna_identity.tsv` |
| star | each sample's GBRS counts with every library in a STAR count matrix (`STAR_COUNTS`) | the best match is the sample's own library (`star_library`) | `star_identity.tsv` |

`dna` and `star` run only when their file is set; `summary.json` then says `skipped`. Samples whose mouse or library
is missing from the file are listed as `unchecked`, not failed.

### Genotype overlap

The duplicates and DNA checks score two sets of founder probabilities by their overlap: at each autosomal grid
point, the sum over the eight founders of the smaller of the two probabilities, averaged over points (every tenth
grid point by default). 1 means identical. Values from a Bennett lab DO liver run of October 2026:

| Pair | Overlap |
| --- | --- |
| A library's RNA genotype and its own mouse's DNA genotype (four pilot libraries) | 0.81 to 0.83 |
| A library's RNA genotype and the closest other mouse's DNA genotype | 0.52 to 0.59 |
| RNA genotypes of two different DO mice | about 0.28 to 0.30 |
| Two libraries holding the same reads | 1.0 |

RNA and DNA genotypes never agree fully: GBRS sees only expressed genes, and the two use different markers. Two
libraries are counted as one mouse at `SAME_MOUSE_SCORE` (0.9). Siblings and inbred mice score higher than
unrelated DO mice; give mice expected to share a genotype (one CC strain, or one mouse sequenced twice) the same
`genotype_group` in the sample list.

### Sex

In DO liver, males had 4 to 7 Xist counts per million and about 200 for the four Y genes; females about 1,600
Xist and none for the Y genes. A library with at least 100 Xist and 20 Y-gene counts per million is called `mixed`,
which suggests reads from two mice. Other tissues express Xist at different levels; read the numbers, not only the
call.

### STAR identity

Log counts per million of the 2,000 genes that vary most among the STAR libraries, centred on each gene's mean
(on each side's own mean when the run has at least ten samples, which removes the genes' overall difference between
the two pipelines). In four DO pilot libraries, the own library scored r = 0.75 to 0.82 and the next
library at most 0.53.

## When a check fails

- A sample that matches another mouse or library: do not rename or relabel anything to make the check pass. Record
  the result, check for the same finding in the STAR counts and the DNA, and ask the lab. A library can hold another
  mouse's RNA from the start (a tube or plate mix-up), or a processing step can give it another sample's reads (for
  example, joining lane files by a name prefix, so that `MOUSE_24` picks up `MOUSE_240`'s files).
- Two samples with one genotype: compare their read files. Identical read counts mean the same reads were used twice.
- Sex disagreement: check the study's records for that mouse and the DNA check. A wrong sex in the sample list
  changes the X-chromosome reconstruction; fix the list and rerun that sample (`./run.sh submit --array=N`).
- An unfinished sample: read `slurmout/gbrs_<job>_<task>.out` and the `<sample>.failed_<time>/` folder, fix the
  cause and rerun only that row.

## Comparing with earlier results

When the libraries were processed before, or when a setting or GBRS version changes, compare sample by sample:

```bash
./run.sh compare --rule=either --name=dec2025 \
    --other='/path/to/earlier/{sample}/{sample}_gbrs.quantified.diploid.genes.expected_read_counts'
```

The job writes `comparison/<rule>_vs_<name>.tsv` (per sample: total reads on each side, correlation of log
counts and of founder counts, share of genes with the same diplotype) and a `.json` summary listing samples missing
on either side and samples that are identical. Report the agreement as measured; there is no pass mark. In
a Bennett lab DO study, rerunning the 2025 method from the same reads gave identical results; `both` against `either` gave
founder-count correlations near 0.975 and the same diplotype for about 91% of genes.

## Other things to look at

- `combined/<rule>/qc.tsv`: reads and aligned reads per mate. About 80% of reads aligned in the DO pilot.
- The genome plot (`<sample>.plotted.genome.pdf`) should show two founders at most positions, changing in blocks.
