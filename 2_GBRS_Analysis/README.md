# 2_GBRS_Analysis

R Markdown script for allele-specific expression (ASE) visualization and interpretation of GBRS pipeline output.

## Scripts

### `GBRS_ASE.Rmd`
Visualizes and analyzes allele-specific expression results from `1_GBRS_Pipeline/`. Includes strain-level expression summaries and ASE plots.

## Usage

```r
rmarkdown::render("GBRS_ASE.Rmd")
```

## Notes

- Transition probability matrices and CC/CCRIX implementation files required for genotype reconstruction are in `3_CC_Implementation/`
- Input count data is each sample's `<sample>.diploid.genes.expected_read_counts` from `1_GBRS_Pipeline/gbrs_sample.slurm`, joined by `./run.sh combine`
