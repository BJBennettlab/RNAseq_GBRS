# 4_CCRIX_Implementation

GBRS implementation and validation for CC-RIX (CCRIX) mice — reciprocal intercrosses between inbred Collaborative Cross (CC) strains. This directory documents the approach, validation workflow, and results for extending GBRS-based genotype reconstruction and allele-specific expression (ASE) quantification to CCRIX samples.

## Overview

Generating CCRIX-specific transition probability files proved intractable. Following advice from Dan Gatti (Churchill Lab, Jackson Laboratory), **DO.G0 transition probability files** (from the GBRS reference data Zenodo archive; Lloyd 2025) are used instead for GBRS steps 5 and 6 (genotype reconstruction and ASE quantification). Validation via cosine similarity against JAX consensus microarray diplotype probabilities confirmed this approach is highly concordant.

**Key finding:** The DO.G0 transprob file (most conservative — permits least recombination) performed best among all tested generations (G0, G8, G9, G10, G100), though generation had minimal overall effect on concordance.

## Contents

| File | Description |
|------|-------------|
| `CCRIX_GBRS_Testing.Rmd` | Main validation notebook — cosine similarity comparisons of GBRS-derived CCRIX diplotypes against JAX microarray-based diplotype probabilities |
| `CCDiet_genotypes_founderProb.Rmd` | Secondary validation — converts CCDiet project genotypes (CC001, CC019, CC040) into founder probability matrices and cross-validates against GBRS reconstruction results |

## CCRIX Crosses

Six reciprocal intercrosses from three CC parental strains:

| Cross |
|-------|
| CC001 × CC019 |
| CC019 × CC001 |
| CC019 × CC040 |
| CC040 × CC019 |
| CC040 × CC001 |
| CC001 × CC040 |

## Validation Workflow

Full analysis is in `CCRIX_GBRS_Testing.Rmd`. The workflow consists of four steps:

### 1. Interpolate Microarray Markers to Genes

The JAX consensus CC diplotype file (`CC_diplotype_probs.rds`) is defined by microarray markers (143,259 markers; `gm_uwisc_v4.csv` annotation). GBRS uses GRCm39 v105 Ensembl gene annotations (55,414 genes), so markers must be mapped to gene space before comparison. Interpolation code provided by Dan Gatti.

- Input: `CC_diplotype_probs.rds` + `gm_uwisc_v4.csv`
- Output: `CC_gene_interpolated_probs` — 53,751 genes with interpolated diplotype probabilities

### 2. CC Baseline Concordance

Before evaluating CCRIX samples, establish the baseline concordance for the three parental CC strains (CC001, CC019, CC040). Cosine similarity between GBRS-derived CC diplotypes and `CC_gene_interpolated_probs`.

- **Result: ~95% concordance for CC samples**

### 3. Generate CCRIX Combined Parental Matrices

CCRIX mice are F1 offspring of two inbred CC strains, so their expected diplotype cannot be compared directly to either parent. Both parents' `CC_gene_interpolated_probs` are averaged at each gene to produce combined matrices. Male CCRIX pups are hemizygous on the X chromosome (inherited from mother), which is accounted for in the final **"Combined matrices by Sex"**.

### 4. CCRIX Concordance Evaluation

Cosine similarity of GBRS-derived CCRIX diplotypes vs. "Combined matrices by Sex". Tested multiple DO transprob generations — DO.G0 consistently performed best.

- **Result: ~90% concordance for CCRIX samples**

## Recommended Pipeline Usage

For CCRIX samples, use the DO.G0 transition probability files from the GBRS Zenodo reference archive (Lloyd 2025):
keep the default `TPROB_TEMPLATE` (`transition_probabilities/tranprob.DO.G{generation}.{sex}.npz`) and give every
CCRIX sample `generation` 0 in the sample list.

## Large Output Files (Not in Git)

Quantification outputs are stored on AFRI storage:

| Description | Path |
|-------------|------|
| CCRIX GBRS quantification (DO.G0) — best results | `GBRS_development/Development/quantify_ASE_G0/` |
| Final aggregated CCRIX counts | `GBRS_development/Development/Final_counts_CCRIX/` |
| Interpolated genotype probability TSVs (176 samples) | `GBRS_development/../genotype_probs/` |

## Collaborators

- **Mike Lloyd** — Churchill Lab, Jackson Laboratory (GBRS development, Zenodo reference data)
- **Dan Gatti** — Churchill Lab, Jackson Laboratory (interpolation code, transprob strategy advice)

## Reference

> Lloyd, Michael. (2025). *Reference Data for EMASE and GBRS*. Zenodo. https://doi.org/10.5281/zenodo.15091310
