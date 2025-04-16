28.3.2025

# Gender Differences in Global Antimicrobial Resistance

This repository contains a cleaned demo dataset (`tse_AMRdemo.rds`) and supporting code used in the analysis of gender-related differences in global antimicrobial resistance (AMR), as illustrated in Figure 5 of the manuscript.

## Study Information

-   Title: Gender differences in global antimicrobial resistance

-   Authors: Mahkameh Salehi, Ville Laitinen, Shivang Bhanushali, Johan Bengtsson-Palme, Peter Collignon, John J Beggs, Katariina Pärnänen, Leo Lahti

-   DOI: <https://doi.org/10.5281/zenodo.14909582>

## Contents

-   `tse_AMRdemo.rds` – Cleaned TreeSummarizedExperiment object containing AMR gene count and metadata.

-   `data.R` – Preprocessing script: filters metadata, calculates ARG loads, top antibiotic class abundances, and computes diversity indices.

-   `funcs.R` – Utility functions used for formatting and significance annotations.

-   `model.R` – Runs the Bayesian models for the manuscript (slow, runs brms).

-   `ComparisonModels.R` – Comparison of log-normal and log-transformed normal models.

-   `StratifiedModel.R` – Stratified model fitting by subgroups such as age and gender.

-   `HierarchicalModel.R` – Hierarchical model with nested random effects.

-   `FrequentistCompare.R` – Simulation-based comparison of Bayesian and frequentist models for ARG and Shannon.

-   `report.qmd` – Quarto document that loads precomputed model results and generates Figure 5.

-    `main.R` – Wrapper script to execute all of the above in correct order.

## How to Use and Reproduce Figure 5

To load the dataset in R:

``` r
# Load the TSE object
tse <- readRDS("tse_AMRdemo.rds")

# Check the structure of the dataset
tse
```

To reproduce the analysis and generate Figure 5 from the manuscript:

1.  [Download all necessary files](https://github.com/microbiome/data/blob/main/Salehi2025) from the repository:

    - `tse_AMRdemo.rds`
    - `data.R`
    - `funcs.R`
    - `model.R`
    - `ComparisonModels.R`
    - `StratifiedModel.R`
    - `HierarchicalModel.R`
    - `FrequentistCompare.R`
    - `report.qmd`
    - `main.R`

2.  Open R or RStudio in the same directory where these files are located.

3. Run the following line in your R console to execute all required preprocessing and modeling steps:

```r
source("main.R")
```
