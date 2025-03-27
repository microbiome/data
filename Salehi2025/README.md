27.3.2025

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

-   `report.qmd` – Quarto document that loads precomputed model results and generates Figure 5.

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

    -   `tse_AMRdemo.rds`
    -   `data.R`
    -   `funcs.R`
    -   `model.R`
    -   `report.qmd`

2.  Open R or RStudio in the same directory where these files are located.

3.  Open the file `report.qmd` and follow the instructions inside to run the analysis and generate the figure.
