19.3.2025


## Using the TreeSummarizedExperiment R data object

Just download tse_AMRdemo.rds

This contains the data in a readily processed format.

# Demo Dataset for TreeSummarizedExperiment

## Study Information
- Title: Gender differences in global antimicrobial resistance
- Authors: Mahkameh Salehi, Ville Laitinen, Shivang Bhanushali, Johan Bengtsson-Palme, Peter Collignon, John J Beggs, Katariina Pärnänen, Leo Lahti
- Description: This dataset contains antimicrobial resistance (AMR) data in a TreeSummarizedExperiment (TSE) format. The dataset has been processed to facilitate further analyses of AMR patterns, particularly in relation to gender differences on a global scale.
- [DOI: https://doi.org/10.5281/zenodo.14909582](https://doi.org/10.5281/zenodo.14909582)

## Contents
- `tse_AMRdemo.rds` – The full TreeSummarizedExperiment object in R’s native `.rds` format.
- `Figure5.qmd` – A Quarto document that performs data processing and generates Figure 5 from the manuscript.
- `Figure5models.qmd` – A separate Quarto script for fitting the statistical models used in the figure (runs slower).

## How to Use
To load the dataset in R:
```r
# Load the TSE object
tse <- readRDS("tse_AMRdemo.rds")

# Check the structure of the dataset
tse
