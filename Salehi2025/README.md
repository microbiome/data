19.3.2025


## Using the TreeSummarizedExperiment R data object

Just download tse_AMRdemo.rds

This contains the data in a readily processed format.

# Demo Dataset for TreeSummarizedExperiment

## Study Information
- Title: Gender differences in global antimicrobial resistance
- Authors: Mahkameh Salehi, Ville Laitinen, Shivang Bhanushali, Johan Bengtsson-Palme, Peter Collignon, John J Beggs, Katariina Pärnänen, Leo Lahti
- Description: This dataset contains antimicrobial resistance (AMR) data in a TreeSummarizedExperiment (TSE) format. The dataset has been processed to facilitate further analyses of AMR patterns, particularly in relation to gender differences on a global scale.

## Contents
- `tse_AMRdemo.rds` – The full TreeSummarizedExperiment object in R’s native `.rds` format.

## How to Use
To load the dataset in R:
```r
# Load the TSE object
tse <- readRDS("tse_AMRdemo.rds")

# Check the structure of the dataset
tse
