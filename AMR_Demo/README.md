18.3.2025


## Using the TreeSummarizedExperiment R data object

Just download TSE_filtered.rds

This contains the data in a readily processed format.

# Demo Dataset for TreeSummarizedExperiment

## Overview
This repository contains a cleaned and publicly shareable demo dataset.
It includes a TreeSummarizedExperiment (TSE) object representing the gender study subset of the Global AMR data.

## Contents
- `tse_AMRdemo.rds` – The full TreeSummarizedExperiment object in R’s native `.rds` format.

## How to Use
To load the dataset in R:
```r
# Load the TSE object
tse <- readRDS("tse_AMRdemo.rds")

# Check the structure of the dataset
tse
