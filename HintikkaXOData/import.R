# Load the mia R package
library(mia)


### Defining FILE PATHS ###

## assays
url_to_microbiome_assay <- url()
url_to_metabolites_assay <- url()
url_to_biomarkers_assay <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Aggregated_humanization2.biom")

## coldata
url_to_microbiome_coldata <- url()
url_to_metabolites_coldata <- url()
url_to_biomarkers_coldata <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Mapping_file_ADHD_aggregated.csv")

## rowdata
url_to_microbiome_rowdata <- url()
url_to_metabolites_rowdata <- url()
url_to_biomarkers_rowdata <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Data_humanization_phylo_aggregation.tre")
