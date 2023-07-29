# Load the mia R package
library(mia)


### Defining FILE PATHS ###

## assays
url_to_microbiome_assay <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_microbiome_assay.csv")
url_to_metabolites_assay <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_metabolites_assay.csv")
url_to_biomarkers_assay <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_biomarkers_assay.csv")

## rowdata
url_to_microbiome_rowdata <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_microbiome_rowdata.csv")
url_to_metabolites_rowdata <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_metabolites_rowdata.csv")
url_to_biomarkers_rowdata <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_biomarkers_rowdata.csv")

## samples
url_to_sample_data <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_sample_data.csv")
url_to_sample_map <- url("https://github.com/microbiome/data/raw/main/HintikkaXOData/hintikka_sample_map.csv")

### Constructing TSEs ###

## microbiome experiment

microbiome_assay <- read.csv(url_to_microbiome_assay, row.names = 1)
microbiome_rowdata <- read.csv(url_to_microbiome_rowdata, row.names = 1)

microbiome_tse <- TreeSummarizedExperiment(assays = list(counts = microbiome_assay),
                                           rowData = microbiome_rowdata)

## metabolites experiment

metabolites_assay <- read.csv(url_to_metabolites_assay, row.names = 1)
metabolites_rowdata <- read.csv(url_to_metabolites_rowdata, row.names = 1, header = FALSE)
metabolites_rowdata$V2 <- NULL

metabolites_tse <- TreeSummarizedExperiment(assays = list(nmr = metabolites_assay),
                                            rowData = metabolites_rowdata)

## biomarkers experiment

biomarkers_assay <- read.csv(url_to_biomarkers_assay, row.names = 1)
biomarkers_rowdata <- read.csv(url_to_biomarkers_rowdata, row.names = 1, header = FALSE)
biomarkers_rowdata$V2 <- NULL

biomarkers_tse <- TreeSummarizedExperiment(assays = list(signals = biomarkers_assay),
                                           rowData = biomarkers_rowdata)


### Constructing MAE ###

sample_data <- read.csv(url_to_sample_data, row.names = 1)
sample_map <- read.csv(url_to_sample_map)

mae <- MultiAssayExperiment(experiments = ExperimentList(microbiota = microbiome_tse,
                                                         metabolites = metabolites_tse,
                                                         biomarkers = biomarkers_tse),
                            colData = sample_data,
                            sampleMap = sample_map)


### SAVE into RDS ###
saveRDS(mae, file = "mae.rds")
