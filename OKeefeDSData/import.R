# Load the mia R package
library(mia)


### Defining FILE PATHS ###

url_to_assay <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_assay.csv")
url_to_coldata <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_coldata.csv")
url_to_rowdata <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_rowdata.csv")


### Importing FILES ###

okeefe_assay <- read.csv(url_to_assay, row.names = 1)
colnames(okeefe_assay) <- gsub("\\.", "-", colnames(okeefe_assay))

okeefe_coldata <- read.csv(url_to_coldata, row.names = 1)
okeefe_rowdata <- read.csv(url_to_rowdata, row.names = 1)


### Constructing TSE ###

tse <- TreeSummarizedExperiment(assays = list(counts = okeefe_assay),
                                colData = okeefe_coldata,
                                rowData = okeefe_rowdata)


### SAVE into RDS ###
saveRDS(tse, file = "tse.rds")
