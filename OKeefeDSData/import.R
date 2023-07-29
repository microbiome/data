# Load the mia R package
library(mia)


### Defining FILE PATHS ###

url_to_assay <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_assay.csv")
url_to_coldata <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_coldata.csv")
url_to_rowdata <- url("https://github.com/microbiome/data/raw/main/OKeefeDSData/okeefe_rowdata.csv")

### Importing FILES ###

okeefe_assay <- read.csv(url_to_assay)
okeefe_coldata <- read.csv(url_to_coldata)
okeefe_rowdata <- read.csv(url_to_rowdata)

### Constructing TSE ###

tse <- TreeSummarizedExperiment(assays = list(counts = okeefe_assay),
                                colData = okeefe_coldata,
                                rowData = okeefe_rowdata)