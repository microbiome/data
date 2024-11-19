# Load the mia R package
library(mia)

### Defining FILE PATHS ###

url_phyloseq_obj <- url("https://github.com/microbiome/data/raw/main/MovingPictures/MovingPictures.rda")

### TSE from Phyloseq ###

# Read phyloseq from url
phyloseq_obj <- load(url_phyloseq_obj)

# Convert phyloseq to TreeSE
tse <- convertFromPhyloseq(phyloseq_obj)

### SAVE into RDS ###
saveRDS(tse, file="MovingPictures.rds")
