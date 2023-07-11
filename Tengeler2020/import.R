# Load the mia R package
library(mia)


### Defining FILE PATHS ###

## Biom file (taxonomic profiles)
url_to_biom <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Aggregated_humanization2.biom")
## Sample phenodata
url_to_sample_meta <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Mapping_file_ADHD_aggregated.csv")
## Phylogenetic tree
url_to_tree <- url("https://github.com/microbiome/data/raw/main/Tengeler2020/Data_humanization_phylo_aggregation.tre")


### TSE from BIOM ###

makeTreeSEFromBiom(url)

# Import the data into SummarizedExperiment container
se <- loadFromBiom(url_to_biom)

# Convert this data to TreeSE container (no direct importer exists)
tse <- as(se, "TreeSummarizedExperiment")

# We notice that the rowData fields do not have descriptibve names.
# Hence, let us rename the columns in rowData
names(rowData(tse)) <- c("Kingdom", "Phylum", "Class",
                         "Order", "Family", "Genus")

# We also notice that the taxa names are of form "c__Bacteroidia" etc.
# Goes through the whole DataFrame. Removes '.*[kpcofg]__' from strings, where [kpcofg] 
# is any character from listed ones, and .* any character.
rowdata_modified <- BiocParallel::bplapply(rowData(tse), 
                                           FUN = stringr::str_remove, 
                                           pattern = '.*[kpcofg]__')

# Genus level has additional '\"', so let's delete that also
rowdata_modified <- BiocParallel::bplapply(rowdata_modified, 
                                           FUN = stringr::str_remove, 
                                           pattern = '\"')

# rowdata_modified is a list, so convert this back to DataFrame format. 
# and assign the cleaned data back to the TSE rowData
rowData(tse) <- DataFrame(rowdata_modified)


### SAMPLE PHENODATA ###

# It seems like a comma separated file and it does not include headers
# Let us read the file; note that sample names are in the first column
sample_meta <- read.table(url_to_sample_meta, sep=",", header=FALSE, row.names=1)

# Add headers for the columns (as they seem to be missing)
colnames(sample_meta) <- c("patient_status", "cohort",
                           "patient_status_vs_cohort", "sample_name")

# Add this sample data to colData of the taxonomic data object
# Note that the data must be given in a DataFrame format (required for our purposes)
colData(tse) <- DataFrame(sample_meta)


### PHYLOGENETIC TREE ###

# Read the tree file
tree <- ape::read.tree(url_to_tree)

# Add tree to rowTree
rowTree(tse) <- tree

### SAVE into RDS ###
saveRDS(tse, file="tse.rds")

