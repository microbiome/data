# This R script includes the data processing steps required for the analysis of Figure 5.
# It prepares and cleans the input data

# Libraries
library(scater)
library(ggpubr)
library(rstatix)
library(viridis)
library(ggsignif)
library(Matrix)
library(RColorBrewer)
library(caret)
library(randomForest)
library(xgboost)
library(brms)
library(glmnet)
library(tidyverse)
library(vegan)
library(TreeSummarizedExperiment)
library(dplyr)
library(loo)
library(cmdstanr)
library(broom)
library(tibble)

# Load TSE object
tse <- readRDS("tse_AMRdemo.rds")

# Filter samples with complete gender and age data
non_na_samples <- !is.na(colData(tse)$gender) & !is.na(colData(tse)$host_age_years)
TSE_filtered <- tse[, non_na_samples]

# Extract metadata
adult_metadata <- as.data.frame(colData(TSE_filtered))

# Total read counts per sample
counts_mat <- assay(TSE_filtered, "counts")
adult_metadata$readcount_2 <- colSums(counts_mat)


# This section extracts ARG class information for each gene, reshapes the count matrix to long format, 
# and computes total counts for each antibiotic class per sample. 
# Only the top 5 most common AB classes are retained for further analysis.

# Extract count matrix and metadata
counts_mat <- assay(tse, "counts")
row_df <- as.data.frame(rowData(tse))  # Gene annotations
col_df <- as.data.frame(colData(tse)) %>% mutate(sample = acc)  # Sample metadata

# Set sample names as column names
colnames(counts_mat) <- col_df$acc 

# Merge gene-level counts with their corresponding antibiotic resistance class (long format)
counts_df <- as.data.frame(counts_mat) %>%
  mutate(gene = rownames(counts_mat)) %>%
  left_join(row_df[, c("GENE", "Class")], by = c("gene" = "GENE")) %>%
  pivot_longer(cols = -c(gene, Class), names_to = "sample", values_to = "count")

# Aggregate counts at sample and antibiotic class level
agg_counts <- counts_df %>%
  group_by(sample, Class) %>%
  summarize(class_total = sum(count, na.rm = TRUE), .groups = "drop")

# Keep only the top 5 most frequent antibiotic classes
top_5_AB_classes <- c(
  "Tetracycline", 
  "Beta-lactam", 
  "Macrolide, Lincosamide, Streptogramin B", 
  "Aminoglycoside", 
  "Amphenicol"
)

agg_counts <- agg_counts %>% filter(Class %in% top_5_AB_classes)


# This section prepares AB class count data for modeling by: 
# 1. Adding a small pseudocount to avoid log(0) 
# 2. Converting long format to wide (samples x AB classes) 
# 3. Merging with metadata 
# 4. Applying log transformation to the top 5 AB classes

# Add pseudocount to avoid log(0)
agg_counts$class_total <- ifelse(agg_counts$class_total == 0,
                                 min(agg_counts$class_total[agg_counts$class_total != 0]),
                                 agg_counts$class_total)

# Reshape data: one row per sample, columns = AB classes
agg_counts_wide <- agg_counts %>%
  spread(key = "Class", value = "class_total")

# Filter only samples that exist in metadata
agg_counts_wide <- agg_counts_wide %>%
  filter(sample %in% adult_metadata$acc)

# Merge wide-format AB class counts with metadata
adult_metadata <- full_join(adult_metadata %>% mutate(sample = acc),
                            agg_counts_wide, 
                            by = "sample")

# Log-transform the top 5 antibiotic class counts
top5_log <- adult_metadata[, top_5_AB_classes] %>% log

# Rename columns to indicate log-transformed values
colnames(top5_log) <- paste0("log_", gsub(" ", "_", colnames(top5_log)))

# Add log-transformed AB class values to metadata
adult_metadata <- cbind(adult_metadata, top5_log)


# This section prepares the data for alpha diversity analysis and 
# computes Shannon diversity index from relative abundance profiles

# Filter metadata to include only samples with GDP and antibiotic usage information
adult_metadata <- adult_metadata %>%
  mutate(
    continent = factor(geo_loc_name_country_continent_calc)
  ) %>%
  filter(!is.na(Usage), !is.na(GDP_per_head))

# Keep only samples with non-missing GDP and antibiotic usage in the TSE object
non_na_conditions <- !is.na(colData(TSE_filtered)$GDP_per_head) &
  !is.na(colData(TSE_filtered)$Usage)
TSE_adult_filtered <- TSE_filtered[, non_na_conditions]

# Update the sample metadata in the filtered TSE object
colData(TSE_adult_filtered) <- DataFrame(adult_metadata)

# Extract relative abundance data
assay_data_clean <- assay(TSE_adult_filtered, "relabundance")

# Remove samples that contain only zero abundances
non_empty_samples <- colSums(assay_data_clean > 0, na.rm = TRUE) > 0
assay_data_clean <- assay_data_clean[, non_empty_samples]

# Subset metadata to match samples with non-zero data
adult_metadata <- colData(TSE_adult_filtered)[non_empty_samples, ]

# Calculate Shannon diversity for each sample (based on species richness and evenness)
adult_metadata$shannon_diversity <- diversity(t(assay_data_clean), index = "shannon")


# This section creates dummy variables (one-hot encoding) for categorical 
# predictors to be used in the GLM and Bayesian regression models

# Copy metadata to a temporary dataframe
temp_df <- adult_metadata

# Set the correct factor level order for age categories
temp_df$age_category <- factor(temp_df$age_category,
                               levels = c("Middle-Aged Adult", "Infant",
                                          "Toddler", "Children", "Teenager",
                                          "Young Adult", "Older Adult", "Oldest Adult"))

# Create binary variable: High income (1) vs. other (0)
temp_df$income_group_HIC <- ifelse(temp_df$World_Bank_Income_Group == "High income", 1, 0)

# Create binary variable for antibiotic usage: high (1) vs. low (0)
temp_df$Usage_high <- ifelse(temp_df$Usage < 10, 0, 1)

# Convert gender to binary numeric variable: Men = 0, Women = 1
temp_df$sex_num_Men <- ifelse(temp_df$gender == "Men", 0, 1)

# Perform dummy (one-hot) encoding for selected predictors
temp_df_one_hot <- model.matrix(~ 0 + 
                                  sex_num_Men +
                                  continent + 
                                  age_category + 
                                  income_group_HIC + 
                                  Usage_high,
                                data = temp_df) %>%
  as.data.frame()

# Clean column names by replacing spaces and dashes with underscores
colnames(temp_df_one_hot) <- gsub(" |-", "_", colnames(temp_df_one_hot))

# Store the dummy variable names for later model specification
dummy_var_names <- colnames(temp_df_one_hot)

# Combine dummy variables with original metadata
adult_metadata <- cbind(adult_metadata, temp_df_one_hot)

# Clean up temporary dataframes
rm(temp_df, temp_df_one_hot)


# Create a balanced subset to make the hierarchical model faster to fit
# We sample 20 observations per group defined by income group, gender, and age category
subset_data <- as.data.frame(adult_metadata) %>%
  group_by(age_category, gender) %>%
  group_modify(~ slice_sample(.x, n = min(20, nrow(.x)))) %>%
  ungroup()

# Set "Middle-Aged Adult" as the reference category for age
subset_data$age_category <- relevel(factor(subset_data$age_category), ref = "Middle-Aged Adult")

# Create a readable categorical variable for income group (HIC vs. LMIC)
subset_data <- subset_data %>%
  mutate(income_group = if_else(income_group_HIC == 1, "HIC", "LMIC"))

# Create base data frame for modeling and set "Middle-Aged Adult" as the reference category for age
adult_df <- as.data.frame(adult_metadata) %>%
  mutate(
    age_category = factor(age_category,
                          levels = c("Infant", "Toddler", "Children", "Teenager",
                                     "Young Adult", "Middle-Aged Adult",
                                     "Older Adult", "Oldest Adult")),
    age_category = relevel(age_category, ref = "Middle-Aged Adult")
  )

# Create a readable categorical variable for income group (HIC vs. LMIC)
adult_df <- adult_df %>%
  mutate(income_group = if_else(income_group_HIC == 1, "HIC", "LMIC"))

# LMIC group: continent base level Asia
adult_df_LMIC <- adult_df %>%
  filter(income_group == "LMIC") %>%
  mutate(
    continent = relevel(factor(continent), ref = "Asia"),
    age_category = relevel(factor(age_category), ref = "Middle-Aged Adult"),
    gender = relevel(factor(gender), ref = "Men"),
    Usage_high = factor(Usage_high)
  )

# HIC group: continent base level Europe
adult_df_HIC <- adult_df %>%
  filter(income_group == "HIC") %>%
  mutate(
    continent = relevel(factor(continent), ref = "Europe"),
    age_category = relevel(factor(age_category), ref = "Middle-Aged Adult"),
    gender = relevel(factor(gender), ref = "Men"),
    Usage_high = factor(Usage_high)
  )


