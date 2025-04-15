# Data processing
source("data.R")

# Utility functions
source("funcs.R")

# Modeling (slow, only needs to be done once)
source("model.R")

# Model comparison: log-normal vs. log-transformed normal (slow, only needs to be done once)
source("ComparisonModels.R")

# Stratified model fitting (by subgroups like Income/gender) (slow, only needs to be done once)
source("StratifiedModel.R")

# Hierarchical model with nested random effects (slow, only needs to be done once)
source("HierarchicalModel.R")

# Simulation-based comparison: Bayesian vs. Frequentist under small samples (ARG + Shannon)
# (slow, only needs to be done once)
source("FrequentistCompare.R")

# Render the full Quarto report
quarto::quarto_render("report.qmd")

