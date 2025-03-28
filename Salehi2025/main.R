# Data processing
source("data.R")

# Utility functions
source("funcs.R")

# Modeling (slow, only needs to be done once)
source("model.R")

# Render the Quarto report and generate Figure 5
quarto::quarto_render("report.qmd")
