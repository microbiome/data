## Helpers
# This section defines utility functions used later in the analysis 
# for formatting results and interpreting significance.

# Format numbers: Convert values below 0.001 to "<0.001"
to_neat <- function(x, threshold = 0.001) {
  ifelse(x > threshold, as.character(x), "<0.001")
}

# Assign significance stars for confidence intervals (log-scale)
get_stars_CI_log <- function(lower, upper) {
  ifelse(lower > 1 | upper < 1, "*", "ns")
}

# Assign significance stars for confidence intervals (linear scale)
get_stars_CI <- function(lower, upper) {
  ifelse(lower > 0 | upper < 0, "*", "ns")
}
