## Enrichment analysis: High ARG prevalence across age categories #######################

incomes <- c("LMIC", "HIC")
fit_enrichment_age <- list()

for (ic in incomes) {
  
  message("Fitting enrichment model for: ", ic)
  
  # Filter dataset for current income group
  my_data <- enrichment_df %>%
    filter(income_group == ic)
  
  # Fit Bayesian logistic regression with age category as predictor
  fit <- brm(
    high_ARG ~ age_category,
    data = my_data,
    family = bernoulli(),
    prior = set_prior("normal(0, 1)", class = "b"),
    chains = 2,
    iter = 4000,
    control = list(adapt_delta = 0.99),
    cores = parallel::detectCores()
  )
  
  fit_enrichment_age[[ic]] <- fit
}

# Save model
saveRDS(fit_enrichment_age, "fit_enrichment_age.rds")



## Enrichment analysis: High ARG prevalence across gender categories #######################

incomes <- c("LMIC", "HIC")
fit_enrichment_gender <- list()

for (ic in incomes) {
  
  message("Fitting enrichment model for: ", ic)
  
  # Filter dataset for current income group
  my_data <- enrichment_df %>%
    filter(income_group == ic)
  
  # Fit Bayesian logistic regression with gender as predictor
  fit <- brm(
    high_ARG ~ gender,
    data = my_data,
    family = bernoulli(),
    prior = set_prior("normal(0, 1)", class = "b"),
    chains = 2,
    iter = 4000,
    control = list(adapt_delta = 0.99),
    cores = parallel::detectCores()
  )
  
  fit_enrichment_gender[[ic]] <- fit
}

# Save model
saveRDS(fit_enrichment_gender, "fit_enrichment_gender.rds")





## Enrichment analysis: High ARG prevalence across continent categories #######################

fit_enrichment_continent <- list()

fit <- brm(
  high_ARG ~ continent,
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 4000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_continent[["all"]] <- fit

# Save model
saveRDS(fit_enrichment_continent, "fit_enrichment_continent.rds")
