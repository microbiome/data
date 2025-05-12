## Enrichment analysis: High ARG prevalence across age categories

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



## Enrichment analysis: High ARG prevalence across gender categories 

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





## Enrichment analysis: High ARG prevalence across continent categories 

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





## Enrichment analysis: High ARG prevalence across age categories (hierarchical model)

fit_enrichment_age_hier <- list()

# Fit hierarchical Bayesian logistic regression with random slopes for income_group
fit <- brm(
  high_ARG ~ age_category + (age_category | income_group),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 2000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_age_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_age_hier, "fit_enrichment_age_hier.rds")






## Enrichment analysis: High ARG prevalence across gender categories (hierarchical model)

fit_enrichment_gender_hier <- list()

# Fit hierarchical Bayesian logistic regression with random slopes for income_group
fit <- brm(
  high_ARG ~ gender + (gender | income_group),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 2000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_gender_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_gender_hier, "fit_enrichment_gender_hier.rds")




## Enrichment analysis: High ARG prevalence across age and gender categories (hierarchical model)

fit_enrichment_age_gender_hier <- list()

# Fit hierarchical Bayesian logistic regression with random intercepts and slopes for age_category within income_group and gender
fit <- brm(
  high_ARG ~ age_category + (age_category | income_group/gender),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 2000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_age_gender_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_age_gender_hier, "fit_enrichment_age_gender_hier.rds")



## Enrichment analysis: High ARG prevalence across age and continent categories (hierarchical model)

fit_enrichment_age_continent_hier <- list()

# Fit hierarchical Bayesian logistic regression with random intercepts and slopes for age_category within continent
fit <- brm(
  high_ARG ~ age_category + (age_category | continent),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 2000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_age_continent_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_age_continent_hier, "fit_enrichment_age_continent_hier.rds")




## Enrichment analysis: High ARG prevalence across age, gender, and continent categories (hierarchical model)

fit_enrichment_income_age_gender_continent_hier <- list()

# Fit hierarchical Bayesian logistic regression with random intercepts for age_category within continent, income_group, and gender
fit <- brm(
  high_ARG ~ age_category + (1 | continent/income_group/gender),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 3000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_income_age_gender_continent_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_income_age_gender_continent_hier, "fit_enrichment_income_age_gender_continent_hier.rds")





## Enrichment analysis:  High ARG prevalence across age, gender, and continent categories with slopes (hierarchical model)

fit_enrichment_all_hier <- list()

# Fit hierarchical Bayesian logistic regression with random intercepts and slopes for age_category within income_group and gender
fit <- brm(
  high_ARG ~ age_category + (age_category | continent/income_group/gender),
  data = enrichment_df,
  family = bernoulli(),
  prior = set_prior("normal(0, 1)", class = "b"),
  chains = 2,
  iter = 2000,
  control = list(adapt_delta = 0.99),
  cores = parallel::detectCores()
)

fit_enrichment_all_hier[["all"]] <- fit

# Save hierarchical model
saveRDS(fit_enrichment_all_hier, "fit_enrichment_all_hier.rds")

