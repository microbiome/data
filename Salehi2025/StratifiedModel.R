# Model for ARG load stratified by income and sex

fit_ARG_load_age <- list()

incomes <- 0:1
sex_levels <- c("Women", "Men")

# Set reference level for age_category
adult_metadata$age_category <- relevel(factor(adult_metadata$age_category), ref = "Middle-Aged Adult")

for (ic in incomes) {
  for (sex_group in sex_levels) {
    
    # Filter data by income group and gender
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic, gender == sex_group)
    
    # Use only age_category as predictor
    my_formula <- as.formula("ARG_load ~ age_category")
    
    print(paste("Fitting ARG_load model for:", ic, "-", sex_group))
    
    my_fit <- brm(
      formula = my_formula,
      data = my_data,
      family = lognormal(),
      prior = c(
        set_prior("normal(0, 1)", class = "b"),
        set_prior("normal(0, 1)", class = "Intercept")
      ),
      chains = 2,
      iter = 5000,
      control = list(adapt_delta = 0.99, max_treedepth = 11),
      cores = parallel::detectCores()
    )
    
    fit_ARG_load_age[[paste0(ic, "_", sex_group)]] <- my_fit
  }
}

# Save model
saveRDS(fit_ARG_load_age, "fit_ARG_load_age.rds")



# Model for shannon_diversity stratified by income and sex

fit_shannon_diversity_age <- list()

for (ic in incomes) {
  for (sex_group in sex_levels) {
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic, gender == sex_group)
    
    # Use only age_category as predictor
    my_formula <- as.formula("shannon_diversity ~ age_category")
    
    print(paste("Fitting shannon_diversity model for:", ic, "-", sex_group))
    
    my_fit <- brm(
      formula = my_formula,
      data = my_data,
      family = gaussian(),
      prior = c(
        set_prior("normal(0, 1)", class = "b"),
        set_prior("normal(0, 1)", class = "Intercept")
      ),
      chains = 2,
      iter = 5000,
      control = list(adapt_delta = 0.99, max_treedepth = 11),
      cores = parallel::detectCores()
    )
    
    fit_shannon_diversity_age[[paste0(ic, "_", sex_group)]] <- my_fit
  }
}

# Save model
saveRDS(fit_shannon_diversity_age, "fit_shannon_diversity_age.rds")
