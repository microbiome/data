# Model for shannon_diversity stratified by income and sex

fit_shannon_diversity_age <- list()

incomes <- c("LMIC", "HIC")
sex_levels <- c("Women", "Men")

for (ic in incomes) {
  for (sex_group in sex_levels) {
    
    my_data <- adult_df %>%
      filter(income_group == ic, gender == sex_group)
    
    my_formula <- shannon_diversity ~ age_category
    
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
      iter = 2000,
      control = list(adapt_delta = 0.99, max_treedepth = 11),
      cores = parallel::detectCores()
    )
    
    fit_shannon_diversity_age[[paste(ic, sex_group, sep = "_")]] <- my_fit
  }
}

saveRDS(fit_shannon_diversity_age, "fit_shannon_diversity_age.rds")

