#### Models for Comparison of Log-normal and Log-transformed Normal Models ####

# Define income groups for stratified modeling
incomes <- c("LMIC", "HIC")
set.seed(123123)

## Model for ARG load (log-normal)
fit_ARG_lognormal <- list()

for(ic in incomes) {  # Loop over income levels for ARG_load
  
  # Select data with correct pre-set base level
  my_data <- if (ic == "LMIC") adult_df_LMIC else adult_df_HIC
  
  # Define fixed effects model
  my_formula <- ARG_load ~ continent + age_category + gender + Usage_high
  
  print(paste("Fitting ARG_load Model (lognormal) for", ic))
  my_fit <- brm(
    formula = my_formula,
    data = my_data,
    family = lognormal(),
    prior = c(
      set_prior("normal(0, 1)", class = "b"),
      set_prior("normal(0, 1)", class = "Intercept")
    ),
    chains = 2,
    iter = 2000,
    control = list(adapt_delta = 0.99, max_treedepth = 11),
    cores = parallel::detectCores()
  )
  
  fit_ARG_lognormal[[ic]] <- my_fit
}

# Save ARG_load lognormal model
saveRDS(fit_ARG_lognormal, "fit_ARG_lognormal.rds")