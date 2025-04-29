#### Models for Comparison of Log-normal and Log-transformed Normal Models ####

# Define income groups for stratified modeling
incomes <- c("LMIC", "HIC")
set.seed(123123)

## Model for ARG load (log-normal)
fit_ARG_lognormal <- list()

for(ic in incomes) {  # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>%
    filter(income_group == ic)
  
  # Set reference levels for factors
  continent_base <- ifelse(ic == "LMIC", "Asia", "Europe")
  my_data <- my_data %>%
    mutate(
      continent = relevel(factor(continent), ref = continent_base)
    )
  
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