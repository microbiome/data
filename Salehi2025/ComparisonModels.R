#### Models for Comparison of Log-normal and Log-transformed Normal Models ####

# First, we define the income groups (LMIC, HIC) used for stratified statistical modeling.
incomes <- c("LMIC", "HIC")
set.seed(123123)


## Model for ARG load (log-normal)

fit_ARG_lognormal <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group == ic)
  
  if(ic == "LMIC") {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == "HIC") {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  my_formula <- my_formula %>% as.formula()
  
  print("ARG_load Model (lognormal)")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = lognormal(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 2000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit_ARG_lognormal[[ic]] <- my_fit
}

# Save ARG_load lognormal model
saveRDS(fit_ARG_lognormal, "fit_ARG_lognormal.rds")




## Model for ARG load (Gaussian on log-transformed data)

fit_ARG_gaussian <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data and log-transform the data
  my_data <- adult_df %>% 
    filter(income_group == ic) %>%
    mutate(log_ARG_load = log(ARG_load))  # Log-transform the ARG load
  
  if(ic == "LMIC") {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("log_ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == "HIC") {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("log_ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  my_formula <- my_formula %>% as.formula()
  
  print("ARG_load Model (Gaussian on log-transformed data)")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 2000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit_ARG_gaussian[[ic]] <- my_fit
}

# Save ARG_load gaussian model
saveRDS(fit_ARG_gaussian, "fit_ARG_gaussian.rds")
