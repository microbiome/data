# First, we define the income groups (LMIC, HIC) used for stratified statistical modeling.
incomes <- 0:1
set.seed(123123)



## Model 1 for ARG load **************************** ####

fit1_ARG_load <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>%
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  my_formula <- my_formula %>% as.formula()
  
  print("ARG_load Model 1")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = lognormal(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit1_ARG_load[[ic+1]] <- my_fit
}

# Save ARG_load model
saveRDS(fit1_ARG_load, "fit1_ARG_load.rds")


## Model 1 for shannon_diversity ****************************

fit1_shannon_diversity <- list()

for(ic in incomes) {    # Loop over income levels for shannon_diversity
  
  # Filter data
  my_data <- adult_df %>%
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  my_formula <- my_formula %>% as.formula()
  
  print("shannon_diversity Model 1")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit1_shannon_diversity[[ic+1]] <- my_fit
}

# Save shannon_diversity model
saveRDS(fit1_shannon_diversity, "fit1_shannon_diversity.rds")






## Model 2 for ARG_load **************************** ####

fit2_ARG_load <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add random intercept for bioproject
  my_formula <- paste0(my_formula, " + (1 | bioproject)") %>% as.formula()
  
  print("ARG_load Model 2")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = lognormal(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept"),
                  set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit2_ARG_load[[ic+1]] <- my_fit
}

# Save ARG_load Model 2
saveRDS(fit2_ARG_load, "fit2_ARG_load.rds")



## Model 2 for shannon_diversity ****************************

fit2_shannon_diversity <- list()

for(ic in incomes) {    # Loop over income levels for shannon_diversity
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add random intercept for bioproject
  my_formula <- paste0(my_formula, " + (1 | bioproject)") %>% as.formula()
  
  print("shannon_diversity Model 2")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept"),
                  set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit2_shannon_diversity[[ic+1]] <- my_fit
}

# Save shannon_diversity Model 2
saveRDS(fit2_shannon_diversity, "fit2_shannon_diversity.rds")











## Model 3 for ARG_load **************************** ####

fit3_ARG_load <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add log(readcount) as covariate
  my_formula <- paste0(my_formula, " + log(readcount)") %>% as.formula()
  
  print("ARG_load Model 3")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = lognormal(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit3_ARG_load[[ic+1]] <- my_fit
}

# Save ARG_load Model 3
saveRDS(fit3_ARG_load, "fit3_ARG_load.rds")



## Model 3 for shannon_diversity ****************************

fit3_shannon_diversity <- list()

for(ic in incomes) {    # Loop over income levels for shannon_diversity
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add log(readcount) as covariate
  my_formula <- paste0(my_formula, " + log(readcount)") %>% as.formula()
  
  print("shannon_diversity Model 3")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit3_shannon_diversity[[ic+1]] <- my_fit
}

# Save shannon_diversity Model 3
saveRDS(fit3_shannon_diversity, "fit3_shannon_diversity.rds")







## Model 4 for ARG_load **************************** ####

fit4_ARG_load <- list()

for(ic in incomes) {    # Loop over income levels for ARG_load
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("ARG_load ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add log(readcount) as covariate and random intercept for bioproject
  my_formula <- paste0(my_formula, " + log(readcount) + (1 | bioproject)") %>% as.formula()
  
  print("ARG_load Model 4")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = lognormal(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept"),
                  set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit4_ARG_load[[ic+1]] <- my_fit
}

# Save ARG_load Model 4
saveRDS(fit4_ARG_load, "fit4_ARG_load.rds")


## Model 4 for shannon_diversity ****************************

fit4_shannon_diversity <- list()

for(ic in incomes) {    # Loop over income levels for shannon_diversity
  
  # Filter data
  my_data <- adult_df %>% 
    filter(income_group_HIC == ic)
  
  if(ic == 0) {
    # Remove Africa, North America and Oceania for LMIC
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  if(ic == 1) {
    # Remove Africa, South America, set Europe to base level
    my_formula <- paste0("shannon_diversity ~ ",
                         paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                       x = dummy_var_names)],
                                collapse = "+"))
  }
  
  # Add log(readcount) as covariate and random intercept for bioproject
  my_formula <- paste0(my_formula, " + log(readcount) + (1 | bioproject)") %>% as.formula()
  
  print("shannon_diversity Model 4")
  my_fit <- brm(formula = my_formula,
                data = my_data,
                family = gaussian(),
                prior = c(
                  set_prior("normal(0, 1)", class = "b"),
                  set_prior("normal(0, 1)", class = "Intercept"),
                  set_prior("normal(0, 1)", class = "sd", group = "bioproject")
                ),
                chains = 2,
                iter = 5000,
                control = list(adapt_delta = 0.99, max_treedepth = 11),
                cores = parallel::detectCores())
  
  fit4_shannon_diversity[[ic+1]] <- my_fit
}

# Save shannon_diversity Model 4
saveRDS(fit4_shannon_diversity, "fit4_shannon_diversity.rds")
