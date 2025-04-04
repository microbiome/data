
# First, we define the response variables (ARG_load, shannon_diversity) and 
# income groups (LMIC, HIC) used for stratified statistical modeling.

responses <- c("ARG_load", "shannon_diversity")
incomes <- 0:1
set.seed(123123)



## Model 1 **************************** ####

fit_list1 <- list()

for(r in responses) {     # Loop over responses
  for(ic in incomes) {    # Loop over income levels
    
    # Filter data
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    if(ic == 0) {
      
      # Remove Africa, North America and Oceania since no observation for LMIC
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "Oceania|North_America|Africa|Asia|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
      
    }
    
    if(ic == 1) {
      # Remove Africa, South America, set Europe to base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl(pattern = "South_America|Africa|Europe|HIC",
                                                         x = dummy_var_names)],
                                  collapse = "+"))
    }
    
    my_formula <- my_formula %>% as.formula()
    
    if(r == "shannon_diversity") {
      print(r)
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
      
    } else if(r == "ARG_load") {
      print(r)
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
      
    }
    
    
    fit_list1[[r]][[ic+1]] <- my_fit
    
  }
}

# Save model 1
saveRDS(fit_list1, "fit1.rds")




## Model 2 **************************** ####

fit_list2 <- list()

for (r in responses) {  # Loop over responses
  for (ic in incomes) {  # Loop over income levels
    
    # Filter data based on income group
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    # Define formula based on income group
    if (ic == 0) {
      # Remove Africa, North America, and Oceania (no LMIC observations)
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("Oceania|North_America|Africa|Asia|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    if (ic == 1) {
      # Remove Africa, South America; set Europe as base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("South_America|Africa|Europe|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    # Add random intercept for bioproject
    my_formula <- paste0(my_formula, " + (1 | bioproject)") %>% as.formula()
    
    # Fit model based on response variable
    if (r == "shannon_diversity") {
      print(r)
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
      
    } else if (r == "ARG_load") {
      print(r)
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
      
    }
    
    # Store model in list
    fit_list2[[r]][[ic+1]] <- my_fit
  }
}

# Save model 2
saveRDS(fit_list2, "fit2.rds")




## Model 3 **************************** ####

fit_list3 <- list()

for (r in responses) {  # Loop over responses
  for (ic in incomes) {  # Loop over income levels
    
    # Filter data based on income group
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    # Define formula based on income group
    if (ic == 0) {
      # Remove Africa, North America, and Oceania (no LMIC observations)
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("Oceania|North_America|Africa|Asia|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    if (ic == 1) {
      # Remove Africa, South America; set Europe as base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("South_America|Africa|Europe|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    # Add read count as covariate
    my_formula <- paste0(my_formula, " + log(readcount)") %>% as.formula()
    
    # Fit model based on response variable
    if (r == "shannon_diversity") {
      print(r)
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
      
    } else if (r == "ARG_load") {
      print(r)
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
      
    }
    
    # Store model in list
    fit_list3[[r]][[ic+1]] <- my_fit
  }
}

# Save model 3
saveRDS(fit_list3, "fit3.rds")


## Model 4 **************************** ####

fit_list4 <- list()

for (r in responses) {  # Loop over responses
  for (ic in incomes) {  # Loop over income levels
    
    # Filter data based on income group
    my_data <- adult_metadata %>%
      as.data.frame() %>%
      filter(income_group_HIC == ic)
    
    # Define formula based on income group
    if (ic == 0) {
      # Remove Africa, North America, and Oceania (no LMIC observations)
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("Oceania|North_America|Africa|Asia|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    if (ic == 1) {
      # Remove Africa, South America; set Europe as base level
      my_formula <- paste0(r, " ~ ",
                           paste0(dummy_var_names[!grepl("South_America|Africa|Europe|HIC", 
                                                         x = dummy_var_names)], 
                                  collapse = "+"))
    }
    
    # Add readcount as covariate and random intercept for bioproject
    my_formula <- paste0(my_formula, " + log(readcount) + (1 | bioproject)") %>% as.formula()
    
    # Fit model based on response variable
    if (r == "shannon_diversity") {
      print(r)
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
      
    } else if (r == "ARG_load") {
      print(r)
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
      
    }
    
    # Store model in list
    fit_list4[[r]][[ic+1]] <- my_fit
  }
}

# Save model 4
saveRDS(fit_list4, "fit4.rds")