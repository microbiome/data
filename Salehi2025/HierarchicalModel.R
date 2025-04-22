
# Fit the hierarchical model for ARG load
fit_ARG_hierarchical <- brm(
  formula = bf(ARG_load ~ age_category + (1 + age_category | income_group_HIC / gender)),
  data = subset_data,
  family = lognormal(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("normal(0, 1)", class = "Intercept")
  ),
  chains = 2,
  iter = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  cores = parallel::detectCores()
)

# Save the fitted ARG load model
saveRDS(fit_ARG_hierarchical, "fit_ARG_hierarchical.rds")



# Fit the hierarchical model for Shannon diversity
fit_shannon_hierarchical <- brm(
  formula = bf(shannon_diversity ~ age_category + (1 + age_category | income_group_HIC / gender)),
  data = subset_data,
  family = gaussian(),
  prior = c(
    set_prior("normal(0, 1)", class = "b"),
    set_prior("normal(0, 1)", class = "Intercept")
  ),
  chains = 2,
  iter = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  cores = parallel::detectCores()
)

# Save the fitted Shannon diversity model
saveRDS(fit_shannon_hierarchical, "fit_shannon_hierarchical.rds")

