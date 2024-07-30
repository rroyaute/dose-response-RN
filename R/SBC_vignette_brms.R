library(SBC)
library(brms)
library(ggplot2)

# 0. Setup ----

use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", TRUE) # Set to false to use rstan instead
# use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", FALSE) # Set to false to use rstan instead

if(use_cmdstanr) {
  options(brms.backend = "cmdstanr")
} else {
  options(brms.backend = "rstan") 
  rstan::rstan_options(auto_write = TRUE)
}

# Using parallel processing
library(future)
plan(multisession)

# The fits are very fast,
# so we force a minimum chunk size to reduce overhead of
# paralellization and decrease computation time.
options(SBC.min_chunk_size = 5)

# Setup caching of results
if(use_cmdstanr) {
  cache_dir <- "./_brms_SBC_cache"
} else { 
  cache_dir <- "./_brms_rstan_SBC_cache"
}
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

# First option ----

# We need a "template dataset" to let brms build the model.
# The predictor (x) values will be used for data generation,
# the response (y) values will be ignored, but need to be present and 
# of the correct data type
set.seed(213452)
template_data = data.frame(y = rep(0, 15), x = rnorm(15))
priors <- prior(normal(0,1), class = "b") +
  prior(normal(0,1), class = "Intercept") +
  prior(normal(0,1), class = "sigma")
generator <- SBC_generator_brms(y ~ x, 
                                data = template_data, 
                                # backend = "rstan",
                                prior = priors, 
                                thin = 50, 
                                warmup = 10000, refresh = 2000,
                                out_stan_file = file.path(cache_dir, "brms_linreg1.stan")
)

set.seed(22133548)
datasets <- generate_datasets(generator, 100)

backend <- SBC_backend_brms_from_generator(generator, chains = 1, thin = 1,
                                           warmup = 500, iter = 1500,               
                                           inits = 0.1)

# More verbose alternative that results in exactly the same backend:
# backend <- SBC_backend_brms(y ~ x, template_data = template_data, prior = priors, warmup = 500, iter = 1000, chains = 1, thin = 1
#                            init = 0.1)

results <- compute_SBC(datasets, backend,
                       cache_mode = "results", 
                       cache_location = file.path(cache_dir, "first"))

plot_rank_hist(results)
plot_ecdf_diff(results)

# Second option ----
one_sim_generator <- function(N, K) {
  # N - number of datapoints, K number of groups for the varying intercept
  stopifnot(3 * K <= N)
  x <- rnorm(N) + 5
  
  group <- sample(1:K, size = N, replace = TRUE)
  # Ensure all groups are actually present at least twice
  group[1:(3*K)] <- rep(1:K, each = 3)
  
  b_Intercept <- rnorm(1, 5, 1)   
  b_x <- rnorm(1, 0, 1)
  
  sd_group__Intercept <- abs(rnorm(1, 0, 0.75))
  r_group <- matrix(rnorm(K, 0, sd_group__Intercept), 
                    nrow = K, ncol = 1,
                    dimnames = list(1:K, "Intercept"))
  
  sigma <- abs(rnorm(1, 0, 3))
  
  predictor <- b_Intercept + x * b_x + r_group[group]
  y <- rnorm(N, predictor, sigma)
  
  list(
    variables = list(
      b_Intercept = b_Intercept,
      b_x = b_x,
      sd_group__Intercept = sd_group__Intercept,
      r_group = r_group,
      sigma = sigma
    ),
    generated = data.frame(y = y, x = x, group = group)
  )
}

n_sims_generator <- SBC_generator_function(one_sim_generator, N = 18, K = 5)

log_lik_dq_func <- derived_quantities(
  log_lik = sum(dnorm(y, b_Intercept + x * b_x + r_group[group], sigma, log = TRUE))
  # Testing CRPS, probably not worth it
  #, CRPS = mean(scoringRules::crps_norm(y, b_Intercept + x * b_x + r_group[group], sigma))
)

set.seed(12239755)
datasets_func <- generate_datasets(n_sims_generator, 100)

priors_func <- prior(normal(0,1), class = "b") +
  prior(normal(5,1), class = "Intercept") +
  prior(normal(0,5), class = "sigma") +
  prior(normal(0,0.75), class = "sd")


backend_func <- SBC_backend_brms(y ~ x + (1 | group), 
                                 # backend = "cmdstanr", 
                                 prior = priors_func, 
                                 chains = 1,
                                 file_refit = "always",
                                 template_data = datasets_func$generated[[1]],
                                 out_stan_file = file.path(cache_dir, "brms_linreg2.stan"))

results_func <- compute_SBC(datasets_func, backend_func, 
                            dquants = log_lik_dq_func, 
                            cache_mode = "results", 
                            cache_location = file.path(cache_dir, "func"))

plot_rank_hist(results_func)
plot_ecdf_diff(results_func)
