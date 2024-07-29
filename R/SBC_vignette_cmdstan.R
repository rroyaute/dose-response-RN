library(SBC)

use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", TRUE) # Set to false to use rstan instead

if(use_cmdstanr) {
  library(cmdstanr)
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
}

options(mc.cores = parallel::detectCores())

# Enabling parallel processing via future
library(future)
plan(multisession)

# The fits are very fast,
# so we force a minimum chunk size to reduce overhead of
# paralellization and decrease computation time.
options(SBC.min_chunk_size = 5)

# Setup caching of results
if(use_cmdstanr) {
  cache_dir <- "./_basic_usage_SBC_cache"
} else {
  cache_dir <- "./_basic_usage_rstan_SBC_cache"
}
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

cat(readLines("stan/poisson.stan"), sep = "\n")

if(use_cmdstanr) {
  cmdstan_model <- cmdstanr::cmdstan_model("stan/poisson.stan")
} else {
  rstan_model <- rstan::stan_model("stan/poisson.stan")
}

# A generator function should return a named list containing elements "variables" and "generated"

poisson_generator_single <- function(N){  # N is the number of data points we are generating
  lambda <- rgamma(n = 1, shape = 15, rate = 5)
  y <- rpois(n = N, lambda = lambda)
  list(
    variables = list(
      lambda = lambda
    ),
    generated = list(
      N = N,
      y = y
    )
  )
}

set.seed(54882235)
n_sims <- 100  # Number of SBC iterations to run

poisson_generator <- SBC_generator_function(poisson_generator_single, N = 40)
poisson_dataset <- generate_datasets(
  poisson_generator, 
  n_sims)

if(use_cmdstanr) {
  poisson_backend <- SBC_backend_cmdstan_sample(
    cmdstan_model, iter_warmup = 1000, iter_sampling = 1000, chains = 2)
} else {
  poisson_backend <- SBC_backend_rstan_sample(
    rstan_model, iter = 2000, warmup = 1000, chains = 2)  
}

results <- compute_SBC(poisson_dataset, poisson_backend, 
                       cache_mode = "results", 
                       cache_location = file.path(cache_dir, "results"))

results$stats
plot_rank_hist(results)
plot_ecdf(results)
plot_ecdf_diff(results)
