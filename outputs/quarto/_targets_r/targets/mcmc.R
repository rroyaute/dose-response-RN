tar_stan_mcmc_rep_summary(
  name = mcmc,
  stan_files = "model.stan",
  data = simulate_data(), # Runs once per rep.
  batches = 100, # Number of branch targets.
  reps = 10, # Number of model reps per branch target.
  chains = 4, # Number of MCMC chains.
  parallel_chains = 4, # How many MCMC chains to run in parallel.
  iter_warmup = 4e4, # Number of MCMC warmup iterations to run.
  iter_sampling = 4e4, # Number of MCMC post-warmup iterations to run.
  summaries = list(
    # Compute posterior intervals at levels 50% and 95%.
    # The 50% intervals should cover prior predictive parameter draws
    # 50% of the time. The 95% intervals are similar.
    # We also calculate posterior medians so we can compare them
    # directly to the prior predictive draws.
    ~posterior::quantile2(.x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
    # We use Gelman-Rubin potential scale reduction factors to
    # assess convergence:
    rhat = ~posterior::rhat(.x)
  ),
  deployment = "worker"
)
