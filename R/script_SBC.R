library(tidyverse); library(viridis); library(brms); library(tidybayes)
library(patchwork); library(SBC); library(future)

bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

# use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", TRUE) # Set to false to use rstan instead
use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", FALSE) # Set to false to use rstan instead

if(use_cmdstanr) {
  options(brms.backend = "cmdstanr")
} else {
  options(brms.backend = "rstan") 
  rstan::rstan_options(auto_write = TRUE)
}

# Using parallel processing
plan(multisession)

# The fits are very fast,
# so we force a minimum chunk size to reduce overhead of
# paralellization and decrease computation time.
options(SBC.min_chunk_size = 5)

# Setup caching of results
if(use_cmdstanr) {
  cache_dir <- "./_SBC_cache_DR_RN"
} else { 
  cache_dir <- "./_SBC_cache_DR_RN_rstan"
}
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

# Simulate data
D = seq(0, 100, by = 10)
I = 10
N = length(D) * I

sim_generator = function(N, I) {
  # N - number of datapoints, I number of individuals 
  D = seq(0, 100, by = 10)
  Dose = rep(D, I)
  
  ID = rep(1:I, each = length(D))
  
  b_alpha_Intercept = rnorm(1, 100, 20)
  b_beta_Intercept = rexp(1, 10) # positive slope mean = .1, sd = .1
  b_NEC_Intercept = runif(1, 0, 100)
  sigma = rexp(1, 1)
  sd_ID__alpha_Intercept = rexp(1, .05)
  sd_ID__beta_Intercept = rexp(1, 5)
  sd_ID__NEC_Intercept = rexp(1, .1)
  
  # Sample from lkj(4) for a 3x3 correlation matrix
  cor_mat = rethinking::rlkjcorr(1, eta = 4, K = 3)
  
  # Plug in correlations
  cor_ID__alpha_Intercept__beta_Intercept = cor_mat[2,1]
  cor_ID__alpha_Intercept__NEC_Intercept = cor_mat[3,1]
  cor_ID__beta_Intercept__NEC_Intercept = cor_mat[2,3]
  
  r_ID__alpha = matrix(rnorm(I, 0, sd_ID__alpha_Intercept), 
                       nrow = I, ncol = 1,
                       dimnames = list(1:I, "Intercept"))
  r_ID__beta_Intercept = matrix(rnorm(I, 0, sd_ID__beta_Intercept), 
                       nrow = I, ncol = 1,
                       dimnames = list(1:I, "Intercept"))
  r_ID__NEC_Intercept = matrix(rnorm(I, 0, sd_ID__NEC_Intercept), 
                       nrow = I, ncol = 1,
                       dimnames = list(1:I, "Intercept"))
  
  # Need to add cor variables somehow
  
  log_yhat = log(b_alpha_Intercept + r_ID__alpha[ID]) - 
    (b_beta_Intercept + r_ID__beta_Intercept[ID]) * 
    (Dose - (b_NEC_Intercept + r_ID__NEC_Intercept[ID])) * 
    (Dose > (b_NEC_Intercept + r_ID__NEC_Intercept[ID]))
  y = rlnorm(N, log_yhat, sigma)
  
  list(
    variables = list(
      b_alpha_Intercept = b_alpha_Intercept,
      b_beta_Intercept = b_beta_Intercept,
      b_NEC_Intercept = b_NEC_Intercept,
      sd_ID__alpha_Intercept = sd_ID__alpha_Intercept,
      sd_ID__beta_Intercept = sd_ID__beta_Intercept,
      sd_ID__NEC_Intercept = sd_ID__NEC_Intercept,
      cor_ID__alpha_Intercept__beta_Intercept = cor_ID__alpha_Intercept__beta_Intercept,
      cor_ID__alpha_Intercept__NEC_Intercept = cor_ID__alpha_Intercept__NEC_Intercept,
      cor_ID__beta_Intercept__NEC_Intercept = cor_ID__beta_Intercept__NEC_Intercept,
      sigma = sigma#,
      # r_ID__alpha = r_ID__alpha,
      # r_ID__beta_Intercept = r_ID__beta_Intercept,
      # r_ID__NEC_Intercept = r_ID__NEC_Intercept
    ),
    generated = data.frame(y = y, Dose = Dose, ID = ID)
  )
}

n_sims_generator = SBC_generator_function(sim_generator, N = N, I = I)

log_lik_dq_func <- derived_quantities(
  log_lik = sum(dnorm(y, 
                      log(b_alpha_Intercept + r_ID__alpha[ID]) - 
                        (b_beta_Intercept + r_ID__beta_Intercept[ID]) * 
                        (Dose - (b_NEC_Intercept + r_ID__NEC_Intercept[ID])) * 
                        (Dose > (b_NEC_Intercept + r_ID__NEC_Intercept[ID])), 
                      sigma, log = TRUE))
  # Testing CRPS, probably not worth it
  #, CRPS = mean(scoringRules::crps_norm(y, b_Intercept + x * b_x + r_group[group], sigma))
)

set.seed(12239755)
datasets_func = generate_datasets(n_sims_generator, 100)

priors.vi.2 = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Random effects priors (Weakly informative ~ 20 % of mean parameter)
  prior(exponential(.05), nlpar = alpha, class = sd, group = ID) +
  prior(exponential(5), nlpar = beta, class = sd, group = ID) +
  prior(exponential(.1), nlpar = NEC, class = sd, group = ID) +
  # Residual prior
  prior(exponential(1), class = sigma) +
  prior(lkj(4), class = cor)

bf.vi = bf(y ~ log(alpha) - beta * (Dose - NEC) * (Dose > NEC), 
           alpha + beta + NEC ~ 1 + (1|c|ID), 
           nl = T, 
           family = lognormal)

backend_func <- SBC_backend_brms(bf.vi,  
                                 prior = priors.vi.2, 
                                 control = list(adapt_delta = .95,
                                                max_treedepth = 12),
                                 chains = 4,
                                 template_data = datasets_func$generated[[1]],
                                 out_stan_file = file.path(cache_dir, "brms_NEC_vi.stan"))

results_func <- compute_SBC(datasets_func, 
                            backend_func, 
                            # dquants = log_lik_dq_func,
                            cache_mode = "results", 
                            cache_location = file.path(cache_dir, "func"))
plot_rank_hist(results_func)
plot_ecdf_diff(results_func)
plot_coverage(results_func)
plot_coverage_diff(results_func)
plot_sim_estimated(results_func, alpha = 0.5)
