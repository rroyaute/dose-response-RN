
# Dose-response function on log(Response) scale ----
DR_logRN_fun = function(Dose, Rmin, Rmax, beta, NEC){
  logyhat = Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)) 
  yhat = exp(logyhat)
  return(yhat)
}

# Function to generate 1 dataset according to prior distribution of parameters ----
data_generator = function(D = seq(0, 100, by = 10), 
                          I = 10, 
                          CV = .05, 
                          CVi = .1,
                          mu_pars = list(Rmin = log(1), 
                                         Rmax = log(100), 
                                         beta = -4.5,
                                         NEC = 10),
                          sigma_par = 1/10,
                          ID_pars = c("Rmax", "NEC", "beta", "all", "cov")) {
  # N - number of datapoints, I number of "individuals" 
  Dose = rep(D, I)
  N = length(D) * I
  
  ID = rep(1:I, each = length(D))
  
  # Assume small variation around mean value
  # CV set to 75 % variation around the mean by default
  b_Rmin_Intercept = mu_pars$Rmin # Fix the minimum value
  # Alternative if we want to sample Rmin from a prior distribution
  # b_Rmin_Intercept = rtruncnorm(1, a = 1, b = Inf, mu_pars$Rmin, mu_pars$Rmin * CV)
  b_Rmax_Intercept = rtruncnorm(1, a = 0, b = Inf,  mu_pars$Rmax, mu_pars$Rmax * CV)
  b_beta_Intercept = rnorm(1, mu_pars$beta, abs(mu_pars$beta * CV))
  b_NEC_Intercept = mu_pars$NEC # Fix the NEC
  # Alternative if we want to sample NEC from prior distribution
  # Assume NEC is anywhere between 0 and maximum dose
  # Set CV = 75 % variation as reasonable 
  # b_NEC_Intercept = rtruncnorm(1, a = min(D), b = max(D), 
  #                              mean(D), mean(D) * .75)
  
  sigma = rexp(1, 1/sigma_par) # Sample residual variance with relatively small value
  # Fix individual variance to mean parameter x CVi
  # For Rmax, and beta the CV on the natural scale is roughly equivalent to sd
  # For beta, we take the absolute value multiplied by the set CV
  sd_ID__Rmax_Intercept = ifelse(ID_pars == "Rmax" | 
                                   ID_pars == "all" | 
                                   ID_pars == "cov",
                                 CVi, 0)
  sd_ID__NEC_Intercept = ifelse(ID_pars == "NEC" | 
                                  ID_pars == "all" | 
                                  ID_pars == "cov",
                                mu_pars$NEC * CVi, 0)
  sd_ID__beta_Intercept = ifelse(ID_pars == "beta" | 
                                   ID_pars == "all" | 
                                   ID_pars == "cov",
                                 CVi, 0)
  
  # Define correlation matrix among varying parameters (TO COMPLETE)
  # Sample from lkj(1) for a 3x3 correlation matrix
  # cor_mat = rethinking::rlkjcorr(1, eta = 1, K = 3)
  # Plug in correlations 
  # cor_ID__Rmax_Intercept__NEC_Intercept = cor_mat[1,2]
  # cor_ID__Rmax_Intercept__beta_Intercept = cor_mat[1,3]
  # cor_ID__beta_Intercept__NEC_Intercept = cor_mat[2,3]
  
  # Store individual deviations from population value
  r_ID__Rmax = matrix(rnorm(I, 0, sd_ID__Rmax_Intercept), 
                      nrow = I, ncol = 1,
                      dimnames = list(1:I, "Intercept"))
  r_ID__beta_Intercept = matrix(rnorm(I, 0, sd_ID__beta_Intercept), 
                                nrow = I, ncol = 1,
                                dimnames = list(1:I, "Intercept"))
  r_ID__NEC_Intercept = matrix(rnorm(I, 0, sd_ID__NEC_Intercept), 
                               nrow = I, ncol = 1,
                               dimnames = list(1:I, "Intercept"))
  
  # Need to add cor variables somehow
  log_yhat = (b_Rmin_Intercept) +
    ((b_Rmax_Intercept + r_ID__Rmax[ID]) - (b_Rmin_Intercept)) *
    exp(-exp(b_beta_Intercept + r_ID__beta_Intercept[ID]) * 
          (Dose - (b_NEC_Intercept + r_ID__NEC_Intercept[ID])) *
          (Dose > (b_NEC_Intercept + r_ID__NEC_Intercept[ID])))
  y = rlnorm(N, log_yhat, sigma)
  
  list(
    variables = list(
      b_Rmin_Intercept = b_Rmin_Intercept,
      b_Rmax_Intercept = b_Rmax_Intercept,
      b_beta_Intercept = b_beta_Intercept,
      b_NEC_Intercept = b_NEC_Intercept,
      sd_ID__Rmax_Intercept = sd_ID__Rmax_Intercept,
      sd_ID__beta_Intercept = sd_ID__beta_Intercept,
      sd_ID__NEC_Intercept = sd_ID__NEC_Intercept,
      # cor_ID__Rmax_Intercept__NEC_Intercept = cor_ID__Rmax_Intercept__NEC_Intercept,
      # cor_ID__Rmax_Intercept__beta_Intercept = cor_ID__Rmax_Intercept__beta_Intercept,
      # cor_ID__beta_Intercept__NEC_Intercept = cor_ID__beta_Intercept__NEC_Intercept,
      sigma = sigma,
      r_ID__Rmax = r_ID__Rmax,
      r_ID__NEC_Intercept = r_ID__NEC_Intercept,
      r_ID__beta_Intercept = r_ID__beta_Intercept
      
    ),
    generated = data.frame(y = y, Dose = Dose, ID = ID)
  )
}

# Generate  n_sim datasets ----

# n_sims_generator = SBC_generator_function(data_generator,
#                                           D = seq(0, 100, by = 10), 
#                                           I = 10, 
#                                           CV = .05, 
#                                           CVi = .1,
#                                           mu_pars = list(Rmin = log(1), 
#                                                          Rmax = log(100), 
#                                                          beta = -4.5,
#                                                          NEC = 10),
#                                           ID_pars = c("Rmax", "NEC", "beta", "all", "cov"))

# Function to define priors to pass on to brms ----
priors_func = function(mu_pars = c(Rmin = log(1), 
                                   Rmax = log(100), 
                                   beta = -4.5,
                                   NEC = 10),
                       CV_prior = .5, 
                       sigma = 1/10
) {
  priors_pop = mu_pars
  priors_sd[1:3] = .5 # set sd to roughly correspond to 50 % variation on the natural scale for Rmin, Rmax & beta
  priors_sd[4] = abs(mu_pars[4]) * CV_prior
  priors_sigma = sigma
  
  prior_params = stanvar(priors_pop, 'priors_pop') + 
    stanvar(priors_sd, 'priors_sd') + 
    stanvar(priors_sigma, 'priors_sigma')
  
  priors.vi = 
    # Intercept priors (50 % variation around mean value)
    prior(normal(priors_pop[1], priors_sd[1]), nlpar = Rmin, class = b, lb = 0) +
    prior(normal(priors_pop[2], priors_sd[2]), nlpar = Rmax, class = b, lb = 0) +
    prior(normal(priors_pop[3], priors_sd[3]), nlpar = beta, class = b) +
    prior(normal(priors_pop[4], priors_sd[4]), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
    # Random effects priors
    prior(exponential(1), nlpar = NEC, class = sd, group = ID) +
    prior(exponential(1), nlpar = Rmax, class = sd, group = ID) +
    prior(exponential(1), nlpar = beta, class = sd, group = ID) +
    # Residual prior
    prior(exponential(10), class = sigma)
  
  priors.pop = priors.vi %>% filter(class != "sd")
  prior.list = list(priors.vi = priors.vi, priors.pop = priors.pop)
  return(prior.list)
}



# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint


# Store bf model list ----
bf_list = list(
  bf.vi.Rmax = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
                  Rmin + NEC + beta ~ 1,
                  Rmax ~ 1 + (1|ID),
                  family = lognormal,
                  nl = T),
  
  bf.vi.NEC = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
                  Rmin + Rmax + beta ~ 1,
                  NEC ~ 1 + (1|ID),
                  family = lognormal,
                  nl = T),
  
  bf.vi.beta = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
                 Rmin + Rmax + NEC ~ 1,
                 beta ~ 1 + (1|ID),
                 family = lognormal,
                 nl = T),
  
  bf.vi.all = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
                 Rmin ~ 1,
                 Rmax + NEC + beta ~ 1 + (1|ID),
                 family = lognormal,
                 nl = T),
  
  bf.vi.cov = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
                 Rmin ~1,
                 Rmax + NEC + beta ~ 1 + (1|c|ID),
                 family = lognormal,
                 nl = T))
  
  


# Function to fit the models to simulated data ----
n_sim = 100
n_id = 20
CVi = c(0, .1, .2, .3, .4, .5)
case = c("Rmax", "NEC", "beta", "all") # "cov" not yet implemented
mod = c("pop", "vi")

n_sims_generator = SBC_generator_function(f = data_generator,
                       D = seq(0, 100, by = 10), 
                       I = 10, 
                       CV = .05, 
                       CVi = .1,
                       mu_pars = list(Rmin = log(1), 
                                      Rmax = log(100), 
                                      beta = -4.5,
                                      NEC = 10),
                       sigma_par = 1/10,
                       ID_pars = c("Rmax", "NEC", "beta", "all", "cov"))

datasets_list = list(
  df_list_Rmax = list(
    CVi.0 = generate_datasets(
      n_sims_generator(ID_pars = case[1], I = n_id, CV = .05,
                       CVi = CVi[1]), n_sims = n_sim)
      
    
    
  )
)
  
  
  
  generate_datasets(
  , 
                                  n_sims = n_sims)


datasets_func = generate_datasets(data_generator(), n_sims)

data_list = list(
  data_list_Rmax = generate_datasets(n_sims_generator, n_sims = n_sims)
    
  datasets_func$generated[[1]]
)



backend_func_Rmax = function(mod = c("pop", "vi")){
  
}


backend_func_Rmax <- SBC_backend_brms(bf.vi,  
                                 prior = priors.vi.2, 
                                 control = list(adapt_delta = .95,
                                                max_treedepth = 12),
                                 chains = 4,
                                 template_data = datasets_func$generated[[1]],
                                 out_stan_file = file.path(cache_dir, "brms_NEC_vi.stan"))
