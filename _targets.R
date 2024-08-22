library(targets)
library(brms)
library(SBC)

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
  b_NEC_Intercept = mu_pars$Rmin # Fix the NEC
  # Alternative if we want to sample NEC from prior distribution
  # Assume NEC is anywhere between 0 and maximum dose
  # Set CV = 75 % variation as reasonable 
  # b_NEC_Intercept = rtruncnorm(1, a = min(D), b = max(D), 
  #                              mean(D), mean(D) * .75)
  
  sigma = rexp(1, 1/sigma_par) # Sample residual variance with relatively small value
  # Fix individual variance to mean parameter x CVi
  sd_ID__Rmax_Intercept = ifelse(ID_pars == "Rmax" | 
      ID_pars == "all" | 
      ID_pars == "cov",
      mu_pars$Rmax * CVi, 0)
  sd_ID__NEC_Intercept = ifelse(ID_pars == "NEC" | 
           ID_pars == "all" | 
           ID_pars == "cov",
         mu_pars$NEC * CVi, 0)
  sd_ID__beta_Intercept = ifelse(ID_pars == "beta" | 
           ID_pars == "all" | 
           ID_pars == "cov",
           abs(mu_pars$beta) * CVi, 0)
  
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
      sigma = sigma#,
      # r_ID__Rmax = r_ID__Rmax,
      # r_ID__NEC_Intercept = r_ID__NEC_Intercept
      # r_ID__beta_Intercept = r_ID__beta_Intercept,
      
    ),
    generated = data.frame(y = y, Dose = Dose, ID = ID)
  )
}

# Function to store n_sim datasets ----
n_sims_generator = SBC_generator_function(data_generator,
                                          D = seq(0, 100, by = 10), 
                                          I = 10, 
                                          CV = .05, 
                                          CVi = .1,
                                          mu_pars = list(Rmin = log(1), 
                                                         Rmax = log(100), 
                                                         beta = -4.5,
                                                         NEC = 10),
                                          ID_pars = c("Rmax", "NEC", "beta", "all", "cov"))
datasets_func = generate_datasets(n_sims_generator, nsims = 100)

# Function to define priors to pass on to brms
priors_func = function(mu_pars = c(Rmin = log(1), 
                                   Rmax = log(100), 
                                   beta = -4.5,
                                   NEC = 10),
                       CV_prior = .5, 
                       sigma = 1/10
                       ) {
  priors_pop = mu_pars
  priors_sd = abs(mu_pars) * CV_prior
  priors_sd[1] = .5 # set sd to roughly correspond to 50 % variation on the natural scale for Rmin
  priors_sd[2] = .5 # set sd to roughly correspond to 50 % variation on the natural scale for Rmax
  priors_sigma = sigma
  
  prior_params = stanvar(priors_pop, 'priors_pop') + 
    stanvar(priors_sdi, 'priors_sdi') + 
    stanvar(priors_sigma, 'priors_sigma')
  
  priors.vi = 
    # Intercept priors (50 % variation around mean value)
    prior(normal(priors_pop[1], priors_sd[1]), nlpar = Rmin, class = b, lb = 0, ub = mu_pars[2]) +
    prior(normal(4.648762, 2.324381), nlpar = Rmax, class = b, lb = 0, ub = 4.671519) +
    prior(normal(0, .5), nlpar = beta, class = b) +
    prior(normal(50, 25), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
    # Random effects priors
    prior(exponential(1), nlpar = NEC, class = sd, group = ID) +
    prior(exponential(1), nlpar = Rmax, class = sd, group = ID) +
    # Residual prior
    prior(exponential(10), class = sigma)
}
  
  prior(normal(0,1), class = "b") +
  prior(normal(5,1), class = "Intercept") +
  prior(normal(0,5), class = "sigma") +
  prior(normal(0,0.75), class = "sd")


# Set target-specific options such as packages:
# tar_option_set(packages = "utils") # nolint

# Build the targets ----
# End this file with a list of target objects.
list(
  tar_target(data, data.frame(x = sample.int(100), y = sample.int(100))),
  tar_target(data_summary, summarize_data(data)) # Call your custom functions.
)
