library(tidyverse); library(viridis); library(brms); library(tidybayes)
library(patchwork); library(SBC)
bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

# Simulate data
D = seq(0, 100, by = 10)
I = 10
N = length(Dose) * ID

sim_generator = function(N, I) {
  # N - number of datapoints, I number of individuals 
  D = seq(0, 100, by = 10)
  Dose = rep(D, I)
  
  ID = rep(1:I, each = length(D))
  
  b_alpha_Intercept = rnorm(1, 100, 20)
  b_beta_Intercept = rnorm(0, 1)
  b_NEC_Intercept = runif(1, 0, 100)
  sigma = rexp(1, 1)
  sd_ID__alpha_Intercept = rexp(1, 1)
  sd_ID__beta_Intercept = rexp(1, 1)
  sd_ID__NEC_Intercept = rexp(1, 1)
  cor_ID__alpha_Intercept__beta_Intercept = runif(1, -1, 1)
  cor_ID__alpha_Intercept__NEC_Intercept = runif(1, -1, 1)
  cor_ID__beta_Intercept__NEC_Intercept = runif(1, -1, 1)
  
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
  
  yhat = log(b_alpha_Intercept + r_ID__alpha[ID]) - 
    (b_beta_Intercept + r_ID__beta_Intercept[ID]) * 
    (Dose - (b_NEC_Intercept + r_ID__NEC_Intercept[ID])) * 
    step((Dose - NEC))
  y = rnorm(N, yhat, sigma)
  
  list(
    variables = list(
      b_alpha_Intercept = b_alpha_Intercept,
      b_beta_Intercept = b_beta_Intercept,
      b_NEC_Intercept = b_NEC_Intercept,
      sd_ID__alpha_Intercept = sd_ID__alpha_Intercept,
      sd_ID__beta_Intercept = sd_ID__beta_Intercept,
      sd_ID__NEC_Intercept = sd_ID__NEC_Intercept,
      
      
      
      r_group = r_group,
      sigma = sigma
    ),
    generated = data.frame(y = y, x = x, group = group)
  )
}

n_sims_generator = SBC_generator_function(sim_generator, N = 18, ID = 5)

set.seed(12239755)
datasets_func = generate_datasets(n_sims_generator, 1)
