# Test functions
library(tidyverse)
library(here)
library(truncnorm)

source(here("R/funs/script_funs.R"))

# 1. Testing data_generator() function (passed)
# Global test parameters
I = 1000
CVi = .5
nsim = 100
## Test vi_Rmax ----
# Test case with extreme individual variation
df = data_generator(ID_pars = "Rmax", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig = df %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  geom_line() +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", size = 2) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig
  
hist(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
mu_Rmax = mean(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
sd_Rmax = sd(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
CV = sd_Rmax / mu_Rmax * 100
CV

## Test vi_NEC ----
# Test case with extreme individual variation
df = data_generator(ID_pars = "NEC", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig = df %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  geom_line() +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", size = 2) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig

hist(vars$b_NEC_Intercept + vars$r_ID__NEC)
mu_NEC = mean(vars$b_NEC_Intercept + vars$r_ID__NEC)
sd_NEC = sd(vars$b_NEC_Intercept + vars$r_ID__NEC)
CV = sd_NEC / mu_NEC * 100
CV


## Test vi_beta ----
# Test case with extreme individual variation
df = data_generator(ID_pars = "beta", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig = df %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  geom_line() +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", size = 2) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig

hist(exp(vars$b_beta_Intercept + vars$r_ID__beta))
mu_beta = mean(exp(vars$b_beta_Intercept + vars$r_ID__beta))
sd_beta = sd(exp(vars$b_beta_Intercept + vars$r_ID__beta))
CV = sd_beta / mu_beta * 100
CV


## Test vi_all ----
# Test case with extreme individual variation
df = data_generator(ID_pars = "all", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig = df %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  geom_line() +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", size = 2) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig


hist(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
mu_Rmax = mean(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
sd_Rmax = sd(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
CV = sd_Rmax / mu_Rmax * 100
CV

hist(vars$b_NEC_Intercept + vars$r_ID__NEC)
mu_NEC = mean(vars$b_NEC_Intercept + vars$r_ID__NEC)
sd_NEC = sd(vars$b_NEC_Intercept + vars$r_ID__NEC)
CV = sd_NEC / mu_NEC * 100
CV

hist(exp(vars$b_beta_Intercept + vars$r_ID__beta))
mu_beta = mean(exp(vars$b_beta_Intercept + vars$r_ID__beta))
sd_beta = sd(exp(vars$b_beta_Intercept + vars$r_ID__beta))
CV = sd_beta / mu_beta * 100
CV




# 2. Testing n_sims_generator
datasets_func = generate_datasets(n_sims_generator, n_sims = nsim) # passed

