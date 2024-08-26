# Test functions
library(tidyverse)
library(here)
library(truncnorm)
library(patchwork)
library(brms)
library(SBC)

theme_set(theme_bw(12))

source(here("R/funs/script_funs.R"))

# 1. Testing data_generator() function (passed) ----
# Global test parameters
I = 1000
CVi = .5
nsim = 100
## Test vi_Rmax ----
# Test case with extreme individual variation
# Dose = seq(0,100, length.out = 200)
df = data_generator(ID_pars = "Rmax", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig.Rmax = df %>% 
  ggplot(aes(y = y, x = Dose)) +
  # Add individual trends
  map(1:I, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      # Dose = seq(0, 100, length.out = 200),
      Rmin = log(exp(vars$b_Rmin_Intercept)),
      Rmax = log(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))[.x],
      beta = (vars$b_beta_Intercept + vars$r_ID__beta_Intercept)[.x],
      NEC =( vars$b_NEC_Intercept + vars$r_ID__NEC_Intercept)[.x]),
    alpha = 0.2, 
    linewidth = .25)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", linewidth = 1.25) +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none") +
  ggtitle("Parameter Rmax", 
          subtitle = "with sd = 0.5")
# fig.Rmax
  
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

fig.NEC = df %>% 
  ggplot(aes(y = y, x = Dose)) +
  # Add individual trends
  map(1:I, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      # Dose = seq(0, 100, length.out = 200),
      Rmin = log(exp(vars$b_Rmin_Intercept)),
      Rmax = log(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))[.x],
      beta = (vars$b_beta_Intercept + vars$r_ID__beta_Intercept)[.x],
      NEC =( vars$b_NEC_Intercept + vars$r_ID__NEC_Intercept)[.x]),
    alpha = 0.2, 
    linewidth = .25)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", linewidth = 1.25) +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none") +
  ggtitle("Parameter NEC", 
          subtitle = "with sd = 0.5")
# fig.NEC

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

fig.beta = 
  df %>% 
  ggplot(aes(y = y, x = Dose)) +
  # Add individual trends
  map(1:I, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      # Dose = seq(0, 100, length.out = 200),
      Rmin = log(exp(vars$b_Rmin_Intercept)),
      Rmax = log(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))[.x],
      beta = (vars$b_beta_Intercept + vars$r_ID__beta_Intercept)[.x],
      NEC =( vars$b_NEC_Intercept + vars$r_ID__NEC_Intercept)[.x]),
    alpha = 0.2, 
    linewidth = .25)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", linewidth = 1.25) +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none") +
  ggtitle(expression(paste("Parameter ", beta)), 
          subtitle = "with sd = 0.5")
# fig.beta

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

fig.all = df %>% 
  ggplot(aes(y = y, x = Dose)) +
  # Add individual trends
  map(1:I, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      # Dose = seq(0, 100, length.out = 200),
      Rmin = log(exp(vars$b_Rmin_Intercept)),
      Rmax = log(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))[.x],
      beta = (vars$b_beta_Intercept + vars$r_ID__beta_Intercept)[.x],
      NEC =( vars$b_NEC_Intercept + vars$r_ID__NEC_Intercept)[.x]),
    alpha = 0.2, 
    linewidth = .25)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", linewidth = 1.25) +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none") +
  ggtitle("All parameters", 
          subtitle = "with sd = 0.5")
# fig.all

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



## Recap figure ----
fig.params = (fig.Rmax + fig.NEC + fig.beta ) / fig.all
# fig.params
fig.params.2 = (fig.Rmax + fig.NEC) / (fig.beta + fig.all)
# fig.params.2
fig.params.3 = fig.params.2 & ylim(0, 510)
# fig.params.3

ggsave("outputs/figs/fig.param.sim.jpeg",
       fig.params, 
       width = 18, 
       scale = 1.2,
       units = "cm")

ggsave("outputs/figs/fig.param.sim.2.jpeg",
       fig.params.2, 
       width = 18, 
       scale = 1.2,
       units = "cm")

ggsave("outputs/figs/fig.param.sim.3.jpeg",
       fig.params.3, 
       width = 18, 
       scale = 1.2,
       units = "cm")

# 2. Testing dataset storage function ----
# n_sims_generator = SBC_generator_function(
#   data_generator,
#   D = seq(0, 100, by = 10),
#   I = 10,
#   CV = .05,
#   CVi = .1,
#   mu_pars = list(Rmin = log(1),
#                  Rmax = log(100),
#                  beta = -4.5,
#                  NEC = 10),
#   ID_pars = c("Rmax", "NEC", "beta", "all", "cov"))
n_sim = 2
data_list = generate_datasets(
  n_sims_generator, 
  n_sims = n_sim)


# 3. Apply data generator to list of cases ----
n_sim = 2
n_id = 5
n_cases = 2 # 2 parameters to loop over: Rmax & NEC
n_cvi = 2 # 2 CVi values to loop over: 0 and 0.1
dlist = vector("list", n_cases)

# dlist = list( 
#   Rmax = list(
#     CVi.0 = list(),
#     CVi.0.1 = list()),
#   NEC = list(
#     CVi.0 = list(),
#     CVi.0.1 = list()))

case_list = list(n_sim = 1:n_sim,
                 CVi = c(0, .1, .2, .3, .4, .5),
                 case = c("Rmax", "NEC", "beta", "all"), # "cov" not yet implemented
                 mod = c("pop", "vi")) 

for (i in 1:length(dlist)) {
  # loop over parameter cases
  dlist[[i]] = generate_datasets(
    SBC_generator_function(
      data_generator,
      D = seq(0, 100, by = 10),
      I = n_id,
      CV = .05,
      CVi = case_list$CVi[1],
      mu_pars = list(Rmin = log(1),
                     Rmax = log(100),
                     beta = -4.5,
                     NEC = 10),
      ID_pars = case_list$case[i]), 
    n_sims = n_sim)}
   
  for (j in 1:length(dlist[[1]])) {
    # loop over cases
    dlist[[i]][j] = 
      generate_datasets(
        SBC_generator_function(
          data_generator,
          D = seq(0, 100, by = 10),
          I = n_id,
          CV = .05,
          CVi = case_list$CVi[i],
          mu_pars = list(Rmin = log(1),
                         Rmax = log(100),
                         beta = -4.5,
                         NEC = 10),
          ID_pars = case_list$case[j]), 
        n_sims = n_sim)
  }
}

                                      
for (i in 1:length(case_list$CVi)) {
  data_list_Rmax[i] = generate_datasets(
    SBC_generator_function(
      data_generator,
      D = seq(0, 100, by = 10),
      I = n_id,
      CV = .05,
      CVi = case_list$CVi[i],
      mu_pars = list(Rmin = log(1),
                     Rmax = log(100),
                     beta = -4.5,
                     NEC = 10),
      ID_pars = case_list$case[1]), 
    n_sims = n_sim)
  names(data_list_Rmax[i]) = paste("CVi =", case_list$CVi[i])
}