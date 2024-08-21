# Code for Dose-Response simulations incorporating VI
library(tidyverse); library(viridis); library(brms); library(tidybayes)
library(patchwork); library(truncnorm)
bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

# Functions for linear dose-response with NEC ----
# y ~ bot + (top - bot) * exp(-exp(beta) * (x - nec) * step(x - nec)) 

# Dose-response with NEC and exponential decay bounded between max and min response values
DR_RN_fun = function(Dose, Rmin, Rmax, beta, NEC){
  yhat = Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)) 
  return(yhat)
}

ECx_RN_fun = function(Rmin, Rmax, beta, NEC, Rx){
  ECx = NEC - 1 / exp(beta) * log((Rx - Rmin) / (Rmax - Rmin)) 
  return(ECx)
}

DR_logRN_fun = function(Dose, Rmin, Rmax, beta, NEC){
  logyhat = Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)) 
  yhat = exp(logyhat)
  return(yhat)
}

# TO DEFINE
# ECx_logRN_fun = function(Rmin, Rmax, beta, NEC, Rx){
#   ECx = NEC - 1 / exp(beta) * log((Rx - Rmin) / (Rmax - Rmin)) 
#   return(ECx)
# }


# Test functions: population trend ----
Rmin = 0
Rmax = 100
beta = -4.5
NEC = 10
Dose = seq(0, 100, by = 10)
sigma = .1
nreps = 10

# EC10 = ECx_nlin_fun(Rmin, Rmax, beta, NEC, Rx = 90)
# EC50 = ECx_nlin_fun(Rmin, Rmax, beta, NEC, Rx = 50)
# EC90 = ECx_nlin_fun(Rmin, Rmax, beta, NEC, Rx = 10)

plot(Dose, DR_RN_fun(Dose, Rmin, Rmax, beta, NEC), type = "l", 
     lwd = 3, ylim = c(0, 100))

plot(Dose, DR_logRN_fun(Dose, Rmin, log(Rmax), beta, NEC), type = "l", 
     lwd = 3, ylim = c(0, 100))

set.seed(42)
df.pop = crossing(rep = 1:nreps, Dose) %>% 
  mutate(yhat = DR_RN_fun(Dose, Rmin, Rmax, beta, NEC)) %>% 
  mutate(y = rnorm(n(), yhat, sigma))
df.pop %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_RN_fun, 
                args = list(
                  Rmin = Rmin,
                  Rmax = Rmax,
                  beta = beta,
                  NEC = NEC),
                color = "red", size = 1) +
  ylim(0, 110) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)

# lognormal formulation
set.seed(42)
df.pop = crossing(rep = 1:nreps, Dose) %>% 
  mutate(yhat = DR_logRN_fun(Dose, Rmin, log(Rmax), beta, NEC)) %>% 
  mutate(y = rlnorm(n(), log(yhat), sigma))
df.pop %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_logRN_fun, 
                args = list(
                  Rmin = Rmin,
                  Rmax = log(Rmax),
                  beta = beta,
                  NEC = NEC),
                color = "red", size = 1) +
  ylim(0, 110) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)
df.pop %>% 
  ggplot(aes(y = log(y), x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  labs(x = "Dose", y = "log(y) response") +
  theme_bw(14)

## Prior predictive checks on population trend ----
### Define formula ---- 
bf.pop = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
            Rmin + Rmax + beta + NEC ~ 1, 
            family = lognormal,
            nl = T)

### Set priors ----
lb_R = quantile(df.pop$y, probs = .1)
ub_R = quantile(df.pop$y, probs = .9)
sd_y = sd(df.pop$y)
sd_Dose = sd(df.pop$Dose)
sd_y_prior = 2.5 * sd_y
sd_Dose_prior = 2.5 * sd_Dose

priors.pop = 
  # Intercept priors (50 % variation around mean value)
  prior(normal(1.824934, 0.912467), nlpar = Rmin, class = b, lb = 0, ub = 4.671519) +
  prior(normal(4.648762, 2.324381), nlpar = Rmax, class = b, lb = 0, ub = 4.671519) +
  prior(normal(0, .5), nlpar = beta, class = b) +
  prior(normal(50, 25), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
  prior(exponential(10), class = sigma)

# Plot priors
priors.pop %>% 
  parse_dist() %>% 
  # filter(class == "b") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Intercepts") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 
### Fit to prior  ----
brm.pop.prior = brm(data = df.pop, 
                    bf.pop, 
                    backend = "cmdstan",
                    prior = priors.pop, 
                    sample_prior = "only", 
                    file_refit = "always",
                    save_pars = save_pars(all = TRUE),
                    seed = 42)
brm.pop.prior
pp_check(brm.pop.prior, ndraws = 100)
conditional_effects(brm.pop.prior)
conditional_effects(brm.pop.prior, spaghetti = T, ndraws = 200)

## Posterior predictive checks ----
brm.pop = brm(data = df.pop, 
              bf.pop,
              backend = "cmdstan",
              prior = priors.pop,
              file_refit = "always",
              save_pars = save_pars(all = TRUE),
              warmup = 4000,
              iter = 5000,
              seed = 42,
              cores = 4,
              threads = 3,
              control = list(adapt_delta = .95,
                             max_treedepth = 12),
              stan_model_args=list(stanc_options = list("O1")))
brm.pop
pp_check(brm.pop, ndraws = 100)
plot(conditional_effects(brm.pop), points = T)
conditional_effects(brm.pop, spaghetti = T, ndraws = 200)

# Simulate individual differences in NEC ----
Rmin = 1
Rmax = 100
beta = -4.5
Dose = seq(0, 100, by = 10)
NEC = 10
sigma = .1

n_id = 20
CV_Rmax  = .1 # 10 % variation around average (on log scale?)
CV_NEC  = .75 # 75 % variation around average (on log scale?)
sigma_i_Rmax = CV_Rmax * Rmax
sigma_i_NEC = CV_NEC * NEC

set.seed(42)
Rmax_i = rtruncnorm(n_id, mean = Rmax, sd = sigma_i_Rmax, a = 0)
NEC_i = rtruncnorm(n_id,  mean = NEC, sd = sigma_i_NEC, a = 0)
ID = data.frame(NEC_i = NEC_i,
                Rmax_i = Rmax_i) %>% 
  arrange(NEC_i) %>% 
  mutate(ID = 1:n_id)

# lognormal formulation
set.seed(42)
df.vi = ID %>%
  expand(nesting(ID, Rmax_i, NEC_i), 
         Dose = Dose) %>%
  mutate(yhat = DR_logRN_fun(Dose, log(Rmin), log(Rmax_i), beta, NEC_i)) %>% 
  mutate(y = rlnorm(n(), log(yhat), sigma))

fig = df.vi %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  # Add individual trends
  map(1:n_id, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      Rmax = log(Rmin),
      Rmax = log(ID$Rmax_i[.x]),
      beta = beta,
      NEC = ID$NEC_i[.x]),
    alpha = 0.5, 
    size = .8)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(Rmin), 
                            Rmax = log(Rmax), 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 3) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig


## Prior predictive checks ----
### Define formula ---- 
bf.vi = bf(y  ~ Rmin + (Rmax - Rmin) * exp(-exp(beta) * (Dose - NEC) * (Dose > NEC)), 
           Rmin  + beta ~ 1,
           Rmax + NEC ~ 1 + (1|ID),
            family = lognormal,
            nl = T)

### Set priors ----
lb_R = quantile(df.vi$y, probs = .1)
ub_R = quantile(df.vi$y, probs = .9)
sd_y = sd(df.vi$y)
sd_Dose = sd(df.vi$Dose)
sd_y_prior = 2.5 * sd_y
sd_Dose_prior = 2.5 * sd_Dose

priors.vi = 
  # Intercept priors (50 % variation around mean value)
  prior(normal(1.824934, 0.912467), nlpar = Rmin, class = b, lb = 0, ub = 4.671519) +
  prior(normal(4.648762, 2.324381), nlpar = Rmax, class = b, lb = 0, ub = 4.671519) +
  prior(normal(0, .5), nlpar = beta, class = b) +
  prior(normal(50, 25), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
  # Random effects priors
  prior(exponential(1), nlpar = NEC, class = sd, group = ID) +
  prior(exponential(1), nlpar = Rmax, class = sd, group = ID) +
  # Residual prior
  prior(exponential(10), class = sigma)

# Plot priors
priors.vi %>% 
  parse_dist() %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

### Fit to prior  ----
brm.vi.prior = brm(data = df.vi, 
                    bf.vi, 
                    backend = "cmdstan",
                    prior = priors.vi, 
                    sample_prior = "only", 
                    file_refit = "always",
                    save_pars = save_pars(all = TRUE),
                    seed = 42)
brm.vi.prior
pp_check(brm.vi.prior, ndraws = 100)
conditional_effects(brm.vi.prior)
conditional_effects(brm.vi.prior, spaghetti = T, ndraws = 200)


## Posterior predictive checks ----
brm.vi = brm(data = df.vi, 
              bf.vi,
              backend = "cmdstan",
              prior = priors.vi,
              file_refit = "always",
              save_pars = save_pars(all = TRUE),
              warmup = 4000,
              iter = 5000,
              seed = 42,
              cores = 4,
              threads = 3,
              control = list(adapt_delta = .95,
                             max_treedepth = 12),
              stan_model_args=list(stanc_options = list("O1")))
brm.vi
pp_check(brm.vi, ndraws = 100)
plot(conditional_effects(brm.vi), points = T)
plot(conditional_effects(brm.vi, re_formula = NULL), points = T)
conditional_effects(brm.vi, re_formula = NULL, 
                    spaghetti = T, ndraws = 200)

re = crossing(Dose = seq(min(df.vi$Dose), 
                         max(df.vi$Dose),
                         length.out=100),
              ID = unique(df.vi$ID)) %>% 
  add_epred_draws(brm.vi, re_formula = NULL, 
                  scale = "response", ndraws = 20)

re %>% 
  ggplot(aes(y = .epred, x = Dose)) +
  geom_line(aes(y = .epred, x = Dose, group = .draw), size = .5, alpha = .5) +
  geom_point(data = df.vi, aes(y = y, x = Dose, color = ID)) +
  facet_wrap(~ID, nrow = 6, ncol = 5) + 
  scale_color_viridis() +
  ylab("y response") + 
  xlab("Dose") +
  theme_bw(12) +
  theme(legend.position = "none")

re = brm.vi %>%
  spread_draws(# Population values
    #b_Rmin_Intercept, 
    b_Rmax_Intercept,
    #b_beta_Intercept, 
    b_NEC_Intercept, 
    # Individual offsets
    #r_ID__Rmin[ID,Intercept], 
    r_ID__Rmax[ID,Intercept], 
    #r_ID__beta[ID,Intercept], 
    r_ID__NEC[ID,Intercept],
    # Individual variances
    #sd_ID__Rmin_Intercept, 
    sd_ID__Rmax_Intercept,
    #sd_ID__beta_Intercept, 
    sd_ID__NEC_Intercept,
    sigma) %>% 
  # Individual offsets converted onto the original length scale (in micrometers)
  mutate(#Rmin_i = (b_Rmin_Intercept + r_ID__Rmin),
         Rmax_i = (b_Rmax_Intercept + r_ID__Rmax),
         #beta_i = (b_beta_Intercept + r_ID__beta),
         NEC_i = b_NEC_Intercept + r_ID__NEC) %>% 
  # Population averge distribution
  mutate(#Rmin_dist = rnorm(n(), b_Rmin_Intercept, sd_ID__Rmin_Intercept),
         Rmax_dist = rtruncnorm(n(), 
                                mean = b_Rmax_Intercept, 
                                sd = sd_ID__Rmax_Intercept, 
                                a = 0),
         #beta_dist = rnorm(n(), b_beta_Intercept, sd_ID__beta_Intercept),
         NEC_dist = rtruncnorm(n(), mean = b_NEC_Intercept, 
                               sd = sd_ID__NEC_Intercept, 
                               a = 0))

re.mean = re %>% 
  select(.chain, .iteration, .draw, Rmax_i, NEC_i) %>% 
  mean_qi(Rmax_i, NEC_i)
# Summarize individual values into mean, lower and upper 95 % quantiles

Rmax_dist = re %>% 
  ggplot(aes(x = Rmax_dist)) +
  stat_histinterval(slab_color = "gray45", 
                    outline_bars = TRUE) +
  labs(x = "", y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)
NEC_dist = re %>% 
  ggplot(aes(x = NEC_dist)) +
  stat_histinterval(slab_color = "gray45", 
                    outline_bars = TRUE) +
  labs(x = expression(NEC), y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)
