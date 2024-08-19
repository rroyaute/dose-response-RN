# Code for Dose-Response simulations incorporating VI
library(tidyverse); library(viridis); library(brms); library(tidybayes)
library(patchwork); library(bayesnec)

bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

# Functions for linear dose-response with NEC ----
DR_lin_fun = function(Dose, Rmax, beta, NEC){
  yhat = Rmax - exp(beta) * (Dose - NEC) * (Dose > NEC)
  return(yhat)
}

ECx_lin_fun = function(Rmax, beta, NEC, Rx){
  ECx = (Rmax - Rx) / exp(beta) + NEC
  return(ECx)
}

# Test functions: population trend ----
Rmax = 100
beta = .05
Dose = seq(0, 100, by = 10)
NEC = 10
sigma = 3
nreps = 10

EC10 = ECx_lin_fun(Rmax, beta, NEC, Rx = 90)
EC50 = ECx_lin_fun(Rmax, beta, NEC, Rx = 50)
EC90 = ECx_lin_fun(Rmax, beta, NEC, Rx = 10)


plot(Dose, DR_lin_fun(Dose, Rmax, beta, NEC), type = "l", 
     lwd = 3, ylim = c(0, 100))

set.seed(42)
df = crossing(rep = 1:nreps, Dose) %>% 
  mutate(yhat = DR_lin_fun(Dose, Rmax, beta, NEC)) %>% 
  mutate(y = rnorm(n(), yhat, sigma))
df %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_lin_fun, 
                args = list(
                  Rmax = Rmax,
                  beta = beta,
                  NEC = NEC),
                color = "red", size = 1) +
  ylim(0, 110) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)


# Simulate individual differences in NEC ----
Rmax = 100
beta = .05
Dose = seq(0, 100, by = 10)
NEC = 10
sigma = 3

n_id = 20
CV = .1 # 10 % variation around average
sigma_i_Rmax = CV * Rmax
sigma_i_NEC = CV * NEC


set.seed(42)
Rmax_i = rnorm(n_id, Rmax, sigma_i_Rmax)
NEC_i = rnorm(n_id, NEC, sigma_i_NEC)
ID = data.frame(NEC_i = NEC_i,
                Rmax_i = Rmax_i) %>% 
  arrange(NEC_i) %>% 
  mutate(EC10_i = ECx_lin_fun(Rmax_i, beta, NEC_i, Rx = 90),
         EC50_i = ECx_lin_fun(Rmax_i, beta, NEC_i, Rx = 50),
         EC90_i = ECx_lin_fun(Rmax_i, beta, NEC_i, Rx = 10),
         ID = 1:n_id)

df = ID %>%
  expand(nesting(ID, Rmax_i, NEC_i), 
         Dose = Dose) %>%
  mutate(yhat = DR_lin_fun(Dose, Rmax_i, beta, NEC_i)) %>% 
  mutate(y = rnorm(n(), yhat, sigma)) %>% 
  mutate(yhat = ifelse(y >= 0, y, 0),
         y = ifelse(y >= 0, y, 0)) # make negative values 0

fig = df %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  # Add individual trends
  map(1:n_id, ~geom_function(
    fun = DR_lin_fun,
    args = list(
      Rmax = ID$Rmax_i[.x],
      beta = beta,
      NEC = ID$NEC_i[.x]),
    alpha = 0.5, 
    size = 1)) +
  # Add population average trend
  geom_function(fun = DR_lin_fun, 
                args = list(Rmax = Rmax, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 3) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig


# Fit the model to brms ----
## Define formulas ---- 
bf.pop = bf(y ~ Rmax - exp(beta) * (Dose - NEC) * (Dose > NEC), 
            Rmax + beta + NEC ~ 1, 
            nl = T)
bf.vi.nocorr = bf(y ~ Rmax - exp(beta) * (Dose - NEC) * (Dose > NEC), 
                  Rmax + NEC ~  1 + (1|ID),
                  beta ~ 1,
           nl = T)

## Set priors ----
ub_Rmax = max(df$y)
sd_y = sd(df$y)
sd_Dose = sd(df$Dose)

sd_y_prior = 2.5 * sd_y
sd_Dose_prior = 2.5 * sd_Dose
priors.pop = 
  # Intercept priors
  prior(normal(100, 84.5603), nlpar = Rmax, class = b, lb = 0, ub = 150.1633) +
  prior(normal(0, 5), nlpar = beta, class = b) +
  # prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  prior(normal(50, 79.23723), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
  prior(exponential(1), class = sigma)

priors.vi = 
  # Intercept priors
  prior(normal(100, 84.5603), nlpar = Rmax, class = b, lb = 0, ub = 150.1633) +
  prior(normal(0, 5), nlpar = beta, class = b) +
  # prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  prior(normal(50, 79.23723), nlpar = NEC, class = b, lb = 0, ub = 100) + # Residual prior
  # Random effects priors
  prior(exponential(1), nlpar = NEC, class = sd, group = ID) +
  # Random effects priors
  prior(exponential(1), nlpar = Rmax, class = sd, group = ID) +
  # Residual prior
  prior(exponential(1), class = sigma)


# Plot priors
priors.pop %>% 
  parse_dist() %>% 
  filter(class == "b") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Intercepts") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

# Fit to prior  ----
brm.pop.prior = brm(data = df, 
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
conditional_effects(brm.pop.prior, spaghetti = T, ndraws = 50)

brm.vi.prior = brm(data = df, 
                   bf.vi.nocorr, 
                   backend = "cmdstan",
                   prior = priors.vi, 
                   sample_prior = "only", 
                   file_refit = "always",
                   save_pars = save_pars(all = TRUE),
                   seed = 42)

brm.vi.prior
pp_check(brm.vi.prior, ndraws = 100)
conditional_effects(brm.vi.prior)
conditional_effects(brm.vi.prior, spaghetti = T, ndraws = 50)


# Fit to data ----
brm.pop = brm(data = df, 
                   bf.pop, 
                   backend = "cmdstan",
                   prior = priors.pop, 
                   file_refit = "always",
                   save_pars = save_pars(all = TRUE),
                   seed = 42)

brm.pop
pp_check(brm.pop, ndraws = 100)
conditional_effects(brm.pop)
conditional_effects(brm.pop, spaghetti = T, ndraws = 50)

brm.vi = brm(data = df, 
              bf.vi.nocorr, 
              backend = "cmdstan",
              prior = priors.vi, 
              file_refit = "always",
              save_pars = save_pars(all = TRUE),
              seed = 42)

brm.vi
pp_check(brm.vi, ndraws = 100)
conditional_effects(brm.vi, re_formula = NULL)
conditional_effects(brm.vi, spaghetti = T, ndraws = 50, re_formula = NULL)


# bayesnec checks ----
df$ID = as.factor(df$ID)

bnec.test <- bnec(y ~ crf(Dose, model = "nec4param"),
              data = df, backend = "cmdstan",
              cores = 4, threads = 3)
autoplot(bnec.test)
check_priors(bnec.test)
pull_prior(bnec.test)

bnec.test.vi <- bnec(y ~ crf(Dose, model = "nec4param") +
                       (top + nec | ID),
                  data = df, backend = "cmdstan",
                  cores = 4, threads = 3)
autoplot(bnec.test.vi)
check_priors(bnec.test.vi)
pull_prior(bnec.test.vi)

new_prior = pull_prior(bnec.test.vi)

