# Code for Dose-Response simulations incorporating VI
library(tidyverse); library(viridis); library(brms); library(tidybayes)
library(patchwork)
bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

# 1. Simulation on the Normal scale ----
# Parameter list
alpha = 100
beta = .08
Dose = seq(0, 100, by = 10)
NEC = 25
sigma = 1
nreps = 10

DR_fun = function(Dose, alpha , beta, NEC){
  yhat = ifelse(Dose < NEC, alpha, alpha * exp(- beta * (Dose - NEC)))
  return(yhat)
}

plot(Dose, DR_fun(Dose, alpha , beta, NEC))

set.seed(42)
df = crossing(rep = 1:nreps, Dose) %>% 
  mutate(yhat = DR_fun(Dose, alpha, beta, NEC)) %>% 
  mutate(y = rnorm(n(), yhat, sigma))
df %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_fun, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 1) +
  ylim(0, 100) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)

# 2. Simulations using a lognormal likelihood ----
DR_fun_log = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - beta * (Dose - NEC))
  return(log_yhat)
}

DR_fun_log_exp = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - beta * (Dose - NEC))
  yhat = exp(log_yhat)
  return(yhat)
}

plot(Dose, exp(DR_fun_log(Dose, alpha , beta, NEC)))

sigma = .08

set.seed(42)
df = crossing(rep = 1:nreps, Dose) %>% 
  mutate(log_yhat = DR_fun_log(Dose, alpha, beta, NEC)) %>% 
  mutate(y = rlnorm(n(), log_yhat, sigma))
df %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_fun_log_exp, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 1) +
  # ylim(0, 100) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)

#  3. Simulate individuals with different sensitivities ----
Dose = seq(0, 100, by = 10)
n_id = 20
sigma = .1
rho_1_2 = -.5 # Individual with stronger initial response have lower decay rate
rho_1_3 = .5 # Individual with stronger initial response have higher NEC
rho_2_3 = -.5 # Individual with decay rate response lower NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * .1, beta * .1, NEC * .1) # 10 % CV around the mean
rho_mat = matrix(c(1, rho_1_2, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID = 1:n_id)

# Simulate individual growth
df = ID %>%
  expand(nesting(ID, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>% 
  mutate(y = rlnorm(n(), log_yhat, sigma))

# Sort individuals by alpha_i and create a mapping
ID_sorted <- ID %>%
  arrange(desc(alpha_i)) %>%
  mutate(color_index = row_number())

# Create a viridis color palette for n_id colors
viridis_colors <- viridis(n_id)

fig = df %>% 
  left_join(ID_sorted, by = "ID") %>%  # Join to get the color_index
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  # Add individual trends
  map(1:n_id, ~geom_function(
    fun = DR_fun_log_exp,
    args = list(
      alpha = ID_sorted$alpha_i[.x],
      beta = ID_sorted$beta_i[.x],
      NEC = ID_sorted$NEC_i[.x]
    ),
    color = viridis_colors[ID_sorted$color_index[.x]],  # Use sorted index for colors
    alpha = 0.5, 
    size = 1
  )) +
  # Add population average trend
  geom_function(fun = DR_fun_log_exp, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 3) +
  scale_color_viridis_d() +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig


# 4. Fit the model to brms ----
## 4.1 Define formulas ---- 
bf.pop = bf(y ~ log(alpha) - beta * (Dose - NEC) * (Dose > NEC), 
            alpha + beta + NEC ~ 1, 
            nl = T, 
            family = lognormal)
bf.vi = bf(y ~ log(alpha) - beta * (Dose - NEC) * (Dose > NEC), 
           alpha + beta + NEC ~ 1 + (1|c|ID), 
           nl = T, 
           family = lognormal)
bf.vi.nocorr = bf(y ~ log(alpha) - beta * (Dose - NEC) * (Dose > NEC), 
                  alpha + beta + NEC ~ 1 + (1|ID), 
                  nl = T, 
                  family = lognormal)
bf.vi.NEC = bf(y ~ log(alpha) - beta * (Dose - NEC) * (Dose > NEC), 
               alpha + beta ~ 1,
               NEC ~ 1 + (1|ID), 
               nl = T, 
               family = lognormal)
## 4.2 Set priors ----
priors.pop = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Residual prior
  prior(exponential(50), class = sigma)

priors.vi = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Random effects priors
  prior(exponential(.05), nlpar = alpha, class = sd, group = ID) +
  prior(exponential(5), nlpar = beta, class = sd, group = ID) +
  prior(exponential(.1), nlpar = NEC, class = sd, group = ID) +
  # prior(normal(0, 5), nlpar = beta, class = sd, group = ID, lb = 0) +
  # prior(normal(0, 5), nlpar = alpha, class = sd, group = ID, lb = 0) +
  # prior(normal(0, 10), nlpar = NEC, class = sd, group = ID, lb = 0) +
  # # Residual prior
  prior(exponential(50), class = sigma) +
  prior(lkj(4), class = cor)

priors.vi.2 = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Random effects priors
  prior(exponential(1), nlpar = alpha, class = sd, group = ID) +
  prior(exponential(1), nlpar = beta, class = sd, group = ID) +
  prior(exponential(1), nlpar = NEC, class = sd, group = ID) +
  # Residual prior
  prior(exponential(1), class = sigma) +
  prior(lkj(4), class = cor)

priors.vi.nocorr = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Random effects priors
  prior(exponential(.05), nlpar = alpha, class = sd, group = ID) +
  prior(exponential(5), nlpar = beta, class = sd, group = ID) +
  prior(exponential(.1), nlpar = NEC, class = sd, group = ID) +
  # prior(normal(0, 5), nlpar = beta, class = sd, group = ID, lb = 0) +
  # prior(normal(0, 5), nlpar = alpha, class = sd, group = ID, lb = 0) +
  # prior(normal(0, 10), nlpar = NEC, class = sd, group = ID, lb = 0) +
  # # Residual prior
  prior(exponential(50), class = sigma)
priors.vi.NEC = 
  # Intercept priors
  prior(normal(100, 20), nlpar = alpha, class = b, lb = 0) +
  prior(exponential(10), nlpar = beta, class = b, lb = 0) +
  prior(uniform(0, 100), nlpar = NEC, class = b, lb = 0, ub = 100) + 
  # Random effects priors
  prior(exponential(.1), nlpar = NEC, class = sd, group = ID) +
  # prior(normal(0, 10), nlpar = NEC, class = sd, group = ID, lb = 0) +
  # # Residual prior
  prior(exponential(50), class = sigma)

# Plot priors
p1 = priors.vi %>% 
  parse_dist() %>% 
  filter(class == "b") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Intercepts") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

p2 = priors.vi %>% 
  parse_dist() %>% 
  filter(class == "sd") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Among-individual variance") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

p3 = priors.vi %>% 
  parse_dist() %>% 
  filter(class == "sigma") %>% 
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Residual variance") +
  xlab("Value") + ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90)) 

(p1 + p2 + p3) + plot_layout(ncol = 1)

## 4.3 Prior predictive checks ----
### 4.3.1 With the population model ----
brm.pop.prior= brm(data = df, 
                   bf.pop, 
                   backend = "cmdstan",
                   prior = priors.pop, 
                   sample_prior = "only", 
                   file_refit = "always",
                   save_pars = save_pars(all = TRUE),
                   seed = 42, 
                   file = "mods/brm.pop.prior")

brm.pop.prior
pp_check(brm.pop.prior, ndraws = 100)
conditional_effects(brm.pop.prior)
cond_plot = conditional_effects(brm.pop.prior, ndraws = 100, spaghetti = T, plot = F)

# ggplot(cond_plot[[1]], aes(x = effect1__, y = estimate__)) +
#   geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "gray50", alpha = 0.3) +
#   geom_line(color = "black") +
#   # theme_darkgray() +
#   labs(x = names(cond_plot)[1], y = "Response")

### 4.3.2 With the individual model ----
brm.vi.prior= brm(data = df, 
                  bf.vi, 
                  backend = "cmdstan",
                  prior = priors.vi, 
                  sample_prior = "only", 
                  file_refit = "always",
                  save_pars = save_pars(all = TRUE),
                  seed = 42, 
                  file = "mods/brm.vi.prior")

brm.vi.prior
# pp_check(brm.vi.prior, ndraws = 100)
conditional_effects(brm.vi.prior)

## 4.3.3 With the individual model - no correlation ----
brm.vi.nocorr.prior = brm(data = df, 
                          bf.vi.nocorr, 
                          backend = "cmdstan",
                          prior = priors.vi.nocorr, 
                          sample_prior = "only", 
                          file_refit = "always",
                          save_pars = save_pars(all = TRUE),
                          seed = 42, 
                          cores = 4,
                          threads = 3,
                          control = list(adapt_delta = .95,
                                         max_treedepth = 12), 
                          file = "mods/brm.vi.nocorr.prior")
brm.vi.nocorr.prior
# pp_check(brm.vi.nocorr.prior, ndraws = 100)
conditional_effects(brm.vi.nocorr.prior, re_formula = NULL)


## 4.3.4 With the individual model - NEC only ----
brm.vi.NEC.prior = brm(data = df, 
                       bf.vi.NEC, 
                       backend = "cmdstan",
                       prior = priors.vi.NEC, 
                       sample_prior = "only", 
                       file_refit = "always",
                       save_pars = save_pars(all = TRUE),
                       seed = 42, 
                       cores = 4,
                       threads = 3,
                       control = list(adapt_delta = .95,
                                      max_treedepth = 12), 
                       file = "mods/brm.vi.NEC.prior")
brm.vi.NEC.prior
# pp_check(brm.vi.NEC.prior, ndraws = 100)
conditional_effects(brm.vi.NEC.prior, re_formula = NULL)

## 4.4 Fitting to data ----
## 4.4.1 With the population model ----
brm.pop= brm(data = df, 
             bf.pop, 
             backend = "cmdstan",
             prior = priors.pop, 
             sample_prior = "yes", 
             file_refit = "always",
             save_pars = save_pars(all = TRUE),
             seed = 42, 
             cores = 4,
             threads = 3,
             control = list(adapt_delta = .95,
                            max_treedepth = 12), 
             file = "mods/brm.pop")
brm.pop
pp_check(brm.pop, ndraws = 100)
conditional_effects(brm.pop, ndraws = 100, spaghetti = T)

## 4.4.2 With the individual model ----
brm.vi = brm(data = df, 
             bf.vi, 
             backend = "cmdstan",
             prior = priors.vi, 
             sample_prior = "yes", 
             file_refit = "always",
             save_pars = save_pars(all = TRUE),
             seed = 42, 
             cores = 4,
             threads = 3,
             control = list(adapt_delta = .95,
                            max_treedepth = 12), 
             file = "mods/brm.vi")
brm.vi
pp_check(brm.vi, ndraws = 100)
conditional_effects(brm.vi, re_formula = NULL)
conditional_effects(brm.vi, re_formula = NULL, spaghetti = T, ndraws = 100)


## 4.4.3 With the individual model - no correlation ----
brm.vi.nocorr = brm(data = df, 
                    bf.vi.nocorr, 
                    backend = "cmdstan",
                    prior = priors.vi.nocorr, 
                    sample_prior = "yes", 
                    file_refit = "always",
                    save_pars = save_pars(all = TRUE),
                    seed = 42, 
                    cores = 4,
                    threads = 3,
                    control = list(adapt_delta = .95,
                                   max_treedepth = 12), 
                    file = "mods/brm.vi.nocorr")
brm.vi.nocorr
pp_check(brm.vi.nocorr, ndraws = 100)
conditional_effects(brm.vi.nocorr, re_formula = NULL)

## 4.4.4 With the individual model - NEC only ----
brm.vi.NEC = brm(data = df, 
                       bf.vi.NEC, 
                       backend = "cmdstan",
                       prior = priors.vi.NEC, 
                       sample_prior = "yes", 
                       file_refit = "always",
                       save_pars = save_pars(all = TRUE),
                       seed = 42, 
                       cores = 4,
                       threads = 3,
                       control = list(adapt_delta = .95,
                                      max_treedepth = 12), 
                       file = "mods/brm.vi.NEC")
brm.vi.NEC
pp_check(brm.vi.NEC, ndraws = 100)
conditional_effects(brm.vi.NEC, re_formula = NULL)


# 5. Plot individual tendencies -----
re = crossing(Dose = seq(min(df$Dose), 
                         max(df$Dose),
                         length.out=100),
              ID = unique(df$ID)) %>% 
  add_epred_draws(brm.vi, re_formula = NULL, 
                  scale = "response", ndraws = 20)

re %>% 
  ggplot(aes(y = .epred, x = Dose)) +
  geom_line(aes(y = .epred, x = Dose, group = .draw), size = .5, alpha = .5) +
  geom_point(data = df, aes(y = y, x = Dose, color = ID)) +
  facet_wrap(~ID, nrow = 6, ncol = 5) + 
  scale_color_viridis() +
  ylab("y response") + 
  xlab("Dose") +
  theme_bw(12) +
  theme(legend.position = "none")


re = brm.vi %>%
  spread_draws(# Population values
    b_alpha_Intercept, b_beta_Intercept, b_NEC_Intercept, 
    # Individual offsets
    r_ID__alpha[ID,Intercept], r_ID__beta[ID,Intercept], r_ID__NEC[ID,Intercept],
    # Individual variances
    sd_ID__alpha_Intercept, sd_ID__beta_Intercept, sd_ID__NEC_Intercept,
    sigma) %>% 
  # Individual offsets converted onto the original length scale (in micrometers)
  mutate(alpha_i = (b_alpha_Intercept + r_ID__alpha),
         beta_i = (b_beta_Intercept + r_ID__beta),
         NEC_i = b_alpha_Intercept + r_ID__NEC) %>% 
  # Population averge distribution
  mutate(alpha_dist = rnorm(n(), b_alpha_Intercept, sd_ID__alpha_Intercept),
         beta_dist = rnorm(n(), b_beta_Intercept, sd_ID__beta_Intercept),
         NEC_dist = rnorm(n(), b_NEC_Intercept, sd_ID__NEC_Intercept))
re.mean = re %>% 
  select(.chain, .iteration, .draw, alpha_i, beta_i, NEC_i) %>% 
  mean_qi(alpha_i, beta_i, NEC_i)
# Summarize individual values into mean, lower and upper 95 % quantiles

# Plot population average (diagonal elements)
alpha_dist = re %>% 
  ggplot(aes(x = alpha_dist)) +
  stat_histinterval(slab_color = "gray45", 
                    outline_bars = TRUE) +
  labs(x = "", y = expression(alpha)) +
  theme_bw(12) +
  theme(aspect.ratio=1)
beta_dist = re %>% 
  ggplot(aes(x = beta_dist)) +
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

# Plot individual average with CI (lower diagonal elements)
corr1 = re.mean %>% 
  ggplot(aes(x = alpha_i, y = beta_i)) +
  geom_errorbarh(aes(xmin = alpha_i.lower, xmax = alpha_i.upper)) +
  geom_errorbar(aes(ymin = beta_i.lower, ymax = beta_i.upper)) +
  geom_point(alpha = .8, size = 3) +
  labs(x = "", y = expression(beta)) +
  theme_bw(12) +
  theme(aspect.ratio=1)
corr2 = re.mean %>% 
  ggplot(aes(x = beta_i, y = NEC_i)) +
  geom_errorbarh(aes(xmin = beta_i.lower, xmax = beta_i.upper)) +
  geom_errorbar(aes(ymin = NEC_i.lower, ymax = NEC_i.upper)) +
  geom_point(alpha = .8, size = 3) +
  labs(x = expression(beta), y = expression(NEC)) +
  theme_bw(12) +
  theme(aspect.ratio=1)
corr3 = re.mean %>% 
  ggplot(aes(x = beta_i, y = NEC_i)) +
  geom_errorbarh(aes(xmin = beta_i.lower, xmax = beta_i.upper)) +
  geom_errorbar(aes(ymin = NEC_i.lower, ymax = NEC_i.upper)) +
  geom_point(alpha = .8, size = 3) +
  labs(x = expression(beta), y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)

# Plot correlation estimate (upper diagonal elements)
dcorr1 = brm.vi %>% 
  spread_draws(`cor.*`, regex = TRUE) %>% 
  ggplot(aes(x = cor_ID__alpha_Intercept__beta_Intercept)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linewidth = 1, color = "black", linetype = "dashed") +
  xlim(-1, 1) +
  labs(x = "", y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)
dcorr2 = brm.vi %>% 
  spread_draws(`cor.*`, regex = TRUE) %>% 
  ggplot(aes(x = cor_ID__alpha_Intercept__NEC_Intercept)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linewidth = 1, color = "black", linetype = "dashed") +
  xlim(-1, 1) +
  labs(x = "", y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)
dcorr3 = brm.vi %>% 
  spread_draws(`cor.*`, regex = TRUE) %>% 
  ggplot(aes(x = cor_ID__beta_Intercept__NEC_Intercept)) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linewidth = 1, color = "black", linetype = "dashed") +
  xlim(-1, 1) +
  labs(x = "", y = "") +
  theme_bw(12) +
  theme(aspect.ratio=1)

# Arrange plot into 3 x 3 grid
alpha_dist + dcorr1 + dcorr2 +
  corr1 + beta_dist + dcorr3 +
  corr2 + corr3 + NEC_dist +
  plot_layout(ncol = 3, nrow = 3, byrow = T)
