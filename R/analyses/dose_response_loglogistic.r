library(tidyverse)
library(drc)
library(brms)
library(ggdist)
library(truncnorm)
library(distributional)
library(viridis)
library(patchwork)
library(GGally)

# Dose-response form and parameter values ----
sigma <- .1
c <- 0 # bottom of the dose-response curve
d <- 1 # top of the dose-response curve
b <- 3 # slope of the dose-response curve
e <- .3 # EC50 of the dose-response curve
x <- seq(0, 1, by = .01)

mu <- c + (d - c) / (1 + exp(b * log(x / e)))
mu <- d / (1 + exp(b * log(x / e)))

plot(x, mu, type = "l")

# Figure: Prediction vs. Confidence Interval over EC50 ----
## Simulate data for typical Dose-Response Experiment ----
# 5 doses + 1 control

set.seed(42)
df.sim <- data.frame(x = seq(0, 1, length.out = 6)) %>% # Sample 6 values between [0;1] with equal spacing
  mutate(x = case_when(x == 0 ~ 0.001, .default = x)) %>% # Replace dose = 0 with small value to avoid computational issues
  mutate(mu = d / (1 + exp(b * log(x / e)))) %>% # Apply dose-response equation to all x-values
  mutate(y = rlnorm(n(), log(mu), sigma)) %>% # Sample from log-normal distribution to keep y-values > 0
  mutate(log_y = log(y))

plot(y ~ x, df.sim)

## brms estimation ----

### Define brms model ----
# 3-parameter log-logistic function
formula_drc <- bf(
  log(y) ~ log(d / (1 + exp(b * (log(x) - log(e))))),
  b + d + e ~ 1,
  nl = TRUE
)

priors <- c(
  prior(normal(2, 1), nlpar = "b", lb = 0), # slope, positive for decreasing curve
  prior(normal(0, 0.5), nlpar = "d"), # log(upper asymptote); log(1) = 0 for max y ≈ 1
  prior(normal(0.5, 0.2), nlpar = "e", lb = 0), # EC50 prior
  prior(normal(0, 1), class = "sigma")
)

inits <- list(
  b_b_Intercept = 2,
  b_d_Intercept = 1,
  b_e_Intercept = 0.5
)
init_list <- rep(list(inits), 4)

### Plot priors ----
priors %>%
  parse_dist() %>%
  filter(class == "b") %>%
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~nlpar, scales = "free") +
  ggtitle("Intercepts") +
  xlab("Value") +
  ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90))

### Run model on priors only ----
brm.drc.prior <- brm(
  formula = formula_drc,
  data = df.sim,
  family = "gaussian",
  prior = priors,
  sample_prior = "only",
  init = rep(list(inits), 4),
  chains = 4,
  iter = 2000,
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.drc.prior",
  file_refit = "always"
)

plot(brm.drc.prior)
pp_check(brm.drc.prior, ndraws = 200)
conditional_effects(brm.drc.prior, method = "posterior_predict")
conditional_effects(brm.drc.prior, spaghetti = T, ndraws = 100)

### Fit model to simulated data ----
brm.drc <- brm(
  formula = formula_drc,
  data = df.sim,
  family = "gaussian",
  prior = priors,
  sample_prior = "yes",
  init = rep(list(inits), 4),
  chains = 4,
  iter = 4000,
  control = list(adapt_delta = .95, max_treedepth = 15),
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.drc",
  file_refit = "always"
)

plot(brm.drc)
pp_check(brm.drc, ndraws = 200)
conditional_effects(brm.drc, method = "posterior_predict")
CRI <- as.data.frame(posterior_summary(brm.drc))

# Prediction intervals via posterior_predict
x_new <- data.frame(x = seq(0, 1, by = 0.01)) %>%
  mutate(x = case_when(x == 0 ~ 0.001, .default = x))

pred <- posterior_predict(brm.drc, newdata = x_new) # Prediction for the expected (mean) value for computing 95% CRI
epred <- posterior_epred(brm.drc, newdata = x_new) # Predictions including residual error to compute 95% Prediction Intervals

pred_df <- data.frame(
  x = x_new$x,
  y = apply(pred, 2, function(p) exp(median(p))),
  y_low = apply(pred, 2, function(p) exp(quantile(p, 0.025))),
  y_up = apply(pred, 2, function(p) exp(quantile(p, 0.975)))
)

epred_df <- data.frame(
  x = x_new$x,
  y = apply(epred, 2, function(p) exp(median(p))),
  y_low = apply(epred, 2, function(p) exp(quantile(p, 0.025))),
  y_up = apply(epred, 2, function(p) exp(quantile(p, 0.975)))
)

fig_predinterval <- ggplot(data = df.sim, aes(x = x, y = y)) +
  # Fitted line
  geom_line(
    data = pred_df,
    aes(x = x, y = y),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  # 95 % Credible Interval for the mean
  geom_ribbon(
    data = epred_df,
    aes(ymin = y_low, ymax = y_up),
    alpha = .4,
    fill = "dodgerblue"
  ) +
  # 95 % Prediction Interval
  geom_ribbon(
    data = pred_df,
    aes(ymin = y_low, ymax = y_up),
    alpha = .4,
    fill = "dodgerblue"
  ) +
  geom_point(
    data = CRI,
    aes(x = CRI["b_e_Intercept", 1], y = .5),
    color = "tomato2"
  ) +
  geom_errorbarh(
    aes(
      y = .5,
      x = CRI["b_e_Intercept", 1],
      xmin = CRI["b_e_Intercept", 3],
      xmax = CRI["b_e_Intercept", 4],
      width = 0
    ),
    color = "tomato2"
  ) +
  # Data points
  geom_point(
    shape = 21,
    colour = "black",
    fill = "white",
    size = 2,
    stroke = 1
  ) +
  labs(x = "Dose", y = "Phenotype") +
  theme_bw(12)

ggsave(filename = "outputs/figs/fig_predinterval.jpeg", fig_predinterval)

# Distribution of EC50 values
# Extract posterior samples of e and its variance
post <- as_draws_df(mod.brm)

# Simulate new individual EC50s assuming e ~ lognormal(mu_e, sigma_e)
e_dist <- rlnorm(
  4000,
  meanlog = mean(log(post$b_e_Intercept)),
  sdlog = sd(log(post$b_e_Intercept))
)
hist(e_dist)

# Plot as a distribution
ggplot(data.frame(e = e_dist), aes(x = e)) +
  stat_halfeye(fill = "skyblue") +
  theme_bw()


## drc estimation (TODO) ----
mod.drm <- drm(y ~ x, data = df.sim, fct = LN.3()) # log-normal model with c = 0 and d = 1
mod.drm <- drm(log_y ~ x, data = df.sim, fct = LL.3()) # Log-logistic model on log(y)

summary(mod.drm)
plot(mod.drm)
coefs <- as.numeric(coef(mod.drm))
CI <- confint(mod.drm, "e")
x <- data.frame(x = seq(0, 1, by = .01))
pred <- predict(mod.drm, newdata = x, interval = "prediction")
pred <- data.frame(x = x, y = pred[, 1], y_low = pred[, 2], y_up = pred[, 3])

ggplot(data = df.sim, aes(x = x, y = y)) +
  geom_point(
    shape = 21,
    colour = "black",
    fill = "white",
    size = 2,
    stroke = 1
  ) +
  geom_line(
    data = pred,
    aes(x = x, y = y),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  geom_ribbon(
    data = pred,
    aes(ymin = y_low, ymax = y_up),
    alpha = .4,
    fill = "dodgerblue"
  ) +
  geom_point(aes(x = coefs[3], y = .5), color = "tomato2") +
  geom_errorbarh(
    aes(y = .5, x = coefs[3], xmin = CI[, 1], xmax = CI[, 2], width = 0),
    color = "tomato2"
  ) +
  theme_bw()

## Old code ----
x <- seq(0, 1, by = .15)

mu <- c + (d - c) / (1 + exp(b * log(x / e)))
mu <- d / (1 + exp(b * log(x / e)))
plot(x, mu)

set.seed(42)
df <- data.frame(x = seq(0, 1, by = .01)) %>%
  mutate(mu = d / (1 + exp(b * log(x / e)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma)) %>%
  mutate(y_low = qlnorm(p = .025, log(mu), sigma)) %>%
  mutate(y_up = qlnorm(p = .975, log(mu), sigma))

mod.drm <- drm(y ~ x, data = df, fct = LL.3())
summary(mod.drm)
coefs <- as.numeric(coef(mod.drm))
CI <- confint(mod.drm, "e")


df %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = mu), linewidth = 1, color = "dodgerblue") +
  geom_point() +
  geom_ribbon(aes(ymin = y_low, ymax = y_up), alpha = .4, fill = "dodgerblue") +
  geom_point(aes(x = coefs[3], y = .5), color = "tomato2") +
  geom_errorbarh(
    aes(y = .5, x = coefs[3], xmin = CI[, 1], xmax = CI[, 2], width = 0),
    color = "tomato2"
  ) +
  theme_bw()

# Figure: Genotype sensitivity ----
Dose <- seq(0, 1, length.out = 6)
n_g <- 5 # 5 genotypes
CVa <- .1 # 10 % of variation around the mean for all parameters
sigma_d <- d * CVa # Upper bound variation
sigma_b <- b * CVa # Rate variation
sigma_e <- e * CVa # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters
# Genotypes with higher basal expression are less sensitive
# r_dxe = rho: higher basal expression <-> higher EC50
# r_dxb = -rho : higher basal expression <-> shallower slope
# r_bxe = -rho : higher EC50 <-> shallower slope

set.seed(42)
d_g <- rtruncnorm(n_g, mean = d, sd = sigma_d, a = 0)
b_g <- rtruncnorm(n_g, mean = b, sd = sigma_b, a = 0)
e_g <- rtruncnorm(n_g, mean = e, sd = sigma_e, a = 0)

Mu <- c(d, b, e)
sigmas <- c(sigma_d, sigma_b, sigma_e) # 10 % CV around the mean
names <- c("sigma_d", "sigma_b", "sigma_e")
rho_mat <- matrix(c(1, -rho, rho, -rho, 1, -rho, rho, -rho, 1), nrow = 3) # Correlation matrix
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

set.seed(42)
G <- MASS::mvrnorm(n_g, Mu, Sigma) %>%
  data.frame() %>%
  set_names("d_g", "b_g", "e_g") %>%
  mutate(G = 1:n_g) %>%
  arrange(e_g) %>%
  mutate(color_index = row_number())

# Visualize the among genotype correlations
G %>%
  dplyr::select(d_g:e_g) %>%
  GGally::ggpairs() +
  theme_bw()

# Visualize genotype-specific dose-responses
df.sim.lines <- G %>%
  expand(
    nesting(G, color_index, d_g, b_g, e_g),
    Dose = seq(0, 1, by = 0.01)
  ) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose / e_g))))

fig_genotypes <- ggplot(
  df.sim.lines,
  aes(y = y, x = Dose, color = factor(color_index), fill = factor(color_index))
) +
  geom_line(data = df.sim.lines, aes(y = mu, x = Dose), linewidth = 1) +
  geom_point(
    data = G,
    aes(y = d_g / 2, x = e_g),
    size = 2.5,
    shape = 21,
    fill = "white",
    alpha = .8
  ) +
  scale_color_viridis_d(option = "H", direction = -1, begin = .2, end = .8) +
  # stat_halfeye(aes(y = .5, xdist = dist_normal(e, sigma_e)),
  #              color = "black", fill = "grey", alpha = .6) +
  labs(x = "Dose", y = "Phenotype") +
  theme_bw(12) +
  theme(legend.position = "none")

ggsave(filename = "outputs/figs/fig_genotypes.jpeg", fig_genotypes)

# Figure: Pre-Post exposure reaction norms ----
## Dose-response figure ----
Dose <- seq(0, 1, length.out = 6)
n_id <- 20 # 10 individuals per doses
CVi <- .1 # 10 % of variation around mean for all parameters
sigma_d <- d * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters

set.seed(42)
d_i <- rtruncnorm(n_id, mean = d, sd = sigma_d, a = 0)
b_i <- rtruncnorm(n_id, mean = b, sd = sigma_b, a = 0)
e_i <- rtruncnorm(n_id, mean = e, sd = sigma_e, a = 0)

Mu <- c(d, b, e)
sigmas <- c(sigma_d, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, -rho, rho, -rho, 1, -rho, rho, -rho, 1), nrow = 3)
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

set.seed(42)
ID <- MASS::mvrnorm(n_id * length(Dose[2:6]), Mu, Sigma) %>%
  data.frame() %>%
  set_names("d_i", "b_i", "e_i") %>%
  mutate(ID = 1:(n_id * length(Dose[2:6]))) %>%
  mutate(assigned_dose = rep(Dose[2:6], each = 20)) %>%
  mutate(control_dose = .001) %>%
  mutate(Group = as.factor(assigned_dose)) %>%
  arrange(e_i) %>%
  mutate(color_index = row_number())

# Visualize the among genotype correlations
ID %>%
  dplyr::select(d_i:e_i) %>%
  GGally::ggpairs() +
  theme_bw()

# Visualize genotype-specific dose-responses
df.sim.id <- ID %>%
  pivot_longer(cols = c(assigned_dose:control_dose), values_to = "Dose") %>%
  mutate(mu = d_i / (1 + exp(b_i * log(Dose / e_i)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma))

df.sim.line <- data.frame(Dose = seq(0, 1, by = 0.01)) %>%
  mutate(y = d / (1 + exp(b * log(Dose / e))))

fig_ID <- ggplot(df.sim.id, aes(y = y, x = Dose)) +
  geom_line(linewidth = .5, aes(group = ID), alpha = .15) +
  geom_point(size = 2.5, shape = 21, fill = "white", alpha = .8) +
  geom_line(
    data = df.sim.line,
    aes(y = y, x = Dose),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  labs(x = "Dose", y = "Phenotype") +
  theme_bw(12) +
  theme(legend.position = "none")

ggsave(filename = "outputs/figs/fig_ID.jpeg", fig_ID)


## Among-individual variation in dose-response ----
df.sim.id <- df.sim.id %>%
  mutate(
    Phase = as.factor(case_when(Dose == 0.001 ~ "Pre", .default = "Post"))
  ) %>%
  mutate(Phase = fct_relevel(Phase, "Pre", "Post"))

### Define brms model ----
formula_drc_vi <- bf(
  log(y) ~ log(d / (1 + exp(b * (log(Dose) - log(e))))),
  b + d + e ~ 1 + (1 | c | ID),
  nl = TRUE
)

priors.vi <- c(
  prior(normal(2, 1), nlpar = "b", lb = 0), # slope, positive for decreasing curve
  prior(normal(0, 0.5), nlpar = "d"), # log(upper asymptote); log(1) = 0 for max y ≈ 1
  prior(normal(0.5, 0.2), nlpar = "e", lb = 0), # ED50, unchanged
  prior(exponential(10), class = "sd", nlpar = "b"), # Exponential prior with mean 0.10
  prior(exponential(10), class = "sd", nlpar = "d"),
  prior(exponential(10), class = "sd", nlpar = "e"),
  prior(exponential(10), class = "sigma"),
  prior(lkj_corr_cholesky(2), class = "L")
)

### Plot priors ----
priors.vi %>%
  parse_dist() %>%
  filter(class %in% c("b", "sd", "sigma")) %>%
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~class, scales = "free") +
  xlab("Value") +
  ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90))

### Run model on priors only ----
brm.drc.vi.prior <- brm(
  formula = formula_drc_vi,
  data = df.sim.id,
  family = "gaussian",
  prior = priors.vi,
  sample_prior = "only",
  init = rep(list(inits), 4),
  chains = 4,
  iter = 4000,
  threads = threading(4),
  cores = 4,
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.drc.prior",
  file_refit = "always"
)

# plot(brm.drc.vi.prior)
conditional_effects(brm.drc.vi.prior, spaghetti = T, ndraws = 100)
plot(conditional_effects(brm.drc.vi.prior, re_formula = NULL), points = T)
pp_check(brm.drc.vi.prior, ndraws = 200)

### Fit model to simulated data ----
brm.drc.vi <- brm(
  formula = formula_drc_vi,
  data = df.sim.id,
  family = "gaussian",
  prior = priors.vi,
  sample_prior = "yes",
  init = rep(list(inits), 4),
  chains = 4,
  warmup = 8000,
  iter = 10000,
  control = list(adapt_delta = .95, max_treedepth = 15),
  threads = threading(4),
  cores = 4,
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.drc.vi",
  file_refit = "always"
)

# plot(brm.drc.vi)
plot(
  conditional_effects(
    brm.drc.vi,
    re_formula = NULL,
    method = "posterior_predict"
  ),
  points = T
)
pp_check(brm.drc.vi, ndraws = 200)

## Slope variance estimation and figure ----
### Define brms model ----
formula_lmm_vi <- bf(
  log(y) ~ Dose + (1 + Phase | gr(ID, by = Group))
)
get_prior(formula_lmm_vi, df.sim.id)

### Prior predictive checks ----
priors.vi.lmm <- c(
  prior(normal(0, 1), class = "b"),
  prior(exponential(10), class = "sd"), # Exponential prior with mean 0.10
  prior(exponential(10), class = "sigma"),
  prior(lkj_corr_cholesky(2), class = "L")
)

### Plot priors ----
priors.vi.lmm %>%
  parse_dist() %>%
  filter(class %in% c("b", "sd", "sigma")) %>%
  ggplot(aes(xdist = .dist_obj, y = format(.dist_obj))) +
  stat_dist_halfeye() +
  facet_wrap(~class, scales = "free") +
  xlab("Value") +
  ylab("Density") +
  theme_bw(12) +
  theme(axis.text.y = element_text(angle = 90))

### Run model on priors only ----
brm.lmm.vi.prior <- brm(
  formula = formula_lmm_vi,
  data = df.sim.id,
  family = "gaussian",
  prior = priors.vi.lmm,
  sample_prior = "only",
  chains = 4,
  iter = 4000,
  threads = threading(4),
  cores = 4,
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.lmm.vi.prior",
  file_refit = "always"
)

plot(brm.lmm.vi.prior)
conditional_effects(brm.lmm.vi.prior, spaghetti = T, ndraws = 100)
plot(conditional_effects(brm.lmm.vi.prior, re_formula = NULL), points = T)
pp_check(brm.lmm.vi.prior, ndraws = 200)

### Fit model to simulated data ----
brm.lmm.vi <- brm(
  formula = formula_lmm_vi,
  data = df.sim.id,
  family = "gaussian",
  prior = priors.vi.lmm,
  sample_prior = "yes",
  chains = 4,
  # warmup = 8000,
  iter = 4000,
  # control = list(adapt_delta = .95, max_treedepth = 15),
  threads = threading(4),
  cores = 4,
  seed = 42,
  backend = "cmdstan",
  file = "mods/brm.lmm.vi",
  file_refit = "always"
)

plot(brm.lmm.vi)
plot(
  conditional_effects(
    brm.lmm.vi,
    # re_formula = NULL,
    method = "posterior_predict"
  ),
  points = T
)
pp_check(brm.lmm.vi, ndraws = 100)


### Extract slope variance at each dose ----
post_draws.vi <- as_draws_df(brm.lmm.vi)
post_draws.vi <- post_draws.vi %>%
  dplyr::select(
    `sd_ID__PhasePost:Group0.2`,
    `sd_ID__PhasePost:Group0.4`,
    `sd_ID__PhasePost:Group0.6`,
    `sd_ID__PhasePost:Group0.8`,
    `sd_ID__PhasePost:Group1`
  ) %>%
  pivot_longer(
    cols = starts_with("sd"),
    names_to = 'Treatment',
    names_prefix = "sigma_i",
    values_to = "Estimate"
  ) %>%
  mutate(
    Treatment = case_when(
      str_detect(Treatment, "Group0.2") == T ~ .2,
      str_detect(Treatment, "Group0.4") == T ~ .4,
      str_detect(Treatment, "Group0.6") == T ~ .6,
      str_detect(Treatment, "Group0.8") == T ~ .8,
      str_detect(Treatment, "Group1") == T ~ 1
    )
  )
fig.slope_vi <- post_draws.vi %>%
  ggplot(aes(y = Estimate, x = Treatment)) +
  # facet_wrap(~Sex) +
  stat_halfeye() +
  ylim(0, 3) +
  scale_fill_manual(values = colorpal) +
  theme_bw() +
  theme(legend.position = "none")

# Combine into 1 figure ----
fig_metrics <- (fig_predinterval + fig_genotypes + fig_ID) & ylim(0, 1.5)
ggsave(
  filename = "outputs/figs/fig_metrics.jpeg",
  fig_metrics,
  height = 8,
  width = 16,
  units = "cm"
)
