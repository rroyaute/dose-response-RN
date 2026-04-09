library(tidyverse)
library(drc)
library(brms)
library(ggdist)
library(truncnorm)
library(distributional)
library(viridis)

# Figure: Prediction vs. Confidence Interval over EC50 ----
sigma <- .1
c <- 0
d <- 1
b <- 3
e <- .3
x <- seq(0,1, by = .01)

mu <- c + (d - c) / (1 + exp(b * log(x/e)))
mu <- d / (1 + exp(b * log(x/e)))

plot(x, mu, type = "l")

## Simulate data for typical Dose-Response Experiment ----
# 5 doses + 1 control

set.seed(42)
df.sim <- data.frame(x = seq(0,1, length.out = 6)) %>% 
  mutate(x = case_when(x == 0 ~ 0.001, .default = x)) %>% 
  mutate(mu = d / (1 + exp(b * log(x/e)))) %>% 
  mutate(y = rlnorm(n(), log(mu), sigma)) %>% 
  mutate(log_y = log(y))

plot(y~x, df.sim)

## brms estimation ----

# Define the 3-parameter log-logistic as a custom non-linear formula
formula_drm <- bf(
  log(y) ~ log(d / (1 + exp(b * (log(x) - log(e))))),
  b + d + e ~ 1,
  nl = TRUE
)

priors <- c(
  prior(normal(2, 1),     nlpar = "b", lb = 0),   # slope, positive for decreasing curve
  prior(normal(0, 0.5),   nlpar = "d"),             # log(upper asymptote); log(1) = 0 for max y ≈ 1
  prior(normal(0.5, 0.2), nlpar = "e", lb = 0),    # ED50, unchanged
  prior(normal(0, 1),     class = "sigma")
)

inits <- list(
  b_b_Intercept = 2,
  b_d_Intercept = 1,
  b_e_Intercept = 0.5
)
init_list <- rep(list(inits), 4)

mod.brm.prior <- brm(
  formula  = formula_drm,
  data     = df.sim,
  family   = "gaussian",
  prior    = priors,
  sample_prior = "only",
  init     = rep(list(inits), 4),
  chains   = 4, iter = 2000, backend = "cmdstan"
)

plot(mod.brm.prior)
conditional_effects(mod.brm.prior, method = "posterior_predict")
conditional_effects(mod.brm.prior, spaghetti = T, ndraws = 100)

mod.brm <- brm(
  formula  = formula_drm,
  data     = df.sim,
  family   = "gaussian",
  prior    = priors,
  sample_prior = "yes",
  init     = rep(list(inits), 4), 
  seed = 42,
  chains   = 4, 
  iter = 4000, 
  control = list(adapt_delta = .95),
  backend = "cmdstan")

plot(mod.brm)
conditional_effects(mod.brm, method = "posterior_predict")
CRI <- as.data.frame(posterior_summary(mod.brm))

# Prediction intervals via posterior_predict
x_new  <- data.frame(x = seq(0, 1, by = 0.01)) %>% 
  mutate(x = case_when(x == 0 ~ 0.001, .default = x))

pred   <- posterior_predict(mod.brm, newdata = x_new)
epred   <- posterior_epred(mod.brm, newdata = x_new)

pred_df <- data.frame(
  x     = x_new$x,
  y     = apply(pred, 2, function(p) exp(median(p))),
  y_low = apply(pred, 2, function(p) exp(quantile(p, 0.025))),
  y_up  = apply(pred, 2, function(p) exp(quantile(p, 0.975)))
)

epred_df <- data.frame(
  x     = x_new$x,
  y     = apply(epred, 2, function(p) exp(median(p))),
  y_low = apply(epred, 2, function(p) exp(quantile(p, 0.025))),
  y_up  = apply(epred, 2, function(p) exp(quantile(p, 0.975)))
)

fig_predinterval <- ggplot(data = df.sim, aes(x = x, y = y)) +
  # Fitted line
  geom_line(data = pred_df, 
            aes(x = x, y = y), linewidth = 1, color = "dodgerblue") +
  # 95 % Credible Interval for the mean
  geom_ribbon(data = epred_df,
              aes(ymin = y_low, ymax = y_up), 
              alpha = .4, fill = "dodgerblue") +
  # 95 % Prediction Interval
  geom_ribbon(data = pred_df,
              aes(ymin = y_low, ymax = y_up), 
              alpha = .4, fill = "dodgerblue") +
  geom_point(data = CRI,
             aes(x = CRI["b_e_Intercept",1], y = .5), color = "tomato2") +
  geom_errorbarh(aes(y = .5, x = CRI["b_e_Intercept",1],
                     xmin = CRI["b_e_Intercept",3],
                     xmax = CRI["b_e_Intercept",4],
                     width = 0), color = "tomato2") +
  # Data points
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 1) +
  labs(x = "Dose", y = "Phenotype") +
  theme_bw()

ggsave(filename = "outputs/figs/fig_predinterval.jpeg", fig_predinterval)

# Distribution of EC50 values
# Extract posterior samples of e and its variance
post <- as_draws_df(mod.brm)

# Simulate new individual EC50s assuming e ~ lognormal(mu_e, sigma_e)
e_dist <- rlnorm(4000, meanlog = mean(log(post$b_e_Intercept)), 
                sdlog  = sd(log(post$b_e_Intercept)))
hist(e_dist)

# Plot as a distribution
ggplot(data.frame(e = e_dist), aes(x = e)) +
  stat_halfeye(fill = "skyblue") +
  theme_bw()



## drc estimation (TODO) ----
mod.drm <- drm(y~x, data = df.sim, fct = LN.3()) # log-normal model with c = 0 and d = 1
mod.drm <- drm(log_y~x, data = df.sim, fct = LL.3()) # Log-logistic model on log(y)

summary(mod.drm); plot(mod.drm)
coefs <- as.numeric(coef(mod.drm))
CI <- confint(mod.drm, "e")
x <- data.frame(x = seq(0,1, by = .01))
pred <- predict(mod.drm, newdata = x, interval = "prediction") 
pred <- data.frame(x = x, y = pred[,1], 
                   y_low = pred[,2], y_up = pred[,3])

ggplot(data = df.sim, aes(x = x, y = y)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 1) +
  geom_line(data = pred, 
            aes(x = x, y = y), linewidth = 1, color = "dodgerblue") +
  geom_ribbon(data = pred,
              aes(ymin = y_low, ymax = y_up), 
              alpha = .4, fill = "dodgerblue") +
  geom_point(aes(x = coefs[3], y = .5), color = "tomato2") +
  geom_errorbarh(aes(y = .5, x = coefs[3], 
                     xmin = CI[,1], xmax = CI[,2],
                     width = 0), color = "tomato2") +
  theme_bw()

## Old code ----
x <- seq(0,1, by = .15)

mu <- c + (d - c) / (1 + exp(b * log(x/e)))
mu <- d / (1 + exp(b * log(x/e)))
plot(x, mu)

set.seed(42)
df <- data.frame(x = seq(0,1, by = .01)) %>% 
  mutate(mu = d / (1 + exp(b * log(x/e)))) %>% 
  mutate(y = rlnorm(n(), log(mu), sigma)) %>% 
  mutate(y_low = qlnorm(p = .025, log(mu), sigma)) %>% 
  mutate(y_up = qlnorm(p = .975, log(mu), sigma))

mod.drm <- drm(y~x, data = df, fct = LL.3())
summary(mod.drm)
coefs <- as.numeric(coef(mod.drm))
CI <- confint(mod.drm, "e")


df %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = mu), linewidth = 1, color = "dodgerblue") +
  geom_point() +
  geom_ribbon(aes(ymin = y_low, ymax = y_up), alpha = .4, fill = "dodgerblue") +
  geom_point(aes(x = coefs[3], y = .5), color = "tomato2") +
  geom_errorbarh(aes(y = .5, x = coefs[3], 
                     xmin = CI[,1], xmax = CI[,2],
                     width = 0), color = "tomato2") +
  theme_bw()

# Figure: Genotype sensitivity ----
Dose <- seq(0,1, length.out = 6)
n_g <- 5 # 5 genotypes
CVa <- .1 # 10 % of variation around mean for all parameters
sigma_d <- d * CVa # Upper bound variation
sigma_b <- b * CVa # Rate variation
sigma_e <- e * CVa # EC50 sensitivity variation
rho <- -.4 # negative covariance between upper bound and EC50

set.seed(42)
d_g <- rtruncnorm(n_g, mean = d, sd = sigma_d, a = 0)
b_g <- rtruncnorm(n_g,  mean = b, sd = sigma_b, a = 0)
e_g <- rtruncnorm(n_g,  mean = e, sd = sigma_e, a = 0)

Mu <- c(d, b, e)
sigmas <- c(sigma_d, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, rho, 0,
                   rho, 1, 0,
                   0, 0, 1), 
                 nrow = 3)
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
G <- MASS::mvrnorm(n_g, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("d_g", "b_g", "e_g") %>% 
  mutate(G = 1:n_g) %>% 
  arrange(e_g) %>%
  mutate(color_index = row_number())

ggplot(G, aes(x = e_g, y = d_g/2, color = factor(color_index))) + 
  geom_point() + 
  scale_color_viridis_d(option = "H", direction = -1) +
  xlim(0, 1) + ylim(0,1) +
  theme_bw()

# df.sim.points <- G %>%
#   expand(nesting(G, color_index, d_g, b_g, e_g), 
#          Dose = Dose) %>%
#   mutate(mu = d_g / (1 + exp(b_g * log(Dose/e_g)))) %>% 
#   mutate(y = rlnorm(n(), log(mu), sigma)) %>% 
#   mutate(log_y = log(y))

df.sim.lines <- G %>%
  expand(nesting(G, color_index, d_g, b_g, e_g),
         Dose = seq(0, 1, by = 0.01)) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose/e_g))))

fig_genotypes <- ggplot(df.sim.lines, 
       aes(y = y, x = Dose, 
           color = factor(color_index),
           fill = factor(color_index))) +
  geom_line(data = df.sim.lines, 
            aes(y = mu, x = Dose), linewidth = 1) +
  geom_point(data = G,
             aes(y = d_g/2, x = e_g), 
             size = 2.5, shape = 21, fill = "white", alpha = .8) +
  scale_color_viridis_d(option = "H", direction = -1, 
                        begin = .2, end = .8) +
  # stat_halfeye(aes(y = .5, xdist = dist_normal(e, sigma_e)),
  #              color = "black", fill = "grey", alpha = .6) +
  labs(x = "Dose", y = "Phenotype") +
  theme_bw(14) +
  theme(legend.position = "none")

ggsave(filename = "outputs/figs/fig_genotypes.jpeg", fig_genotypes)

# Figure: Pre-Post exposure reaction norms ----
n_doses <- 5
n_id <- 20 # 10 individuals per doses
#log(y) ~ log(d / (1 + exp(b * (log(x) - log(e))))),
CVi <- .1 # 10 % of variation around mean for all parameters
sigma_d <- d * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- -.4 # negative covariance between upper bound and EC50
