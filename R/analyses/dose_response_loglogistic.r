library(tidyverse)
library(drc)
library(brms)
library(ggdist)

# Figure 1: Prediction vs. Confidence Interval over EC50 ----
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

ggplot(data = df.sim, aes(x = x, y = y)) +
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
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) +
  theme_bw()

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

# Figure 2: Pre-Post exposure reaction norms ----


