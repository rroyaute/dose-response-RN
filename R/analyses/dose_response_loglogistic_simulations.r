library(tidyverse)
library(ggdist)
library(truncnorm)
library(distributional)

source("R/funs/dose_response_functions.R")

# Code for simulating datasets for phenotypic level, among-genotype level and
# individual-level differences in dose responses
# This script is only used to generate the datasets found in the data subfolder

# Dose-response form and parameter values ----
sigma <- .15
ymin <- 0 # bottom of the dose-response curve
ymax <- 100 # top of the dose-response curve
b <- 8 # slope of the dose-response curve
e <- 50 # EC50 of the dose-response curve
x <- seq(0, 100, by = .01)

mu <- ymin + (ymax - ymin) / (1 + exp(b * log(x / e)))
mu <- ymax / (1 + exp(b * log(x / e)))

plot(x, mu, type = "l")

# Phenotypic level ----
CVe <- .025 # error term for average prediction
e_up <- e + 1.96 * e * CVe
e_low <- e - 1.96 * e * CVe

n_obs <- 1:1000
Dose <- seq(0, 100, by = .01)

# Simulate dose-response curve according to parameters
df.sim.vp <- data.frame(n = n_obs) %>%
  mutate(
    ymax_n = rtruncnorm(n(), mean = ymax, sd = ymax * CVe, a = 0),
    e_n = rtruncnorm(n(), mean = e, sd = e * CVe, a = 0),
    b_n = rtruncnorm(n(), mean = b, sd = b * CVe, a = 0)
  ) %>%
  expand(nesting(n_obs, ymax_n, e_n, b_n), Dose) %>% # Sample 6 values between [0;1] with equal spacing
  mutate(Dose = case_when(Dose == 0 ~ 0.001, .default = Dose)) %>% # Replace dose = 0 with small value to avoid computational issues
  mutate(
    mu = ymax / (1 + exp(b * log(Dose / e))),
    yhat = ymax_n / (1 + exp(b_n * log(Dose / e_n))), # predicted average including measurement error on e
    y = rlnorm(n(), log(mu), sigma) # Apply dose-response equation to all x-values
  )

# Store average and sd values for the dose-response
df.sim.vp.pred <- df.sim.vp %>%
  summarise(
    mu = mean(mu), # Average trend (according to true parameter values)
    mu_yhat = mean(yhat), # Average trend (accounting for sampling error on e)
    sd_yhat = sd(yhat), # variation around average trend (accounting for sampling error on e)
    yhat_low = mu_yhat - 1.96 * sd_yhat, # Lower CI of the mean
    yhat_up = mu_yhat + 1.96 * sd_yhat, # Upper CI of the mean
    # sd_mu = sd(mu_y),
    sd_y = sd(y), # Variation around the average trend accounting for the residual variation term sigma
    y_low = mu - 1.96 * sd_y, # Lower prediction interval (PRI)
    y_up = mu + 1.96 * sd_y, # Upper prediction interval (PRI)
    .by = Dose
  )

# Store prediction interval & confidence interval values
df.intervals.vp <- data.frame(
  low.value = rep(NA, 3),
  up.value = rep(NA, 3),
  conf.type = c("pre", "ci", "e")
)

df.intervals.vp[1, 1] <- df.sim.vp.pred %>%
  filter(abs(y_low - 50) <= .1) %>%
  summarise(value = mean(Dose)) # Lower prediction interval value
df.intervals.vp[1, 2] <- df.sim.vp.pred %>%
  filter(abs(y_up - 50) <= .1) %>%
  summarise(value = mean(Dose)) # Upper prediction interval value
df.intervals.vp[2, 1] <- df.sim.vp.pred %>%
  filter(abs(yhat_low - 50) <= .1) %>%
  summarise(value = mean(Dose)) # Lower confidence interval value
df.intervals.vp[2, 2] <- df.sim.vp.pred %>%
  filter(abs(yhat_up - 50) <= .1) %>%
  summarise(value = mean(Dose)) # Upper confidence interval value
df.intervals.vp[3, 1] <- e_low # Lower confidence interval of the EC50 parameter
df.intervals.vp[3, 2] <- e_up # Lower confidence interval of the EC50 parameter

# Export dataframes
write.csv(df.sim.vp, "data/df.sim.vp.csv")
write.csv(df.sim.vp.pred, "data/df.sim.vp.pred.csv")
write.csv(df.intervals.vp, "data/df.intervals.vp.csv")

# Among-genotype level ----
Dose <- seq(0, 100, length.out = 6)
n_g <- 10 # 10 genotypes
CVa <- .1 # 10 % of variation around the mean for all parameters
sigma_ymax <- ymax * CVa # Upper bound variation
sigma_b <- b * CVa # Rate variation
sigma_e <- e * CVa # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters
# Genotypes with higher basal expression are less sensitive
# r_dxe = rho: higher basal expression <-> higher EC50
# r_dxb = rho : higher basal expression <-> higher slope
# r_bxe = rho : higher EC50 <-> higher slope

set.seed(42)
ymax_g <- rtruncnorm(n_g, mean = ymax, sd = sigma_ymax, a = 0)
b_g <- rtruncnorm(n_g, mean = b, sd = sigma_b, a = 0)
e_g <- rtruncnorm(n_g, mean = e, sd = sigma_e, a = 0)

# Store parameter means as vectors and covariances as matrix
Mu <- c(ymax, b, e)
sigmas <- c(sigma_ymax, sigma_b, sigma_e) # 10 % CV around the mean
names <- c("sigma_ymax", "sigma_b", "sigma_e")
rho_mat <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = 3) # Correlation matrix
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

# Store genotype-specific parameters in a dataframe
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

# Apply dose-response formula per genotype over the dose gradient
df.sim.vg <- G %>%
  expand(
    nesting(G, color_index, d_g, b_g, e_g),
    Dose = seq(0, 1, by = 0.01) # Check!
  ) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose / e_g)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma))

# Export dataframes
write.csv(df.sim.vg, "data/df.sim.vg.csv")


# Among-individual level ----
Dose <- seq(0, 100, length.out = 6)
n_id <- 20 # 20 individuals per doses
CVi <- .1 # 10 % of variation around mean for all parameters
sigma_ymax <- ymax * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- .5 # strong correlation among parameters

set.seed(42)
ymax_i <- rtruncnorm(n_id, mean = ymax, sd = sigma_ymax, a = 0)
b_i <- rtruncnorm(n_id, mean = b, sd = sigma_b, a = 0)
e_i <- rtruncnorm(n_id, mean = e, sd = sigma_e, a = 0)

# Store parameter means as vectors and covariances as matrix
Mu <- c(ymax, b, e)
sigmas <- c(sigma_ymax, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, rho, rho, rho, 1, -rho, rho, -rho, 1), nrow = 3)
names <- c("sigma_ymax", "sigma_b", "sigma_e")
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

# Store individual-specific parameters in a dataframe
set.seed(42)
ID <- MASS::mvrnorm(n_id * length(Dose[2:6]), Mu, Sigma) %>%
  data.frame() %>%
  set_names("ymax_i", "b_i", "e_i") %>%
  mutate(ID = 1:(n_id * length(Dose[2:6]))) %>%
  mutate(assigned_dose = rep(Dose[2:6], each = n_id)) %>%
  mutate(control_dose = .001) %>%
  mutate(Group_n = assigned_dose) %>% # Dose group as numeric
  mutate(Group_f = as.factor(assigned_dose)) %>% # Dose group as factor
  arrange(e_i) %>%
  mutate(color_index = row_number())

# Visualize the among genotype correlations
ID %>%
  dplyr::select(ymax_i:e_i) %>%
  GGally::ggpairs() +
  theme_bw()

# Apply dose-response formula per individual over the dose gradient
df.sim.vi <- ID %>%
  pivot_longer(cols = c(assigned_dose:control_dose), values_to = "Dose") %>%
  mutate(mu = ymax_i / (1 + exp(b_i * log(Dose / e_i)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma)) %>%
  mutate(
    Phase = as.factor(case_when(Dose == 0.001 ~ "Pre", .default = "Post"))
  ) %>%
  mutate(Phase = fct_relevel(Phase, "Pre", "Post"))

# Mean-center data
# Y-values standardized to mean = 0 and SD = 1
# Within individual centring around pre-post exposure doses (taking mid-point between pre-post exposure)
df.sim.vi <- df.sim.vi %>%
  mutate(y_sc = as.numeric(scale(y)))

df.sim.vi.long <- df.sim.vi %>% # long format for comparing model to raw values
  dplyr::select(ID, Group_n, Group_f, Phase, y, y_sc) %>%
  pivot_wider(names_from = Phase, values_from = c(y, y_sc))

# Export dataframes
write.csv(df.sim.vi, "data/df.sim.vi.csv")
write.csv(df.sim.vi.long, "data/df.sim.vi.long.csv")
