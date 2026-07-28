library(tidyverse)
library(ggdist)
library(truncnorm)
library(distributional)

# Code for simulating datasets for phenotypic level, among-genotype level and
# individual-level differences in dose responses
# This script is only used to generate the datasets found in the data subfolder

# Import convenience functions
source("R/funs/dose_response_functions.R")

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
## Simulate data for typical Dose-Response Experiment -----
# 5 doses + 1 control

set.seed(42)
df.sim.vp <- data.frame(Dose = seq(0, 100, length.out = 6)) %>% # Sample 6 values between [0;1] with equal spacing
  mutate(Dose = case_when(Dose == 0 ~ 0.001, .default = Dose)) %>% # Replace dose = 0 with small value to avoid computational issues
  mutate(mu = ymax / (1 + exp(b * log(Dose / e)))) %>% # Apply dose-response equation to all x-values
  mutate(y = rlnorm(n(), log(mu), sigma)) %>% # Sample from log-normal distribution to keep y-values > 0
  mutate(log_y = log(y))

plot(y ~ Dose, df.sim.vp)
write.csv(df.sim.vp, "data/df.sim.vp.csv")

## Monte Carlo simulations to recover confidence and prediction intervals ----
cv <- 0.05 # 5 % coefficient of variation on each parameter
n_mc <- 20000
set.seed(42)

dose_seq <- seq(0.001, 100, length.out = 300)

# Sample parameters (lognormal to keep values positive)
ymax_s <- rlnorm(n_mc, meanlog = log(ymax), sdlog = cv)
b_s <- rlnorm(n_mc, meanlog = log(b), sdlog = cv)
e_s <- rlnorm(n_mc, meanlog = log(e), sdlog = cv)

# Evaluate μ for every MC draw × every dose
# Result: n_mc × length(dose_seq) matrix
mu_mat <- vapply(
  seq_len(n_mc),
  function(i) drc_mu(dose_seq, ymax_s[i], b_s[i], e_s[i]),
  numeric(length(dose_seq))
) |>
  t() # → [n_mc × n_dose]

# Confidence interval = quantiles of μ across MC draws
ci_low <- apply(mu_mat, 2, quantile, 0.025)
ci_up <- apply(mu_mat, 2, quantile, 0.975)

# Prediction interval = CI propagated through the lognormal noise ────────
# For each MC draw, simulate one observation per dose, then take quantiles.
# rlnorm(n, log(mu), sigma) on the whole matrix at once:
y_pred_mat <- matrix(
  rlnorm(n_mc * length(dose_seq), meanlog = log(mu_mat), sdlog = sigma),
  nrow = n_mc
)
pi_low <- apply(y_pred_mat, 2, quantile, 0.025)
pi_up <- apply(y_pred_mat, 2, quantile, 0.975)

# Assemble into dataframe
df.sim.lines.vp <- tibble(
  Dose = dose_seq,
  mu = drc_mu(dose_seq, ymax, b, e), # nominal mean line
  ci_low = ci_low,
  ci_up = ci_up,
  pi_low = pi_low,
  pi_up = pi_up
)
write.csv(df.sim.lines.vp, "data/df.sim.lines.vp.csv")

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
    Dose = seq(0, 1, by = 0.01)
  ) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose / e_g)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma))

write.csv(df.sim.vg, "data/df.sim.vg.csv")


# Among-individual level ----
Dose <- seq(0, 100, length.out = 6)
n_id <- 20 # 20 individuals per doses
CVi <- .1 # 10 % of variation around mean for all parameters
sigma_ymax <- ymax * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters

set.seed(42)
ymax_i <- rtruncnorm(n_id, mean = ymax, sd = sigma_ymax, a = 0)
b_i <- rtruncnorm(n_id, mean = b, sd = sigma_b, a = 0)
e_i <- rtruncnorm(n_id, mean = e, sd = sigma_e, a = 0)

# Store parameter means as vectors and covariances as matrix
Mu <- c(ymax, b, e)
sigmas <- c(sigma_ymax, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = 3)
names <- c("sigma_ymax", "sigma_b", "sigma_e")
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

# Store individual-specific parameters in a dataframe
set.seed(42)
ID <- MASS::mvrnorm(n_id * length(Dose[2:6]), Mu, Sigma) %>%
  data.frame() %>%
  set_names("d_i", "b_i", "e_i") %>%
  mutate(ID = 1:(n_id * length(Dose[2:6]))) %>%
  mutate(assigned_dose = rep(Dose[2:6], each = n_id)) %>%
  mutate(control_dose = .001) %>%
  mutate(Group_n = assigned_dose) %>% # Dose group as numeric
  mutate(Group_f = as.factor(assigned_dose)) %>% # Dose group as factor
  arrange(e_i) %>%
  mutate(color_index = row_number())

# Visualize the among genotype correlations
ID %>%
  dplyr::select(d_i:e_i) %>%
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

write.csv(df.sim.vi, "data/df.sim.vi.csv")
write.csv(df.sim.vi.long, "data/df.sim.vi.long.csv")
