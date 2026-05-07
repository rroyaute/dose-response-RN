library(tidyverse)
library(ggdist)
library(viridis)
library(patchwork)
library(distributional)

# Load simulated datasets ----
df.sim.vp <- read.csv("data/df.sim.vp.csv")
df.sim.lines.vp <- read.csv("data/df.sim.lines.vp.csv")
df.sim.vg <- read.csv("data/df.sim.vg.csv")
df.sim.vi <- read.csv("data/df.sim.vi.csv")
df.sim.vi.long <- read.csv("data/df.sim.vi.long.csv")


# Set plotting theme ----
theme_custom <- function() {
  theme_classic(12) +
    theme(legend.position = "none")
}
# Dose-response form and parameter values ----
sigma <- .1
c <- 0 # bottom of the dose-response curve
d <- 1 # top of the dose-response curve
b <- 3 # slope of the dose-response curve
e <- .3 # EC50 of the dose-response curve
# x <- seq(0, 1, by = .01)
# mu <- d / (1 + exp(b * log(x / e)))

# Plot setup
# Upper rows
# title = "Levels of organisation",
# A) subtitle = "Phenotypic level"
# B) subtitle = "Among-genotype level"
# C) subtitle = "Among-individual level"
# Lower rows
# title = "What to report?",
# A) subtitle = "Phenotypic variation in EC50"
# B) subtitle = "Genetic variation in EC50"
# C) subtitle = "Elevation and slope variation with dose"

# Upper pannel: Dose response depending on biological level of organisation ----
## Figure A: Phenotypic-level ----
fig_drc_vp.A.up <- ggplot(data = df.sim.lines.vp, aes(x = Dose, y = mu)) +
  # 95 % Prediction interval (wider, lighter)
  geom_ribbon(
    data = df.sim.lines.vp,
    aes(ymin = pi_low, ymax = pi_up),
    alpha = .20,
    fill = "dodgerblue"
  ) +
  # 95 % Confidence interval on the mean (narrower, darker)
  geom_ribbon(
    data = df.sim.lines.vp,
    aes(ymin = ci_low, ymax = ci_up),
    alpha = .45,
    fill = "dodgerblue"
  ) +
  # Fitted line
  geom_line(
    data = df.sim.lines.vp,
    aes(x = Dose, y = mu),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  # Data points
  geom_point(
    data = df.sim.vp,
    aes(x = Dose, y = y),
    shape = 21,
    colour = "black",
    fill = "white",
    size = 2,
    stroke = 1
  ) +
  ylim(0, 1.5) +
  labs(
    x = "Dose",
    y = "Phenotype",
    title = "A) Dose-response curve\n(Phenotypic level)"
  ) +
  theme_custom()
fig_drc_vp.A.up

## Figure B: Among-genotypes ----
df.sim.lines.vg <- df.sim.vg %>%
  expand(
    nesting(G, color_index, d_g, b_g, e_g),
    Dose = seq(0, 1, by = 0.01)
  ) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose / e_g))))

fig_drc_vg.B.up <- ggplot(
  df.sim.lines.vg,
  aes(y = y, x = Dose, color = factor(color_index), fill = factor(color_index))
) +
  geom_line(data = df.sim.lines.vg, aes(y = mu, x = Dose), linewidth = 1) +
  geom_point(
    data = df.sim.vg,
    aes(y = d_g / 2, x = e_g),
    size = 2.5,
    shape = 21,
    fill = "white",
    alpha = .8
  ) +
  scale_color_viridis_d(option = "H", direction = -1, begin = .2, end = .8) +
  # stat_halfeye(aes(y = .5, xdist = dist_normal(e, sigma_e)),
  #              color = "black", fill = "grey", alpha = .6) +
  ylim(0, 1.5) +
  labs(x = "Dose", y = "Phenotype", title = "B) Among-genotypes") +
  theme_custom() +
  theme(legend.position = "none")
fig_drc_vg.B.up

## Figure C: Among-individuals ----
df.sim.lines.vi <- data.frame(Dose = seq(0, 1, by = 0.01)) %>%
  mutate(y = d / (1 + exp(b * log(Dose / e))))

fig_drc_vi.C.up <- ggplot(df.sim.vi, aes(y = y, x = Dose)) +
  geom_line(linewidth = .5, aes(group = ID), alpha = .15) +
  geom_point(size = 2.5, shape = 21, fill = "white", alpha = .8) +
  geom_line(
    data = df.sim.lines.vi,
    aes(y = y, x = Dose),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  ylim(0, 1.5) +
  labs(x = "Dose", y = "Phenotype", title = "C) Among-individuals") +
  theme_custom() +
  theme(legend.position = "none")
# fig_drc_vi.C.up <- fig_drc_vi.C.up + plot_annotation()
fig_drc_vi.C.up

# Lower pannel: Parameters to report ----
## Figure A: Phenotypic-level (CI and PRI) ----
CV_e <- .05 # Define measurement error as coefficient of variation
n_draws <- 20000 # Number of random draws
# Random draws for each dose-response parameter
d_s <- rlnorm(n_draws, meanlog = log(d), sdlog = CV_e)
b_s <- rlnorm(n_draws, meanlog = log(b), sdlog = CV_e)
e_s <- rlnorm(n_draws, meanlog = log(e), sdlog = CV_e)

# Find doses for which y = d/2 among prediction draws
# Dose_EC50 = e * exp(log(2*exp(epsilon) - 1)/b)
# Where epsilon is random noise parameter taken from rnorm(0, sigma)
epsilon_s <- rnorm(n_draws, 0, sigma)
prod <- 2 * exp(epsilon_s) - 1 # Take only positive values for this product
valid <- prod > 0 # Take only positive values
ec50_pi_vec <- e_s[valid] * exp(log(prod[valid]) / b_s[valid])

df.EC50dist.vp <- data.frame(
  param = c(
    rep("Prediction interval", n_draws),
    rep("Confidence interval", n_draws)
  ),
  EC50_dist = c(ec50_pi_vec, e_s)
)

fig_drc_vp.A.low <- df.EC50dist.vp %>%
  ggplot(aes(x = EC50_dist, y = param)) +
  stat_halfeye() +
  labs(x = "Dose", y = "EC50", title = "Phenotypic variation in EC50") +
  theme_custom() +
  theme(legend.position = "none")

## Figure B: Among-genotypes (genetic variances and correlations) ----
CV_e <- .2 # Define measurement error as coefficient of variation
df.sigma.cor <- data.frame(
  values = c(0.10, 0.30, 0.03, -0.4, 0.4, -0.4),
  type = c(rep("Standard deviations", 3), rep("Correlations", 3)),
  name = c(
    "sigma_Basal",
    "sigma_Rate",
    "sigma_EC50",
    "Basal x Rate",
    "Basal x EC50",
    "Rate x EC50"
  )
) %>%
  mutate(sigma_e = abs(CV_e * values))

df.sigma.cor %>%
  filter(name == "sigma_EC50") %>%
  ggplot(aes(y = name, xdist = dist_normal(values, sigma_e))) +
  stat_halfeye() +
  labs(x = bquote(sigma[EC[50]]), y = "", title = "Genetic variation") +
  theme_custom() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

df.sigma.cor %>%
  filter(type == "Correlations") %>%
  ggplot(aes(y = name, xdist = dist_normal(values, sigma_e))) +
  stat_halfeye() +
  labs(x = bquote(sigma[EC[50]]), y = "", title = "Genetic correlations") +
  theme_custom() +
  theme(legend.position = "none")

## Figure C: Among-individuals (variance partitionning by dose) ----
### Simulation parameters ----
CVi <- .1 # 10 % of variation around mean for all parameters
VR <- .1 # 10% variation at Dose = 0
sigma_d <- d * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters

# Store parameter means as vectors and covariances as matrix
Mu <- c(d, b, e)
sigmas <- c(sigma_d, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, -rho, rho, -rho, 1, -rho, rho, -rho, 1), nrow = 3)
names <- c("sigma_d", "sigma_b", "sigma_e")
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

### Variance decomposition ----
# Calculate variance explained by:
# population slope at each dose (VD)
# Individual slope (Vv x Vx + mu^2 * Vv)
# Individual intercept (Vu + 2 * Cov(u,v))
# within-individual (residual) variance (VR)

# Population slope variance = beta^2 * VD
df.VD <- df.sim.vi %>%
  summarise(
    y_pre = y[Phase == "Pre"],
    y_post = y[Phase == "Post"],
    y_diff = y[Phase == "Post"] - y[Phase == "Pre"],
    .by = c(Group_n, ID)
  ) %>%
  summarise(y_slope = mean(y_diff) / mean(Group_n), .by = Group_n) %>%
  mutate(VD = y_slope^2 * Group_n^2 / 12) %>%
  arrange(Group_n)

# Means for each dose category
df.mu <- df.sim.vi %>%
  summarise(mu = mean(y), .by = Group_n) %>%
  arrange(Group_n)

# Individual slopes variance (VS)
df.VS <- df.sim.vi %>%
  arrange(Dose, Group_n) %>%
  summarise(
    y_pre = y[Phase == "Pre"],
    y_post = y[Phase == "Post"],
    y_diff = (y[Phase == "Post"] - y[Phase == "Pre"]),
    .by = c(ID, Group_n)
  ) %>%
  mutate(y_slope = y_diff / Group_n, V_x = Group_n^2 / 12) %>%
  summarise(
    V_v = var(y_slope),
    Vu = sigma_d, # Variance at x = 0: variance in d parameter from simulations
    .by = c(Group_n, V_x)
  ) %>%
  mutate(Vvx = V_v * V_x) %>% # Variance in slope (not accounting for mean changes)
  left_join(df.mu) %>% # Add average value per group
  mutate(VS = Vvx + mu^2 * V_v)

# Individual intercepts (Vu + 2*mu*C_uv)
# Calculate intercept-slope covariance
df.C_uv <- df.sim.vi %>%
  arrange(Dose, Group_n) %>%
  summarise(
    y_pre = y[Phase == "Pre"],
    y_post = y[Phase == "Post"],
    y_diff = (y[Phase == "Post"] - y[Phase == "Pre"]),
    .by = c(ID, Group_n)
  ) %>%
  mutate(y_slope = y_diff / Group_n, .by = c(ID, Group_n)) %>%
  summarise(C_uv = cov(y_pre, y_slope), .by = c(Group_n))

# Merge into one dataframe
df.var.comp <- left_join(df.VS, df.VD) %>%
  left_join(df.C_uv) %>%
  add_row(
    Group_n = 0,
    Vu = .1,
    mu = 1,
    VS = 0,
    C_uv = 0,
    VD = 0,
    .before = 1
  ) %>% # Add partitioning at Dose = 0
  mutate(
    VI = VS + Vu + 2 * mu * C_uv, # among individual variance over dose gradient
    VR = VR
  ) %>% # Residual (within-individual variance)
  mutate(VP = VD + VI + VR) %>% # Total phenotypic variance
  mutate(
    Rmar = VI / VP, # marginalized repeatability
    R2S = VS / VP
  ) # Variance explained by among-individual slope differences
df.var.comp

df.var.comp %>%
  select(Group_n, VD, VI, VR) %>%
  pivot_longer(!Group_n, names_to = "type", values_to = "values") %>%
  ggplot(aes(x = Group_n, y = values, fill = type)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Dose", y = "Variance", title = "Variance partitioning by dose ") +
  theme_custom() +
  theme(legend.position = "none")


# Combine into 1 figure ----
fig_metrics.up <- (fig_drc_vp.A.up + fig_drc_vg.B.up + fig_drc_vi.C.up)
fig_metrics.up <- fig_metrics.up +
  plot_annotation(
    title = "Levels of Biological organisation",
    theme = theme(plot.title = element_text(size = 20))
  )
fig_metrics.up

ggsave(
  filename = "outputs/figs/fig_metrics.jpeg",
  fig_metrics.up,
  height = 8,
  width = 16,
  units = "cm"
)
