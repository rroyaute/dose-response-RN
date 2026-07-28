library(tidyverse)
library(ggdist)
library(ltc)
library(patchwork)
library(distributional)
library(truncnorm)
library(ggthemes)
library(geomtextpath)


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
sigma <- .15
ymin <- 0 # bottom of the dose-response curve
ymax <- 100 # top of the dose-response curve
b <- 8 # slope of the dose-response curve
e <- 50 # EC50 of the dose-response curve
CVi <- .1 # 10 % of variation around mean for all parameters
sigma_ymax <- ymax * CVi # Upper bound variation
sigma_b <- b * CVi # Rate variation
sigma_e <- e * CVi # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters

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
  # ylim(0, 1.5) +
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
    Dose = seq(0, 100, by = 0.01)
  ) %>%
  mutate(mu = d_g / (1 + exp(b_g * log(Dose / e_g))))

pal <- ltc(heatmap3, n = 10, type = "continuous")

fig_drc_vg.B.up <- ggplot(
  df.sim.lines.vg,
  aes(y = y, x = Dose)
) +
  geom_line(
    data = df.sim.lines.vp,
    aes(y = mu, x = Dose),
    linewidth = 2,
    color = "dodgerblue"
  ) +
  geom_line(
    data = df.sim.lines.vg,
    aes(y = mu, x = Dose, color = factor(color_index)),
    linewidth = .8,
    alpha = .8
  ) +
  geom_point(
    data = df.sim.vg,
    aes(y = d_g / 2, x = e_g, color = factor(color_index)),
    size = 2.5,
    shape = 21,
    fill = "white",
    alpha = .8
  ) +
  # scale_color_viridis_d(option = "H", direction = -1, begin = .2, end = .8) +
  scale_color_manual(values = pal) +
  # ylim(0, 1.5) +
  labs(x = "Dose", y = "Phenotype", title = "B) Among-genotypes") +
  theme_custom() +
  theme(legend.position = "none")
fig_drc_vg.B.up

## Figure C: Among-individuals ----
df.sim.lines.vi <- data.frame(Dose = seq(0, 100, by = 0.01)) %>%
  mutate(y = ymax / (1 + exp(b * log(Dose / e))))

fig_drc_vi.C.up <- ggplot(df.sim.vi, aes(y = y, x = Dose)) +
  geom_line(linewidth = .5, aes(group = ID), alpha = .15) +
  geom_point(size = 2.5, shape = 21, fill = "white", alpha = .8) +
  geom_line(
    data = df.sim.lines.vi,
    aes(y = y, x = Dose),
    linewidth = 1,
    color = "dodgerblue"
  ) +
  # ylim(0, 1.5) +
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
ymax_s <- rlnorm(n_draws, meanlog = log(ymax), sdlog = CV_e)
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
fig_drc_vp.A.low

## Figure B: Among-genotypes (genetic variances and correlations) ----
### Genetic variation and Credible Intervals (simulated) ----
CV_e <- .1 # Define measurement error as coefficient of variation
df.sigma.cor <- data.frame(
  values = c(100.00, 0.64, 25.00, 0.4, 0.4, 0.4),
  type = c(rep("Standard deviations", 3), rep("Correlations", 3)),
  name = c(
    "sigma_ymax",
    "sigma_b",
    "sigma_e",
    "ymax x b",
    "ymax x e",
    "b x e"
  )
) %>%
  mutate(sigma_e = abs(CV_e * values))

fig.sigma_g <- df.sigma.cor %>%
  filter(name == "sigma_e") %>%
  ggplot(aes(y = name, xdist = dist_normal(values, sigma_e))) +
  stat_halfeye() +
  labs(x = bquote(sigma[EC[50]]), y = "", title = "Genetic variation") +
  theme_custom() +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

#### Change in Additive genetic variance (CVa) with dose ----
# Sample large amount of genotypes to recreate genetic variance
sigma <- .05
ymin <- 0 # bottom of the dose-response curve
ymax <- 100 # top of the dose-response curve
b <- 8 # slope of the dose-response curve
e <- 50 # EC50 of the dose-response curve
x <- seq(0, 100, by = .01)

Dose <- seq(0, 100, length.out = 100)
n_g <- 100000 # 10 genotypes
CVa <- .1 # 10 % of variation around the mean for all parameters
sigma_ymax <- ymax * CVa # Upper bound variation
sigma_b <- b * CVa # Rate variation
sigma_e <- e * CVa # EC50 sensitivity variation
rho <- .4 # moderate correlation among parameters

set.seed(42)
ymax_g <- rtruncnorm(n_g, mean = ymax, sd = sigma_ymax, a = 0)
b_g <- rtruncnorm(n_g, mean = b, sd = sigma_b, a = 0)
e_g <- rtruncnorm(n_g, mean = e, sd = sigma_e, a = 0)

# Store parameter means as vectors and covariances as matrix
Mu <- c(ymax, b, e)
sigmas <- c(sigma_ymax, sigma_b, sigma_e) # 10 % CV around the mean
rho_mat <- matrix(c(1, rho, rho, rho, 1, rho, rho, rho, 1), nrow = 3)
names <- c("sigma_ymax", "sigma_b", "sigma_e")
colnames(rho_mat) <- names
rownames(rho_mat) <- names
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

# Store individual-specific parameters in a dataframe
set.seed(123456)
G <- MASS::mvrnorm(n_g * length(Dose[2:100]), Mu, Sigma) %>%
  data.frame() %>%
  set_names("d_g", "b_g", "e_g") %>%
  mutate(G = 1:(n_g * length(Dose[2:100]))) %>%
  mutate(assigned_dose = rep(Dose[2:100], each = n_g)) %>%
  mutate(control_dose = .001) %>%
  mutate(Group_n = assigned_dose) %>% # Dose group as numeric
  mutate(Group_f = as.factor(assigned_dose)) %>% # Dose group as factor
  arrange(e_g) %>%
  mutate(color_index = row_number())

# Apply dose-response formula per individual over the dose gradient
df.sim.vg.2 <- G %>%
  pivot_longer(cols = c(assigned_dose:control_dose), values_to = "Dose") %>%
  mutate(mu = ymax_g / (1 + exp(b_g * log(Dose / e_g)))) %>%
  mutate(y = rlnorm(n(), log(mu), sigma)) %>%
  mutate(yhat = ymax / (1 + exp(b * log(Dose / e)))) %>%
  mutate(
    Phase = as.factor(case_when(Dose == 0.001 ~ "Pre", .default = "Post"))
  ) %>%
  mutate(Phase = fct_relevel(Phase, "Pre", "Post"))

# Mean-center data
# Y-values standardized to mean = 0 and SD = 1
# Within individual centring around pre-post exposure doses (taking mid-point between pre-post exposure)
df.sim.vg.2 <- df.sim.vg.2 %>%
  mutate(y_sc = as.numeric(scale(y)))

df.sim.vg.2.sum <- df.sim.vg.2 %>% # Summary statistics
  select(Dose, yhat, mu) %>%
  summarise(
    yhat = mean(yhat),
    sda = sd(mu),
    delta_mu = (yhat - 100), # % change in mean
    delta_sda = (sda - ymax * CVa) / (ymax * CVa) * 100, # % change in genetic variance
    lnRR = log(yhat) - log(100),
    lnCVa = log(sda) - log(ymax * CVa),
    .by = Dose
  )

df.sim.vg.2.sum.long <- df.sim.vg.2.sum %>%
  pivot_longer(cols = c(yhat:lnCVa), names_to = "type", values_to = "values")

# df.sim.vg.2.sum %>%
#   ggplot(aes(x = Dose, y = delta_mu)) +
#   geom_line()
#
# df.sim.vg.2.sum %>%
#   ggplot(aes(x = Dose, y = delta_sda)) +
#   geom_line()

fig.lncva <- df.sim.vg.2.sum.long %>%
  filter(type %in% c("lnRR", "lnCVa")) %>%
  ggplot(aes(x = Dose, y = values, group = type, color = type)) +
  geom_line(linewidth = 1) +
  scale_color_wsj() +
  geom_textline(aes(label = type), vjust = -0.1) +
  ylab("Log Ratio of CV (LnCV)") +
  theme_custom() +
  theme(
    legend.position = "none"
  )
fig.lncva

#### Plot! -----
fig_drc_vg.B.low <- fig.sigma_g / fig.lncva
fig_drc_vg.B.low

## Figure C: Among-individuals (variance partitionning by dose) ----
### Simulation parameters ----
VR_sim <- .2 * ymax # residual variance = 20% of ymax
# CVR <- .1 # 10% residual variance relative to mean at a given dose

### Variance decomposition ----
# Calculate variance explained by:
# population slope at each dose (VF)
# Individual slope (VS)
# Individual intercept (VI at the mid-dose point)
# within-individual (residual) variance (VR)

# Calculate pre-post changes (slope) at each dose
# Fixing x = 0 at dose for post-exposure phase
df.pre.post <- df.sim.vi %>%
  select(ID, Group_n, Phase, Dose, mu) %>%
  pivot_wider(names_from = Phase, values_from = c(Dose, mu)) %>%
  mutate(
    D_g = Dose_Post - Dose_Pre,
    sd_x = D_g / 2, # SD for two points {0, D_g}
    mu_Post = mu_Post,
    mu_Pre = mu_Pre,
    u0_i = (mu_Post + mu_Pre) / 2, # intercept at the mid-dose point
    u1_i = (mu_Post - mu_Pre) / D_g # individual slope
  )

# Variance breakdown in exposed groups
var_exposed <- df.pre.post %>%
  group_by(Group_n) %>%
  summarise(
    D_g = first(D_g),
    VI = var(u0_i),
    VS = var(u1_i * sd_x),
    mu = mean(mu_Post),
    beta_g = (mean(mu_Post) - mean(mu_Pre)) / first(D_g),
    beta_std = beta_g * sd_x, # Standardized slope
    VF = beta_std^2,
    VR = VR_sim,
    .groups = "drop"
  )

df.ctrl <- df.pre.post %>% filter(Group_n == first(Group_n)) # All individuals from pre-exposure phase
u0_i_ctrl <- df.sim.vi %>%
  filter(Phase == "Pre") %>%
  summarise(u0_i_ctrl = mu) # Intercepts at dose = 0 in pre-exposure phase
u0_i_ctrl <- as.numeric(u0_i_ctrl$u0_i_ctrl)

var_control <- tibble(
  Group_n = 0,
  D_g = 0,
  VI = var(u0_i_ctrl),
  VS = 0, # No slope variation in controls
  mu = mean(df.ctrl$mu_Pre),
  beta_std = 0, # No population slope variation in controls
  VF = 0, # No population slope variation in controls
  VR = VR_sim
  # VR = CVR * 1 # 10% residual variation around average
)

df.var.comp <- bind_rows(var_control, var_exposed) %>%
  arrange(Group_n) %>%
  mutate(
    mu = mu,
    VI = VI, # among-individual variance
    VS = VS, # slope variance
    VF = VF, # Fixed effect variance (attributed to average changes with dose),
    VR = VR, # Residual variance
    VP = VI + VF + VR, # total phenotypic variance
    R2_I = VI / VP * 100, # % variation explained by among-individual variance
    R2_S = VS / VP * 100, # % variation explained by slope variation
    R2_Res = VR / VP * 100, # % Residual variation
    lnRR = log(mu / 100),
    lnCVI = log(sqrt(VI) / mu) - log(sd(u0_i_ctrl) / mu),
    lnCVR = log(sqrt(VR) / mu) - log(sqrt(VR_sim) / mu),
  )
df.var.comp

### Plot! ----
pal <- ltc(expevo)
fig_varcomp.stack <- df.var.comp %>%
  select(Group_n, VI, VS, VF, VR) %>%
  pivot_longer(!Group_n, names_to = "type", values_to = "values") %>%
  mutate(type = fct_relevel(type, "VR", "VF", "VS", "VI")) %>%
  ggplot(aes(x = Group_n, y = values, fill = type)) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  scale_fill_manual(values = pal) +
  labs(
    x = "Dose",
    y = "Proportion of variance explained",
    title = "Variance partitioning by dose "
  ) +
  annotate(
    "text",
    color = "white",
    x = c(0.1, 20, 40, .1),
    y = c(.1, .75, .85, .95),
    label = c("VI", "VS", "VF", "VR")
  ) +
  theme_custom()
fig_varcomp.stack

fig_varcomp.lnCV <- df.var.comp %>%
  select(Group_n, lnCVI, lnCVR, lnRR) %>%
  pivot_longer(!Group_n, names_to = "type", values_to = "values") %>%
  mutate(type = fct_relevel(type, "lnRR", "lnCVR", "lnCVI")) %>%
  ggplot(aes(
    x = Group_n,
    y = values,
    fill = type,
    group = type,
    color = type,
    label = type
  )) +
  geom_point() +
  geom_line() +
  # scale_color_wsj() +
  scale_color_manual(values = pal) +
  geom_textline(aes(label = type), vjust = 0.1, hjust = .7) +
  labs(
    x = "Dose",
    y = "lnCV",
    title = "Magnitude of variance change with dose"
  ) +
  theme_custom()
# theme(legend.position = "none")
fig_varcomp.lnCV

fig_drc_vi.C.low <- fig_varcomp.stack / fig_varcomp.lnCV
fig_drc_vi.C.low

# Combine into 1 figure ----
fig_metrics.up <- (fig_drc_vp.A.up + fig_drc_vg.B.up + fig_drc_vi.C.up)
fig_metrics.up <- fig_metrics.up +
  plot_annotation(
    title = "Levels of Biological organisation",
    theme = theme(plot.title = element_text(size = 20))
  )
fig_metrics.up

fig_metrics.low <- (fig_drc_vp.A.low + fig_drc_vg.B.low + fig_drc_vi.C.low)
fig_metrics.low <- fig_metrics.low +
  plot_annotation(
    title = "What to report?",
    theme = theme(plot.title = element_text(size = 20))
  )
fig_metrics.low

fig_metrics <- fig_metrics.up / fig_metrics.low
fig_metrics

ggsave(
  filename = "outputs/figs/fig_metrics.jpeg",
  fig_metrics,
  height = 37,
  width = 37,
  units = "cm"
)
