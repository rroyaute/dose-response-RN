source("R/funs/dose_response_functions.R")


# Dose-response form and parameter values ----
sigma <- .05
ymin <- 0 # bottom of the dose-response curve
ymax <- 100 # top of the dose-response curve
b <- 8 # slope of the dose-response curve
e <- 50 # EC50 of the dose-response curve
x <- seq(0, 100, by = .01)

mu <- ymin + (ymax - ymin) / (1 + exp(b * log(x / e)))
mu <- ymax / (1 + exp(b * log(x / e)))

plot(x, mu, type = "l")

# Variance change with dose ----
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

# Visualize the among genotype correlations
# ID %>%
#   dplyr::select(d_g:e_g) %>%
#   GGally::ggpairs() +
#   theme_bw()

# Apply dose-response formula per individual over the dose gradient
df.sim.vg <- G %>%
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
df.sim.vg <- df.sim.vg %>%
  mutate(y_sc = as.numeric(scale(y)))

df.sim.vg.long <- df.sim.vg %>% # long format for comparing model to raw values
  dplyr::select(G, Group_n, Group_f, Phase, y, y_sc) %>%
  pivot_wider(names_from = Phase, values_from = c(y, y_sc))

df.test <- df.sim.vg %>%
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

df.test.long <- df.test %>%
  pivot_longer(cols = c(yhat:lnCVa), names_to = "type", values_to = "values")

df.test %>%
  ggplot(aes(x = Dose, y = delta_mu)) +
  geom_line()

df.test %>%
  ggplot(aes(x = Dose, y = delta_sda)) +
  geom_line()

df.test.long %>%
  filter(type %in% c("delta_mu", "delta_sda")) %>%
  ggplot(aes(x = Dose, y = values, group = type, color = type)) +
  geom_line() +
  ylab("Percent difference")

df.test.long %>%
  filter(type %in% c("lnRR", "lnCVa")) %>%
  ggplot(aes(x = Dose, y = values, group = type, color = type)) +
  geom_line() +
  ylab("LnCV")
