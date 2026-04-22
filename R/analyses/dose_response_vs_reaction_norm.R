library(tidyverse)
library(patchwork)

# Set global plotting theme ----
theme_custom <- function() {
  theme_classic(16) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}


# Dose-response figure ----
sigma <- .1
c <- 0 # bottom of the dose-response curve
d <- 1 # top of the dose-response curve
b <- 3 # slope of the dose-response curve
e <- .3 # EC50 of the dose-response curve
x <- seq(0, 1, by = .01)

mu <- c + (d - c) / (1 + exp(b * log(x / e)))
mu <- d / (1 + exp(b * log(x / e)))

plot(x, mu, type = "l")

df.drc <- data.frame(x = x, y = mu)
fig.drc <- df.drc %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(color = "dodgerblue", linewidth = 1.5) +
  geom_vline(
    xintercept = e,
    color = "darkred",
    linewidth = .5,
    linetype = "dashed"
  ) +
  annotate("text", x = e, y = .12, label = "EC50") +
  labs(x = "Dose", y = "Phenotype", title = "Dose Response Approach") +
  theme_custom()
fig.drc
ggsave(fig.drc, "outputs/figs/fig.drc.jpeg")


# Reaction norm figure ----
n_id <- 4
x <- seq(-2, 2, by = .01) # Environmental gradient
sigma_w <- .3 # Within-individual (residual) variance
b0 <- 0 # Average phenotype at the center of the environmental gardient (x = .5)
b1 <- -.5 # Average slope

ID <- data.frame(
  ID = 1:n_id,
  b0_i = c(-2, -1, 1, 2),
  b1_i = c(-.5, -.2, .2, .5)
)

df.rn <- ID %>%
  expand(nesting(ID, b0_i, b1_i), x = x) %>%
  mutate(yhat = (b0 + b0_i) + (b1 + b1_i) * x) %>%
  mutate(y = rnorm(n(), yhat, sigma_w))
df.rn.subset <- df.rn %>%
  filter(x %in% c(-1.5, -.5, .5, 1.5))
df.rn.i <- df.rn %>%
  filter(x == .5)

colorpal <- c("#c72e29", "#016392", "#be9c2e", "#098154")
fig.rn <- df.rn %>%
  ggplot(aes(
    x = x,
    y = yhat,
    group = ID,
    color = factor(ID),
    fill = factor(ID)
  )) +
  # Individual reaction norms
  geom_abline(intercept = b0, slope = b1, linewidth = 1.5) +
  # Individual reaction norms
  geom_line(linewidth = 1.5) +
  # Individual points
  geom_point(
    data = ID,
    aes(x = 0, y = b0 + b0_i),
    size = 4,
    shape = 21,
    color = "black",
    fill = "white"
  ) +
  # Individual points
  geom_point(data = df.rn.subset, aes(x = x, y = y), size = 3, alpha = .4) +
  scale_color_manual(values = colorpal) +
  scale_fill_manual(values = colorpal) +
  labs(
    x = "Environmental gradient",
    y = "Phenotype",
    title = "Reaction Norm Approach"
  ) +
  theme_custom()
fig.rn
ggsave(plot = fig.rn, "outputs/figs/fig.rn.jpeg")

# Both figures side-by-side
fig.drc.vs.rn <- fig.drc + fig.rn
ggsave(
  plot = fig.drc.vs.rn,
  "outputs/figs/fig.drc.vs.rn.jpeg",
  height = 4,
  width = 8
)


# Old code ----
sigma_w <- .2 # within-individual (residual) variance
sigma_i <- 1 # among-individual variance
sigma_s <- .3 # slope variance
rho <- .8 # slope-intercept correlation

mu <- c(b0, b1) # vector of means
sigmas <- c(sigma_i, sigma_s) # vector of SDs
rho_mat <- matrix(
  c(
    1,
    rho, # correlation matrix
    rho,
    1
  ),
  nrow = 2
)
Sigma <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # Covariance matrix

# Generate dataframe for 4 individuals with specific intercept and slope values
set.seed(123456)
ID <- MASS::mvrnorm(n_id, mu, Sigma) %>%
  data.frame() %>%
  set_names("b0_i", "b1_i") %>%
  mutate(ID = 1:n_id)

# Simulate individual values over the environmental gradient
