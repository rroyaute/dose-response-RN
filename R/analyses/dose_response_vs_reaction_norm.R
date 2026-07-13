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
  # Add ymax
  annotate(
    geom = "text",
    x = 0,
    y = .9,
    fontface = "bold",
    label = expression(y[max]),
    size = 5,
    color = "dodgerblue"
  ) +
  # Add ymin
  annotate(
    geom = "text",
    x = 0,
    y = 0,
    fontface = "bold",
    label = expression(y[min]),
    size = 5,
    color = "dodgerblue"
  ) +
  # Add EC50
  annotate(
    geom = "text",
    x = .3,
    y = -.05,
    fontface = "bold",
    label = expression(paste(EC[50])),
    size = 5,
    color = "darkred"
  ) +
  annotate(
    geom = "segment",
    x = .3,
    y = 0,
    yend = .5,
    size = 1,
    color = "darkred"
  ) +
  annotate(
    geom = "segment",
    x = 0,
    xend = .3,
    y = .5,
    yend = .5,
    fontface = "bold",
    label = expression(paste(EC[50])),
    size = 1,
    color = "darkred"
  ) +
  # Add EC50 point
  geom_point(aes(x = .3, y = .5), size = 5, color = "darkred") +
  # Add 50% of y
  annotate(
    geom = "text",
    x = .1,
    y = .4,
    size = 4,
    color = "darkred",
    label = expression(paste(1 / 2(y[max] - y[min])))
  ) +
  # Add Beta parameter
  annotate(
    geom = "text",
    x = .4,
    y = .5,
    label = expression(beta)
  ) +
  annotate(
    geom = "segment",
    x = 0.25,
    xend = 0.35,
    y = .625,
    yend = .375,
    arrow = arrow(
      ends = "both",
      length = unit(0.1, "inches"),
      angle = 30
    ),
    size = 1.2
  ) +
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
  # Individual averages
  geom_point(
    data = ID,
    aes(x = 0, y = b0 + b0_i),
    size = 4,
    shape = 21,
    fill = "white"
  ) +
  # Individual points
  geom_point(data = df.rn.subset, aes(x = x, y = y), size = 3, alpha = .4) +
  # Within-individual variation ribbon (centered around each individual)
  geom_ribbon(
    aes(xmin = -2, xmax = 2, ymin = yhat - 2 * sigma, ymax = yhat + 2 * sigma),
    alpha = .2,
    linewidth = 0
  ) +
  scale_color_manual(values = colorpal) +
  scale_fill_manual(values = colorpal) +
  geom_segment(
    # Among-individual variance
    aes(x = -.1, y = -1.9, yend = 1.9),
    arrow = arrow(ends = "both", length = unit(.15, "inches")),
    color = "purple3"
  ) +
  annotate(
    "text",
    x = 1,
    y = .1,
    label = "Among-individual \n variance (VI)",
    color = "purple3",
    size = 3
  ) +
  labs(
    x = "Environmental gradient",
    y = "Phenotype",
    title = "Reaction Norm Approach"
  ) +
  theme_custom()
fig.rn
ggsave(plot = fig.rn, "outputs/figs/fig.rn.jpeg")

# Both figures side-by-side ----
fig.drc.vs.rn <- fig.drc + fig.rn + plot_annotation(tag_levels = 'A')
ggsave(
  plot = fig.drc.vs.rn,
  "outputs/figs/fig.drc.vs.rn.jpeg",
  height = 5,
  width = 10
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
