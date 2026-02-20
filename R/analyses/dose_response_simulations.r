# Code for Dose-Response simulations incorporating VI
library(here); library(tidyverse); library(viridis)
library(patchwork); library(truncnorm); library(ggthemes)
library(GGally); library(ggExtra); library(corrplot)
library(ggdist); library(distributional)

theme_set(theme_bw(14))

source("R/funs/dose_response_functions.R")

# Global simulation parameters ----
alpha = 100
beta = -3.5
Dose = seq(0, 100, by = 1)
NEC = 25

n_id = 10000
CV_i = .2
sigma = .08

data.frame(Dose = Dose,
           y = DR_fun_log_exp(Dose, 
                              alpha, 
                              beta, 
                              NEC)) %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_line(linewidth = 3)


# Cov alpha x NEC > 0 ---- 
rho_1_2 = 0 # r_alpha_beta
rho_1_3 = .9 # r_alpha_NEC
rho_2_3 = 0 # r_beta_NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * CV_i, beta * 0, NEC * CV_i) # 20 % CV around the mean
rho_mat = matrix(c(1, rho_1_2, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID.pos.cov = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID.pos.cov = 1:n_id)
GGally::ggpairs(ID.pos.cov[,c(1,3)])

# Simulate individual growth
df.pos.cov = ID.pos.cov %>%
  expand(nesting(ID = ID.pos.cov, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>% 
  mutate(yhat = exp(log_yhat)) %>%
  mutate(y = rlnorm(n(), log_yhat, sigma))

fig.pos.cov = df.pos.cov %>% 
  filter(ID %in% sample(unique(ID), 100)) %>% 
  ggplot(aes(y = exp(log_yhat), x = Dose, group = ID)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(alpha = .5, linewidth = .5, color = "#016392") +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.pos.cov

fig.pos.cov.sd = df.pos.cov %>% 
  summarise(y_i_mu = mean(yhat),
            y_i_sd = sd(yhat), 
            .by = c(Dose))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose)) +
  geom_line(linewidth = 2, color = "#016392") +
  labs(x = "Dose", y = "sd y response") +
  ylim(0,30) +
  theme_bw(14) +
  theme(legend.position = "none")
fig.pos.cov.sd

# Cov alpha x NEC < 0 ----
rho_1_2 = 0 # r_alpha_beta
rho_1_3 = -.9 # r_alpha_NEC
rho_2_3 = 0 # r_beta_NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * CV_i, beta * 0, NEC * CV_i) # 10 % CV around the mean
rho_mat = matrix(c(1, rho_1_2, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID.neg.cov = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID.neg.cov = 1:n_id)
GGally::ggpairs(ID.neg.cov[,c(1,3)])

# Simulate individual growth
df.neg.cov = ID.neg.cov %>%
  expand(nesting(ID = ID.neg.cov, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>%
  mutate(yhat = exp(log_yhat)) %>%
  mutate(y = rlnorm(n(), log_yhat, sigma))

fig.neg.cov = df.neg.cov %>% 
  filter(ID %in% sample(unique(ID), 100)) %>% 
  ggplot(aes(y = yhat, x = Dose, group = ID)) +
  geom_line(alpha = .5, linewidth = .5, color = "#c72e29") +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.neg.cov

fig.neg.cov.sd = df.neg.cov %>% 
  summarise(y_i_mu = mean(yhat),
            y_i_sd = sd(yhat), 
            .by = c(Dose))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose)) +
  geom_line(linewidth = 2, color = "#c72e29") +
  ylim(0,30) +
  labs(x = "Dose", y = "sd y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.neg.cov.sd

# Cov alpha x NEC = 0 ----
rho_1_2 = 0 # r_alpha_beta
rho_1_3 = 0 # r_alpha_NEC
rho_2_3 = 0 # r_beta_NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * CV_i, beta * 0, NEC * CV_i) # 10 % CV around the mean
rho_mat = matrix(c(1, rho_1_2, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID.0.cov = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID.0.cov = 1:n_id)
GGally::ggpairs(ID.0.cov[,c(1,3)])

# Simulate individual growth
df.0.cov = ID.0.cov %>%
  expand(nesting(ID = ID.0.cov, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>%
  mutate(yhat = exp(log_yhat)) %>%
  mutate(y = rlnorm(n(), log_yhat, sigma))

fig.0.cov = df.0.cov %>% 
  filter(ID %in% sample(unique(ID), 100)) %>% 
  ggplot(aes(y = yhat, x = Dose, group = ID)) +
  geom_line(alpha = .5, linewidth = .5, color = "#c72e29") +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.0.cov

fig.0.cov.sd = df.0.cov %>% 
  # mutate(yhat_std = (yhat - alpha) / (CV_i * alpha)) %>% # standardised over control
  summarise(y_i_mu = mean(yhat_std),
            y_i_sd = sd(yhat_std), 
            .by = c(Dose))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose)) +
  geom_line(linewidth = 2, color = "#c72e29") +
  # ylim(0,1) +
  ylim(0,30) +
  labs(x = "Dose", y = "sd y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.0.cov.sd


# Combine figures ----
fig.cor.pos <- 
  ggplot(ID.pos.cov,
         aes(x = NEC_i, y = alpha_i)) +
  stat_ellipse(geom = "polygon", 
               fill = "#016392", 
               alpha = .8, color = "black") +
  ylab("Individuality") +
  xlab("Sensitivity") +
  theme_minimal() + 
  theme(aspect.ratio = 1) +
  annotate("text",
           label = "r = 0.9",
           x = 18, y = 130, 
           color = "#016392", 
           size = 3.1)

fig.cor.neg <- 
  ggplot(ID.neg.cov,
         aes(x = NEC_i, y = alpha_i)) +
  stat_ellipse(geom = "polygon", 
               fill = "#c72e29", 
               alpha = .8, color = "black") +
  ylab("Individuality") +
  xlab("Sensitivity") +
  theme_minimal() + 
  theme(aspect.ratio = 1) +
  annotate("text",
           label = "r = -0.9",
           x = 32, y = 130,
           color = "#c72e29",
           size = 3.1)

fig.pos.cov <- 
  fig.pos.cov + 
  inset_element(fig.cor.pos, 
                left = 0.55,   # Start at 55% from left
                bottom = 0.55, # Start at 55% from bottom  
                right = 0.98,  # End at 98% from left
                top = 0.98,    # End at 98% from bottom
                align_to = 'panel')

fig.neg.cov <- 
  fig.neg.cov + 
  inset_element(fig.cor.neg, 
                left = 0.55,   # Start at 55% from left
                bottom = 0.55, # Start at 55% from bottom  
                right = 0.98,  # End at 98% from left
                top = 0.98,    # End at 98% from bottom
                align_to = 'panel')


fig.dr.full <- (fig.pos.cov + fig.neg.cov) /
  (fig.pos.cov.sd + fig.neg.cov.sd)

ggsave("outputs/figs/fig.dr.full.jpeg",
       fig.dr.full, 
       width = 24,
       height = 20,
       units = "cm")
