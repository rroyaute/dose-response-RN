# Code for Dose-Response simulations incorporating VI
library(here); library(tidyverse); library(viridis); library(brms); 
library(tidybayes); library(patchwork); library(truncnorm)
library(SBC); library(bayestestR); library(ggthemes)
library(GGally); library(ggExtra); library(corrplot)

bayesplot::color_scheme_set("darkgray")
theme_set(theme_bw(14))

source("R/funs/script_funs.R")

# Simulate DGP ----
I = 1000
CVi = .5
nsim = 100

# Test case with extreme individual variation
df = data_generator(ID_pars = "all", CVi = CVi, I = I)
vars = df$variables
df = df$generated 

fig.all = df %>% 
  ggplot(aes(y = y, x = Dose)) +
  # Add individual trends
  map(1:I, ~geom_function(
    fun = DR_logRN_fun,
    args = list(
      # Dose = seq(0, 100, length.out = 200),
      Rmin = log(exp(vars$b_Rmin_Intercept)),
      Rmax = log(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))[.x],
      beta = (vars$b_beta_Intercept + vars$r_ID__beta_Intercept)[.x],
      NEC =( vars$b_NEC_Intercept + vars$r_ID__NEC_Intercept)[.x]),
    alpha = 0.2, 
    linewidth = .25)) +
  # Add population average trend
  geom_function(fun = DR_logRN_fun, 
                args = list(Rmin = log(exp(vars$b_Rmin_Intercept)), 
                            Rmax = log(exp(vars$b_Rmax_Intercept)), 
                            beta = vars$b_beta_Intercept, 
                            NEC = vars$b_NEC_Intercept),
                color = "red", linewidth = 1.25) +
  labs(x = "Dose", y = "Response") +
  theme_bw(14) +
  theme(legend.position = "none") #+
  # ggtitle("All parameters", 
  #         subtitle = "with sd = 0.5")
fig.all
ggsave("outputs/figs/fig.all.jpeg",
       fig.all, 
       width = 18, 
       scale = 1.2,
       units = "cm")

hist(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
mu_Rmax = mean(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
sd_Rmax = sd(exp(vars$b_Rmax_Intercept + vars$r_ID__Rmax))
CV = sd_Rmax / mu_Rmax * 100
CV

hist(vars$b_NEC_Intercept + vars$r_ID__NEC)
mu_NEC = mean(vars$b_NEC_Intercept + vars$r_ID__NEC)
sd_NEC = sd(vars$b_NEC_Intercept + vars$r_ID__NEC)
CV = sd_NEC / mu_NEC * 100
CV

hist(exp(vars$b_beta_Intercept + vars$r_ID__beta))
mu_beta = mean(exp(vars$b_beta_Intercept + vars$r_ID__beta))
sd_beta = sd(exp(vars$b_beta_Intercept + vars$r_ID__beta))
CV = sd_beta / mu_beta * 100
CV


# Linear regression by dose class ----
## Sample ID by dose pairs -----
D_eps <- .00001 # for log-scale plotting
df <- df %>% mutate(Dose = case_when(Dose == 0 ~ D_eps, .default = Dose))
df.1 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 10)
df.2 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 20)
df.3 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 30)
df.4 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 40)
df.5 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 50)
df.6 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 60)
df.7 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 70)
df.8 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 80)
df.9 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 90)
df.10 <- df %>% 
  filter(ID %in% sample(unique(ID), 30)) %>% 
  filter(Dose == D_eps | Dose == 100)

## plot ----
p1 <- df.1 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p2 <- df.2 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p3 <- df.3 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p4 <- df.4 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p5 <- df.5 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p6 <- df.6 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p7 <- df.7 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p8 <- df.8 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p9 <- df.9 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)
p10 <- df.10 %>% 
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_line(alpha = .5) +
  xlim(0, 100)

## Stitch figure ----
fig_lin_rn <- 
  (p1 | p2 | p3 | p4 | p5) /
  (p6 | p7 | p8 | p9 | p10) 

ggsave("outputs/figs/fig_lin_rn.jpeg",
       fig_lin_rn, 
       width = 18, 
       scale = 1.2,
       units = "cm")
  



# Manual version ----
## Cov alpha x NEC > 0 ---- 
alpha = 100
beta = .08
Dose = seq(0, 100, by = 10)
NEC = 25

Dose = seq(0, 100, by = 1)
n_id = 10000
CV_i = .2
sigma = .08
rho_1_2 = 0 # r_alpha_beta
rho_1_3 = .9 # r_alpha_NEC
rho_2_3 = 0 # r_beta_NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * CV_i, beta * CV_i, NEC * CV_i) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID.pos.cov = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID.pos.cov = 1:n_id)
GGally::ggpairs(ID.pos.cov[,1:3])

# Simulate individual growth
df.pos.cov = ID.pos.cov %>%
  expand(nesting(ID = ID.pos.cov, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>% 
  mutate(yhat = exp(log_yhat)) %>%
  mutate(y = rlnorm(n(), log_yhat, sigma))

fig.pos.cov = df.pos.cov %>% 
  ggplot(aes(y = exp(log_yhat), x = Dose, group = ID)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(alpha = .5, linewidth = .3, color = "#016392") +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.pos.cov

fig.pos.cov.sd = df.pos.cov %>% 
  summarise(y_i_mu = mean(yhat),
            y_i_sd = sd(yhat), 
            .by = c(Dose))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(linewidth = 2) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.pos.cov.sd

## Cov alpha x NEC < 0 ----
alpha = 100
beta = .08
Dose = seq(0, 100, by = 10)
NEC = 25

Dose = seq(0, 100, by = 1)
n_id = 10000
CV_i = .2
sigma = .08
rho_1_2 = 0 # r_alpha_beta
rho_1_3 = -.9 # r_alpha_NEC
rho_2_3 = 0 # r_beta_NEC


Mu = c(alpha, beta, NEC)
sigmas = c(alpha * CV_i, beta * CV_i, NEC * CV_i) # 10 % CV around the mean
rho_mat = matrix(c(1, 0, rho_1_3,
                   rho_1_2, 1, rho_2_3,
                   rho_1_3, rho_2_3, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID.neg.cov = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID.neg.cov = 1:n_id)
GGally::ggpairs(ID.neg.cov[,1:3])

# Simulate individual growth
df.neg.cov = ID.neg.cov %>%
  expand(nesting(ID = ID.neg.cov, alpha_i, beta_i, NEC_i), 
         Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha_i, beta_i, NEC_i)) %>%
  mutate(yhat = exp(log_yhat)) %>%
  mutate(y = rlnorm(n(), log_yhat, sigma))

fig.neg.cov = df.neg.cov %>% 
  ggplot(aes(y = yhat, x = Dose, group = ID)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(alpha = .5, linewidth = .3, color = "#c72e29") +
labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.neg.cov

fig.neg.cov.sd = df.neg.cov %>% 
  summarise(y_i_mu = mean(yhat),
            y_i_sd = sd(yhat), 
            .by = c(Dose))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(linewidth = 2) +
labs(x = "Dose", y = "sd y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig.neg.cov.sd

## Combine figures ----
fig.dr.ID <- fig.pos.cov / fig.neg.cov
ggsave("outputs/figs/fig.dr.ID.jpeg",
       fig.dr.ID, 
       units = "cm")

fig.pos.cov.sd / fig.neg.cov.sd

df.cov <- rbind(df.pos.cov, df.neg.cov)
df.cov$cov <- c(rep(">0", length(df.pos.cov$Dose)),
                rep("<0", length(df.neg.cov$Dose)))
fig.cov.sd = df.cov %>% 
  summarise(y_i_mu = mean(yhat),
            y_i_sd = sd(yhat), 
            .by = c(Dose, cov))  %>% 
  ggplot(aes(y = y_i_sd, x = Dose, group = cov, color = cov)) +
  # geom_point(alpha = .1, size = 3) +
  geom_line(linewidth = 2) +
  scale_color_wsj() +
  labs(x = "Dose", y = "sd y response") +
  theme_bw(14)
fig.cov.sd
ggsave("outputs/figs/fig.cov.sd.jpeg",
       fig.cov.sd, 
       units = "cm")

p.pos.cov <- ID.pos.cov %>% 
  ggplot(aes(y = alpha_i, x = NEC_i)) +
  geom_point(alpha = .1, size = 3, color = "#016392") +
  ylab("alpha") + xlab("NEC") +
  theme(aspect.ratio=1)
#p.pos.cov <- ggMarginal(p.pos.cov)
p.neg.cov <- ID.neg.cov %>% 
  ggplot(aes(y = alpha_i, x = NEC_i)) +
  geom_point(alpha = .1, size = 3, color = "#c72e29") +
  ylab("alpha") + xlab("NEC") +
  theme(aspect.ratio=1)
#p.neg.cov <- ggMarginal(p.neg.cov)
p.cov <- (p.pos.cov / p.neg.cov)

fig.cov.sd_full <- fig.cov.sd + 
  inset_element(p.cov,
                left = 0.6, bottom = 0.6, 
                right = 1, top = .9)
ggsave("outputs/figs/fig.cov.sd_full.jpeg",
       fig.cov.sd_full, 
       width = 24, 
       scale = 1.2,
       units = "cm")
