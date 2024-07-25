# Code for Dose-Response simulations incorporating VI
library(tidyverse); library(viridis)

# Parameter list
alpha = 100
beta = .08
Dose = seq(0, 100, by = 10)
NEC = 25
sigma = 1
nreps = 10

DR_fun = function(Dose, alpha , beta, NEC){
  yhat = ifelse(Dose < NEC, alpha, alpha * exp(- beta * (Dose - NEC)))
  return(yhat)
}

plot(Dose, DR_fun(Dose, alpha , beta, NEC))

set.seed(42)
df = crossing(rep = 1:nreps, Dose) %>% 
  mutate(yhat = DR_fun(Dose, alpha, beta, NEC)) %>% 
  mutate(y = rnorm(n(), yhat, sigma))
df %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_fun, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 1) +
  ylim(0, 100) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)

# log version
DR_fun_log = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - beta * (Dose - NEC))
  return(log_yhat)
}

DR_fun_log_exp = function(Dose, alpha , beta, NEC){
  log_yhat = ifelse(Dose < NEC, log(alpha), 
                    log(alpha) - beta * (Dose - NEC))
  yhat = exp(log_yhat)
  return(yhat)
}

plot(Dose, exp(DR_fun_log(Dose, alpha , beta, NEC)))

sigma = .08

set.seed(42)
df = crossing(rep = 1:nreps, Dose) %>% 
  mutate(log_yhat = DR_fun_log(Dose, alpha, beta, NEC)) %>% 
  mutate(y = rlnorm(n(), log_yhat, sigma))
df %>% 
  ggplot(aes(y = y, x = Dose)) +
  geom_point(alpha = .2, size = 2.5) +
  geom_function(fun = DR_fun_log_exp, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 1) +
  # ylim(0, 100) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14)

# Simulate individuals with different sensitivities
Dose = seq(0, 100, by = 10)
n_id = 8
rho = 0 # Suppose all parameters are independent

Mu = c(alpha, beta, NEC)
sigmas = c(alpha * .1, beta * .1, NEC * .1) # 10 % CV around the mean
rho_mat = matrix(c(1, rho, rho,
                   rho, 1, rho,
                   rho, rho, 1), 
                 nrow = 3)

Sigma = diag(sigmas) %*% rho_mat %*% diag(sigmas)

set.seed(42)
ID = MASS::mvrnorm(n_id, Mu, Sigma) %>% 
  data.frame() %>% 
  set_names("alpha_i", "beta_i", "NEC_i") %>% 
  mutate(ID = 1:n_id)

# Simulate individual growth
df = ID %>%
  expand(nesting(ID, alpha_i, beta_i, NEC_i), 
                Dose = Dose) %>%
  mutate(log_yhat = DR_fun_log(Dose, alpha, beta, NEC)) %>% 
  mutate(y = rlnorm(n(), log_yhat, sigma))

# Sort individuals by alpha_i and create a mapping
ID_sorted <- ID %>%
  arrange(desc(alpha_i)) %>%
  mutate(color_index = row_number())

# Create a viridis color palette for n_id colors
viridis_colors <- viridis(n_id)

fig = df %>% 
  left_join(ID_sorted, by = "ID") %>%  # Join to get the color_index
  ggplot(aes(y = y, x = Dose, group = ID)) +
  geom_point(aes(color = factor(ID)), alpha = .8, size = 3) +
  # Add individual trends
  map(1:n_id, ~geom_function(
    fun = DR_fun_log_exp,
    args = list(
      alpha = ID_sorted$alpha_i[.x],
      beta = ID_sorted$beta_i[.x],
      NEC = ID_sorted$NEC_i[.x]
    ),
    color = viridis_colors[ID_sorted$color_index[.x]],  # Use sorted index for colors
    alpha = 0.5, 
    size = 1
  )) +
  # Add population average trend
  geom_function(fun = DR_fun_log_exp, 
                args = list(alpha = alpha, 
                            beta = beta, 
                            NEC = NEC),
                color = "red", size = 3) +
  scale_color_viridis_d(guide = guide_legend(title = "ID (ordered by alpha_i)")) +
  labs(x = "Dose", y = "y response") +
  theme_bw(14) +
  theme(legend.position = "none")
fig
