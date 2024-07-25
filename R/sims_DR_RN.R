# Code for Dose-Response simulations incorporating VI
library(tidyverse)

# Parameter list
alpha = 100
beta = .08
Dose = seq(0, 100, by = 1)
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
