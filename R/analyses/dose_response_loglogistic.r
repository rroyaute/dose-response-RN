library(tidyverse)
library(drc)

sigma <- .2
c <- 0
d <- 1
b <- 3
e <- .2
n<- 100
x <- seq(0,1, by = .01)

mu <- c + (d - c) / (1 + exp(b * log(x/e)))
mu <- d / (1 + exp(b * log(x/e)))

plot(x, mu, type = "l")


x <- seq(0,1, by = .15)

mu <- c + (d - c) / (1 + exp(b * log(x/e)))
mu <- d / (1 + exp(b * log(x/e)))
plot(x, mu)

set.seed(42)
df <- data.frame(x = seq(0,1, by = .01)) %>% 
  mutate(mu = d / (1 + exp(b * log(x/e)))) %>% 
  mutate(y = rlnorm(n(), log(mu), sigma)) %>% 
  mutate(y_low = qlnorm(p = .025, log(mu), sigma)) %>% 
  mutate(y_up = qlnorm(p = .975, log(mu), sigma))

mod.drm <- drm(y~x, data = df, fct = LL.3())
summary(mod.drm)
confint(mod.drm, "e")


df %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line(aes(x = x, y = mu), linewidth = 1, color = "dodgerblue") +
  geom_point() +
  geom_ribbon(aes(ymin = y_low, ymax = y_up), alpha = .4, fill = "dodgerblue") +
  geom_point(aes(x = 0.1801391, y = .5), color = "tomato2") +
  geom_errorbarh(aes(y = .5, x = 0.1801391, 
                     xmin = 0.1664766, xmax = 0.1938017,
                     width = 0), color = "tomato2") +
  theme_bw()


