sigma <- .01
VR_val <- sigma^2   # residual variance = 0.01

# Exact reaction norms from TRUE individual means, anchored at the POST dose
df.ind <- df.sim.vi %>%
  select(ID, Group_n, Phase, Dose, mu) %>%
  pivot_wider(names_from = Phase, values_from = c(Dose, mu)) %>%
  filter(!is.na(mu_Pre), !is.na(mu_Post)) %>%
  mutate(
    D_g = Dose_Post - Dose_Pre,        # group dose span
    v_i = (mu_Post - mu_Pre) / D_g,    # individual slope (anchor-invariant)
    u_i = mu_Post                      # elevation at the POST dose (x = 0)
  )

var_exposed <- df.ind %>%
  group_by(Group_n) %>%
  summarise(
    n      = n(),
    D_g    = first(D_g),
    Vu     = var(u_i),                                     # var(response)
    Vv     = var(v_i),
    C_uv   = cov(u_i, v_i),
    beta_g = (mean(mu_Post) - mean(mu_Pre)) / first(D_g),
    .groups = "drop"
  ) %>%
  mutate(
    Vx      = D_g^2 / 4,
    mu_x    = -D_g / 2,              # POST anchor: covariate mean is NEGATIVE
    mu_x_sq = mu_x^2,                # = D_g^2 / 4 = Vx
    VU      = Vu,
    VS_pure = Vx * Vv,               # pure slope-variance piece
    VS_mu   = mu_x_sq * Vv,          # mu_x^2 * Vv piece
    C_cov   = 2 * mu_x * C_uv,       # signed covariance piece
    VS      = VS_pure + VS_mu + C_cov,
    VF      = beta_g^2 * Vx,
    VR      = VR_val
  )

# Control group: first-dose group's PRE only (n = 20)
df.ctrl <- df.ind %>% filter(Group_n == min(Group_n))

var_control <- tibble(
  Group_n = 0, n = nrow(df.ctrl), D_g = 0,
  Vu = var(df.ctrl$mu_Pre), Vv = NA_real_, C_uv = NA_real_, beta_g = NA_real_,
  Vx = 0, mu_x = 0, mu_x_sq = 0,
  VU = var(df.ctrl$mu_Pre), VS_pure = 0, VS_mu = 0, C_cov = 0, VS = 0,
  VF = 0, VR = VR_val
)

df.var.comp <- bind_rows(var_control, var_exposed) %>%
  arrange(Group_n) %>%
  mutate(
    VI = VU + VS,          # total among-individual variance
    VP = VI + VF + VR      # total phenotypic variance
  )
df.var.comp

df.var.comp %>%
  mutate(VI_tot = VI, VPlast = (VI-VS)) %>% 
  select(Group_n, VI) %>%
  pivot_longer(!Group_n, names_to = "type", values_to = "values") %>%
  ggplot(aes(x = Group_n, y = values, fill = type)) +
  geom_bar(
    position = "stack",
    stat = "identity",
    color = "black",
    linewidth = .4
  ) +
  scale_fill_wsj() +
  labs(x = "Dose", y = "Variance", title = "Variance partitioning by dose ")

df.var.comp %>%
  mutate(RI_tot = VI/VP, R2S = VS/VP) %>% 
  select(Group_n, RI_tot, R2S) %>%
  pivot_longer(!Group_n, names_to = "type", values_to = "values") %>%
  ggplot(aes(x = Group_n, y = values, fill = type)) +
  geom_bar(
    position = "fill",
    stat = "identity",
    color = "black",
    linewidth = .4
  ) +
  scale_fill_wsj() +
  labs(x = "Dose", y = "Variance", title = "Variance explained by dose ")
