// generated with brms 2.21.6
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_alpha;  // number of population-level effects
  matrix[N, K_alpha] X_alpha;  // population-level design matrix
  int<lower=1> K_beta;  // number of population-level effects
  matrix[N, K_beta] X_beta;  // population-level design matrix
  int<lower=1> K_NEC;  // number of population-level effects
  matrix[N, K_NEC] X_NEC;  // population-level design matrix
  // covariates for non-linear functions
  vector[N] C_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_alpha_1;
  vector[N] Z_1_beta_2;
  vector[N] Z_1_NEC_3;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector<lower=0>[K_alpha] b_alpha;  // regression coefficients
  vector<lower=0>[K_beta] b_beta;  // regression coefficients
  vector<lower=0,upper=100>[K_NEC] b_NEC;  // regression coefficients
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_alpha_1;
  vector[N_1] r_1_beta_2;
  vector[N_1] r_1_NEC_3;
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_alpha_1 = r_1[, 1];
  r_1_beta_2 = r_1[, 2];
  r_1_NEC_3 = r_1[, 3];
  lprior += normal_lpdf(b_alpha | 100, 20)
    - 1 * normal_lccdf(0 | 100, 20);
  lprior += exponential_lpdf(b_beta | 10);
  lprior += uniform_lpdf(b_NEC | 0, 100)
    - 1 * log_diff_exp(uniform_lcdf(100 | 0, 100), uniform_lcdf(0 | 0, 100));
  lprior += exponential_lpdf(sigma | 1);
  lprior += exponential_lpdf(sd_1 | 1);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 4);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] nlp_alpha = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_beta = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] nlp_NEC = rep_vector(0.0, N);
    // initialize non-linear predictor term
    vector[N] mu;
    nlp_alpha += X_alpha * b_alpha;
    nlp_beta += X_beta * b_beta;
    nlp_NEC += X_NEC * b_NEC;
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_alpha[n] += r_1_alpha_1[J_1[n]] * Z_1_alpha_1[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_beta[n] += r_1_beta_2[J_1[n]] * Z_1_beta_2[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      nlp_NEC[n] += r_1_NEC_3[J_1[n]] * Z_1_NEC_3[n];
    }
    for (n in 1:N) {
      // compute non-linear predictor values
      mu[n] = (log(nlp_alpha[n]) - nlp_beta[n] * (C_1[n] - nlp_NEC[n]) * (C_1[n] > nlp_NEC[n]));
    }
    target += lognormal_lpdf(Y | mu, sigma);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}

