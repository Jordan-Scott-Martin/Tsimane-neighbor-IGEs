
functions {
  //Function to compute the squared exponential covariance matrix
  //Thanks to Nicholas Clark (https://rpubs.com/NickClark47/stan_geostatistical)
  matrix L_cov_exp_quad_dis(
          int N,
          matrix distances,
                real alpha,
                real rho,
                real delta) {
    matrix[N, N] Sigma;
    Sigma = square(alpha) * exp(-0.5 * (square(distances / rho))) +
            diag_matrix(rep_vector(delta, N));
    return cholesky_decompose(Sigma);
  }
}

data {
  int<lower=1> N; //n of obs
  int<lower=1> I; //n of focal women
  int<lower=1> F; //n of birth fathers
  int<lower=1> M; //n of moms (for focal woman)
  int<lower=1> C; //n of coms
  int<lower=1> Y; //n of years
  
  array[N] int focal_id; //focal female index
  array[N] int maternal_id; //mother ID index 
  array[N] int com_id; //community index
  array[N] int yearf_id; //year of first birth index
  array[N] int yearl_id; //year of first birth index
  
  //standardized covariates
  vector[N] yearobs; //length of fertility window
  vector[N] age; //focal female age
  
  int social_dim; //maximum possible social partners
  array[N, social_dim] int social_id; //IDs for social effects
  array[N] int social_n; //number of partners/obs
  
  int father_dim; //maximum possible social partners
  array[N, father_dim] int father_id; //IDs for social effects
  array[N] int father_n; //number of partners/obs
  
  matrix[I, I] A; //relatedness matrix
  matrix[C, C] S; //spatial matrix
  vector[N] fert; //total fertility/community
}
transformed data {
  cholesky_factor_corr[I] LA = cholesky_decompose(A);
  real delta = 1e-9;
}
parameters {
  real W_0; //intercept
  
  //mean fertility effects
  real b_1; //years of sampling offset
  real b_2; //age
  real b_3; //age * age
  real b_4; //offset unknown neighbors
  
  //Gaussian process
  real<lower=0> rho;
  real<lower=0> alpha;
  
  //random effect scales
  real<lower=0> sd_P; //adjusted obs variance (var_G + var_E)
  real<lower=0> sd_M; //maternal effect (on focal woman)
  real<lower=0> sd_Y; //years of observed births
  real<lower=0, upper=1> h2; //var_G / (var_P)
  
  //standardized random effects (z-scores)
  vector[I] z_G;
  vector[C] z_C;
  vector[M] z_M;
  vector[Y] z_Y;
  
}
transformed parameters {
  //scale standardized random effects
  real sd_G = sqrt((sd_P * sd_P) * h2); //h2*var_P = var_G
  real sd_E = sqrt((sd_P * sd_P) * (1-h2));
  vector[I] W_G = LA * z_G * sd_G;
  vector[M] W_M = z_M * sd_M; //maternal effect
  vector[Y] W_Y = z_Y * sd_Y; //year effect
  
  //Calculate the spatial Gaussian process function
  matrix[C, C] L_chol = L_cov_exp_quad_dis(C, S, alpha, rho, delta);
  vector[C] W_C = L_chol * z_C;
  
  vector[N] Y_sum; //summed year effect over births
  for (n in 1 : N) {
    Y_sum[n] = sum(W_Y[yearf_id[n] : yearl_id[n]]);
  }
}
model {
  //initialize linear predictor
  vector[N] W = W_0 + W_G[focal_id] + W_C[com_id] + W_M[maternal_id] + Y_sum +
                b_1 * yearobs + b_2 * age + b_3 * (age .* age);
  
  //likelihood for observed values
  target += normal_lpdf(fert | W, sd_E);
  
  //priors
  W_0 ~ normal(3, 1);
  b_1 ~ normal(0, 1);
  b_2 ~ normal(0, 1);
  b_3 ~ normal(0, 1);
  b_4 ~ normal(0, 1);
  
  rho ~ inv_gamma(5,5);
  alpha ~ exponential(2);
  
  z_G ~ std_normal();
  z_C ~ std_normal();
  z_M ~ std_normal();
  z_Y ~ std_normal();
  
  sd_P ~ exponential(2);
  h2 ~ normal(1.2,1.2);
  sd_M ~ exponential(2);
  sd_Y ~ exponential(2);
}
