functions {
  #include pmfs.stan
  #include infiniteAdaptive.stan
  #include finiteSumFixedCap.stan
  #include likelihood.stan
}
data {
  int<lower=1> D; // number of unique sequence cluster sizes
  array[D] int<lower=0> s; // sequence cluster sizes
  vector[D] n; // size counts
  real<lower=0, upper=1> p; // Pr(mut)
  real<lower=0, upper=1> p_detection;
  real<lower=0> tol;
  real<lower=0> epsilon;
  int<lower=0> max_iter;
  int<lower=0> cap;
  int<lower=0> max_approx_it; // 10000 usually, used for P(S = 0)
  int<lower=1, upper=2> method;
  real<lower=0, upper= 10> R_max;
  real<lower=0> k_max; 
  int<lower=0> N_grid;
  array[N_grid] real<lower=0> Rs;
  array[N_grid] real<lower=0> ks;
}
generated quantities{
  array[N_grid, 3] real lik_table_R;
  array[N_grid, 3] real lik_table_K;
  array[2] real temp_r;
  array[2] real temp_k;
  real ll0_r;
  real lgl_r;
  real ll0_k;
  real lgl_k;
  for(j in 1:N_grid){
    ll0_r = prob_s0(Rs[j] * p, k_max, p_detection, tol, max_iter);
    lgl_r = compute_log_L({Rs[j] * p, k_max, p_detection});
    temp_r = compute_likelihood(D, n, s,
                               Rs[j] * p, k_max, p_detection,
                               ll0_r, lgl_r, 
                               epsilon, max_iter, cap,
                               method);
    lik_table_R[j, ] = {Rs[j], temp_r[1], temp_r[2]};
    //
    ll0_k = prob_s0(R_max * p, ks[j], p_detection, tol, max_iter);
    lgl_k = compute_log_L({R_max * p , ks[j], p_detection});
    temp_k = compute_likelihood(D, n, s,
                               R_max * p, ks[j], p_detection,
                               ll0_k, lgl_k, 
                               epsilon, max_iter, cap,
                               method);
    lik_table_K[j, ] = {ks[j], temp_k[1], temp_k[2]};
  }
}

