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
}
parameters {
  real<lower=0, upper=10> R;
  real<lower=0, upper=10> k; 
}
transformed parameters{
  real R0 = R * p;
  real lgL = compute_log_L({R0, k, p_detection});
  real lp0 = prob_s0(R0, k, p_detection, tol, max_iter);
  // real lp0 = log_marg_prob_s(0, R0, k, p_detection,
  //                            epsilon, max_approx_it, lgL)[1];
  array[2] real the_likelihood = compute_likelihood(D, n, s,
                               R0, k, p_detection,
                               lp0, lgL, 
                               epsilon, max_iter, cap,
                               method);
 real n_iters = the_likelihood[2];                               
}
model{
  target += the_likelihood[1];
}
generated quantities{
  real<lower=0, upper=1> p0 = exp(lp0);
}
