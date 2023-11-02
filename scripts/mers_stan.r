library(cmdstanr)
library(tidyverse)
library(parallel)

R0_model <- cmdstan_model(
  "stan/R0_clusters_binomial_fixedPrObs_MLE.stan",
  include_paths = "stan",
  stanc_options = list("O1")
)

prof_lik <- cmdstan_model(
  "stan/profile_likelihood_fixedPrObs_MLE.stan",
  include_paths = "stan",
  stanc_options = list("O1")
)

fit_once <- function(c_data, ptmut, pdet, psec, method = 1){
  
  cluster.data <- list(
    D = nrow(c_data),
    s = c_data$cluster_size,
    n = c_data$count,
    p = ptmut,
    p_detection = pdet * psec,
    epsilon = 1e-10,
    max_approx_it = 1E4,
    tol = 1e-10,
    method = method,
    cap = 2e3,
    max_iter = 2E4
  )
  
  opt.MLE <- R0_model$optimize(data = cluster.data)
  return(opt.MLE)
}


construct_CI <- function(par.df,
                         log_lik_at_mle){
  ## modified version of get_CI()
  
  chisq_stat_95 <- qchisq(0.95, df = 1)
  chisq_stat_90 <- qchisq(0.90, df = 1)
  chisq_stat_50 <- qchisq(0.50, df = 1)
  
  vec_param <- as.numeric(unlist(par.df[, 1]))
  
  
  ## Log-likelihood profile
  vec_log_lik <- par.df$lik
  
  ## Likelihood ratio
  vec_LR <- 2 * (log_lik_at_mle - vec_log_lik)
  
  param_within_95 <- vec_param[vec_LR < chisq_stat_95]
  param_within_90 <- vec_param[vec_LR < chisq_stat_90]
  param_within_50 <- vec_param[vec_LR < chisq_stat_50]
  
  return(
    c('lower_95' = min(param_within_95),
      'upper_95' = max(param_within_95),
      'lower_90' = min(param_within_90),
      'upper_90' = max(param_within_90),
      'lower_50' = min(param_within_50),
      'upper_50' = max(param_within_50))
  )
}

construct_BCI <- function(theta){
  ## modified version of get_CI()
  return(
    c('lower_95' = as.numeric(quantile(x = theta, probs = .025)),
      'upper_95' = as.numeric(quantile(x = theta, probs = .975)),
      'lower_90' = as.numeric(quantile(x = theta, probs = .05)),
      'upper_90' = as.numeric(quantile(x = theta, probs = .95)),
      'lower_50' = as.numeric(quantile(x = theta, probs = .25)),
      'upper_50' = as.numeric(quantile(x = theta, probs = .75)))
  )
}

fit_all <- function(c_data,
                    ptmut, pdet, psec,
                    method = 1, nboot = 200,
                    ncores = 4,
                    step = 0.01){
  
  main.fit <- fit_once(c_data, ptmut, pdet, psec, method)  
  
  LL.max <- main.fit$mle("the_likelihood")[1]
  
  Rs <- ks <- seq(step, 10, by = step)
  
  prof.data <- list(
    D = nrow(c_data),
    s = c_data$cluster_size,
    n = c_data$count,
    p = ptmut,
    p_detection = pdet * psec,
    R_max = as.numeric(main.fit$mle("R")),
    k_max = as.numeric(main.fit$mle("k")), 
    N_grid = length(Rs),
    Rs = Rs,
    ks = ks,
    epsilon = 1e-10,
    max_approx_it = 1E4,
    tol = 1e-10,
    method = method,
    cap = 2e3,
    max_iter = 2E4
  )
  cat("--- Doing profiling --- \n")
  prof.computations <- prof_lik$sample(data = prof.data,
                                       fixed_param = TRUE,
                                       iter_warmup = 1,
                                       iter_sampling = 1,
                                       chains = 1)
  
  Rstuff <- prof.computations$draws('lik_table_R', format = "array")
  R.d <- data.frame(matrix(as.vector(Rstuff), ncol = 3))
  colnames(R.d) <- c("R", "lik", "n_iter")
  R.df <- tibble(R.d)
  
  Kstuff <- prof.computations$draws('lik_table_K', format = "array")
  K.d <- data.frame(matrix(as.vector(Kstuff), ncol = 3))
  colnames(K.d) <- c("k", "lik", "n_iter")
  K.df <- tibble(K.d)
  
  CI_R <- construct_CI(par.df = R.df, log_lik_at_mle = LL.max)
  CI_k <- construct_CI(par.df = K.df, log_lik_at_mle = LL.max)
  
  ###
  
  cat("--- Doing bootstrapping --- \n")
  
  resamples <- lapply(seq_len(nboot),
                      function(i){
                        newcdata <- c_data
                        N <- sum(c_data$count)
                        newcdata$count <- as.vector(rmultinom(1,
                                                              size = N,
                                                              prob = c_data$count/N ))
                        return(subset(newcdata, count > 0)) # filtering to avoid computing a probability that does not contribute to the likelihood
                      })
  
  
  refits <- parallel::mclapply(
    1:nboot,
    function(i) fit_once(c_data = resamples[[i]],
                         ptmut, pdet, psec, method),
    mc.cores = ncores
  )
  
  ######## Now wrangle output
  
  par <- c("R", "k", "p0", "n_iters")
  ests <- lapply(refits, function(x) x$mle(par) )
  
  est.df <- tibble(as.data.frame(do.call(rbind, ests)))
  
  res.profile <- bind_cols(tibble(mle_estim = main.fit$mle(c('R', 'k')),
                                  param = c('R', 'k'),
                                  p_detect = pdet,
                                  method = "profile"),
                           bind_rows(CI_R, CI_k))
  
  res.bootstrap <- bind_cols(tibble(mle_estim = main.fit$mle(par),
                                    param = par, 
                                    p_detect = pdet,
                                    method = "bootstrap"),
                             bind_rows(construct_BCI(est.df$R),
                                       construct_BCI(est.df$k),
                                       construct_BCI(est.df$p0),
                                       construct_BCI(est.df$n_iters)
                             )
  )
  return(bind_rows(res.profile, res.bootstrap))
}

###############################################
###############################################

## Input files
file_cluster_alloc <- '../data/mers/cluster_alloc_mers.rds'
file_proba_trans_before_mut <- '../results/proba_trans_before_mut/df_p_trans_before_mut_with_uncertainty.rds'

cluster_alloc <- readRDS(file_cluster_alloc)

## Probability that transmission occurs before mutation
p_trans_before_mut <- readRDS(file_proba_trans_before_mut) %>% 
  ungroup() %>% 
  filter(pathogen == 'MERS') %>% 
  select(p_trans_before_mut) %>%
  as.numeric()

## Scenario for the proportion of infection sequenced
n_sequences <- sum(cluster_alloc$df_size_distrib$cluster_size * cluster_alloc$df_size_distrib$count)
n_suspected_cases <- 1646 
prop_cases_sequenced <- n_sequences/n_suspected_cases
vec_p_detect_cases <- c(0.5, 1.0) # Proportion of infections detected (as cases)

c_data <- cluster_alloc$df_size_distrib
ptmut <- p_trans_before_mut
pdet <- vec_p_detect_cases[1]
psec <- prop_cases_sequenced


adaptive.time <- system.time(
  mers.mle.a <- do.call(rbind, lapply(vec_p_detect_cases, function(pdet){
    fit_all(c_data = cluster_alloc$df_size_distrib,
            ptmut = p_trans_before_mut,
            pdet = pdet,
            psec = prop_cases_sequenced,
            method = 1)
  }))
)

subset(mers.mle.a, param %in% c('R', 'k'))
adaptive.time

## Run with a fixed cap
fixed.time <- system.time(
  mers.mle.f <- do.call(rbind, lapply(vec_p_detect_cases, function(pdet){
    fit_all(c_data = cluster_alloc$df_size_distrib,
            ptmut = p_trans_before_mut,
            pdet = pdet,
            psec = prop_cases_sequenced,
            method = 2)
  }))
)

mers.mle.f
fixed.time
