## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", eval = TRUE, warning = FALSE, message = FALSE, echo = TRUE, cache = FALSE, results='hide',
  # eval = FALSE,
  fig.height = 6, fig.width=8)

## ----setup--------------------------------------------------------------------
library(psborrow)

## -----------------------------------------------------------------------------
ss = set_n(ssC = 200, ssE = 180, ssExt = 400)

## -----------------------------------------------------------------------------
covset1 = set_cov(n_cat = 2, n_cont = 1, 
                  mu_int = c(0, 0.5, 0.5), mu_ext = c(0.7, 0.5, 0.9), 
                  var = c(1, 1.2, 1), cov = c(0.5, 0.7, 0.9), 
                  prob_int = c(0.45, 0.55), prob_ext = c(0.65, 0.55))

## -----------------------------------------------------------------------------
covset2 = set_cov(n_cat = 2, n_cont = 3, 
                  mu_int = 0, mu_ext = 0.2, var = 0.7, cov = 0.7, 
                  prob_int = 0.5)

## -----------------------------------------------------------------------------
cov_list <- c(covset1, covset2)
sample_cov <- simu_cov(ssObj = ss, covObj = cov_list, HR = c(0.67, 1), driftHR = c(1, 1.2), nsim = 10, seed = 47)

## ---- results='show'----------------------------------------------------------
head(sample_cov[[1]], 5)

## -----------------------------------------------------------------------------
evt = set_event(event = "pwexp", lambdaC = c(0.0135, 0.02), t_itv = 2,
                beta = c(rep(0.2,3), 0.3, rep(0.5,2), rep(2, 2)))

## ---- eval = FALSE------------------------------------------------------------
#  set_event(event = "pwexp", lambdaC = c(0.0135, 0.02), t_itv = 2,
#            beta = c(0.2, 0.3, 0.5, 2),
#            change = list(c("cov1", "*", "cov3"), c("cov2", "+", "cov3"),
#                          c("cov1", "^", "3")), keep = "cov3")

## ---- eval = FALSE------------------------------------------------------------
#  set_event(event = "weibull", shape = 0.9, lambdaC = 0.0135, beta = 0.5)

## -----------------------------------------------------------------------------
c_int = set_clin(gamma = c(2,3,16), e_itv = c(5,10),  CCOD = "fixed-first", CCOD_t = 45, etaC = c(0.02, 0.03), etaE = c(0.2, 0.3), d_itv = 2)

## -----------------------------------------------------------------------------
c_ext = set_clin(gamma = 10, CCOD = "event", CCOD_t = 150, etaC = 0.05)

## -----------------------------------------------------------------------------
sample_time <- simu_time(dt = sample_cov, eventObj = evt, clinInt = c_int, clinExt = c_ext, seed = 47)

## ---- results='show'----------------------------------------------------------
head(sample_time[[1]], 5)

## -----------------------------------------------------------------------------
pr1 <- set_prior(pred = "all", prior = "cauchy", r0 = 1, alpha = c(0, 0), sigma = 0.03)
pr2 <- set_prior(pred = "all", prior = "gamma", r0 = 1,  alpha = c(0, 0))
pr3 <- set_prior(pred = "all", prior = "no_ext", r0 = 1, alpha = 0)
pr4 <- set_prior(pred = "all", prior = "full_ext", r0 = 1, alpha = 0)

## ----eval = FALSE-------------------------------------------------------------
#  set_prior(pred = c("ps", "cov1", "cov2"), prior = "gamma", r0 = 1,  alpha = c(0, 0))

## ---- eval = FALSE------------------------------------------------------------
#  # Use no covariate. External arm will not be included in the Bayesian model
#  set_prior(pred = "none", prior = "no_ext", r0 = 1, alpha = 0)
#  
#  # Use propensity score calculated using all covariates as predictor. Hazard ratio = 1 for external and internal control arms
#  set_prior(pred = "ps", prior = "full_ext", r0 = 1, alpha = 0)
#  
#  # Use a subset of covaraites as predictors
#  set_prior(pred = c("cov2", "cov5"), prior = "no_ext", r0 = 1, alpha = 0)

## -----------------------------------------------------------------------------
pr_list = c(pr1, pr2, pr3, pr4)
res <- run_mcmc(dt = sample_time, pr_list, n.chains = 2, n.adapt = 10, n.burn = 10, n.iter = 20, seed = 47)
# res <- run_mcmc_p(dt = sample_time, pr_list, n.chains = 2, n.adapt = 100, n.burn = 100, n.iter = 200, seed = 2)

## -----------------------------------------------------------------------------
summ <- get_summary(res)

## -----------------------------------------------------------------------------
head(summ, 5)

## -----------------------------------------------------------------------------
plot_type1error(summ, driftHR = 1.2, pred = "all")
plot_power(summ, HR = 0.67, driftHR = 1, pred = "all")
plot_hr(summ, HR = 0.67, driftHR = 1.2, pred = "all")
plot_bias(summ, HR = 1, driftHR = 1.2, pred = "all")
plot_mse(summ, HR = 0.67, driftHR = 1, pred = "all")

