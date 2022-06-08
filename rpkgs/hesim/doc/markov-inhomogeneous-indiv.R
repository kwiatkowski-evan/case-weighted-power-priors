## ---- warning = FALSE, message = FALSE----------------------------------------
library("hesim")
library("data.table")
library("kableExtra")
library("flexsurv")
library("ggplot2")

## ---- out.width = "700px", echo = FALSE---------------------------------------
knitr::include_graphics("markov-inhomogeneous.png")

## ---- warning = FALSE, message = FALSE----------------------------------------
# Treatment strategies
strategies <- data.table(
  strategy_id = 1:2,
  strategy_name = c("Standard prosthesis", "New prosthesis")
)
n_strategies <- nrow(strategies)

# Patients
n_patients <- 1000
patients <- data.table(
  patient_id = 1:n_patients,
  gender = "Female",
  age = 60
)

# States
states <- data.table( # Non-death health states
  state_id = 1:4,
  state_name = c("Primary THR", "Sucessful primary", "Revision THR", 
                 "Successful revision")
) 
n_states <- nrow(states)

# "hesim data"
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states)
print(hesim_dat)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
print(labs)

## -----------------------------------------------------------------------------
tmat <- rbind(c(NA, 1, NA, NA, 2),
              c(NA, NA, 3, NA, 4),
              c(NA, NA, NA, 5, 6),
              c(NA, NA, 7, NA, 8),
              c(NA, NA, NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- names(labs$state_id)
tmat %>% 
  kable() %>%
  kable_styling()

## -----------------------------------------------------------------------------
ttrrPTHR <- 2 # Time to recovery rate implies mean time of 1/2 years

## -----------------------------------------------------------------------------
# 2 out of 100 patients receiving primary THR died
omrPTHR_shape1 <- 2
omrPTHR_shape2 <- 98

## -----------------------------------------------------------------------------
prob_to_rate <- function(p, t = 1){
  (-log(1 - p))/t
}

## ---- include = FALSE, results = "hide"---------------------------------------
rr_coef <- c(0.3740968, -5.490935, -0.0367022, 0.768536, -1.344474)
names(rr_coef) <- c("lngamma", "cons", "age", "male", "np1")

# Variance-covariance matrix
rr_vcov <- matrix(
  c(0.0474501^2, -0.005691, 0.000000028, 0.0000051, 0.000259,
    -0.005691, 0.207892^2, -0.000783, -0.007247, -0.000642,
    0.000000028, -0.000783, 0.0052112^2, 0.000033, -0.000111,
    0.0000051, -0.007247, 0.000033, 0.109066^2, 0.000184,
    0.000259, -0.000642, -0.000111, 0.000184, 0.3825815^2),
  ncol = 5, nrow = 5, byrow = TRUE
)
rownames(rr_vcov) <- colnames(rr_vcov) <- names(rr_coef)

## -----------------------------------------------------------------------------
print(rr_coef)
print(rr_vcov)

## ---- echo = FALSE------------------------------------------------------------
mort_tbl <- rbind(
  c(35, 45, .00151, .00099),
  c(45, 55, .00393, .0026),
  c(55, 65, .0109, .0067),
  c(65, 75, .0316, .0193),
  c(75, 85, .0801, .0535),
  c(85, Inf, .1879, .1548)
)
colnames(mort_tbl) <- c("age_lower", "age_upper", "male", "female")
mort_tbl <- data.frame(mort_tbl)

## -----------------------------------------------------------------------------
mr <- c(.0067, .0193, .0535, .1548)
mr_times <- c(0, 5, 15, 25)

## -----------------------------------------------------------------------------
ttrRTHR <- 1 # There is no rate, the time is fixed

## -----------------------------------------------------------------------------
# 4 out of 100 patients with a successful revision needed another procedure r
omrRTHR_shape1 <- 4
omrRTHR_shape2 <- 96

## -----------------------------------------------------------------------------
n_samples <- 500

## -----------------------------------------------------------------------------
matrixv <- function(v, n = NULL){
  if (length(v) == 1) v <- rep(v, n_samples) 
  m <- matrix(v)
  colnames(m) <- "cons"
  return(m)
}

## -----------------------------------------------------------------------------
transmod_coef_def <- define_rng({
  omrPTHR <- prob_to_rate(beta_rng(shape1 = omrPTHR_shape1,
                                   shape2 = omrPTHR_shape2))
  mr <- fixed(mr)
  mr_omrPTHR <- omrPTHR + mr
  rr <- multi_normal_rng(mu = rr_coef, Sigma = rr_vcov)
  rrr <- prob_to_rate(beta_rng(shape1 = 4, shape2 = 96))
  
  list(
    log_omrPTHR = matrixv(log(omrPTHR)),
    log_mr = lapply(as.list(log(mr)), matrixv),
    log_ttrrPTHR = matrixv(log(ttrrPTHR)),
    log_mr_omrPTHR = lapply(as.list(log(mr_omrPTHR)), matrixv),
    rr_shape = matrixv(rr$lngamma),
    rr_scale = as.matrix(rr[, -1,]),
    log_rrr = matrixv(log(rrr))
  )
}, n = n_samples)
transmod_coef <- eval_rng(transmod_coef_def)

## -----------------------------------------------------------------------------
head(transmod_coef$log_omrPTHR)
head(transmod_coef$rr_scale)

## -----------------------------------------------------------------------------
transmod_params <- params_surv_list(
  # 1. Primary THR:Successful primary (1:2)
  params_surv(coefs = list(rate = transmod_coef$log_ttrrPTHR), 
              dist = "fixed"),
  
  # 2. Primary THR:Death (1:5)
  params_surv(coefs = list(rate = transmod_coef$log_omrPTHR), 
              dist = "exp"), 
  
  # 3. Successful primary:Revision THR (2:3)
  params_surv(coefs = list(shape = transmod_coef$rr_shape,
                           scale = transmod_coef$rr_scale), 
              dist = "weibullPH"), 
  
  # 4. Successful primary:Death (2:5)
  params_surv(coefs = transmod_coef$log_mr,
              aux = list(time = c(0, 5, 15, 25)),
              dist = "pwexp"),
  
  # 5. Revision THR:Successful revision (3:4)
  params_surv(coefs = list(est = matrixv(ttrRTHR)),
              dist = "fixed"),
  
  # 6. Revision THR:Death (3:5)
  params_surv(coefs = transmod_coef$log_mr_omrPTHR,
              aux = list(time = c(0, 5, 15, 25)),
              dist = "pwexp"),

  # 7. Successful revision:Revision THR (4:3)
  params_surv(coefs = list(rate = transmod_coef$log_rrr),
              dist = "exp"), 
  
  # 8. Successful revision:Death (4:5)
  params_surv(coefs = transmod_coef$log_mr,
              aux = list(time = c(0, 5, 15, 25)),
              dist = "pwexp")
)

## -----------------------------------------------------------------------------
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(0, .85, .3, .75),
             se = c(0, .03, .03, .04)),
  dist = "beta"
)
head(utility_tbl)

## -----------------------------------------------------------------------------
drugcost_tbl <- stateval_tbl(
  data.table(strategy_id = rep(strategies$strategy_id, each = n_states),
             state_id = rep(states$state_id, times = n_strategies),
             est = c(394, 0, 0, 0,
                     579, 0, 0, 0)),
  dist = "fixed"
)

medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = c(0, 0, 5294, 0),
             se = c(0, 0, 1487, 0)),
  dist = "gamma",
)

## -----------------------------------------------------------------------------
transmod_data <- expand(hesim_dat, by = c("strategies", "patients"))
head(transmod_data)

## -----------------------------------------------------------------------------
transmod_data[, cons := 1]
transmod_data[, male := ifelse(gender == "Male", 1, 0)]
transmod_data[, np1 := ifelse(strategy_name == "New prosthesis", 1, 0)]

## -----------------------------------------------------------------------------
# Transition model
transmod <- create_IndivCtstmTrans(transmod_params, 
                                   input_data = transmod_data,
                                   trans_mat = tmat,
                                   clock = "forward",
                                   start_age = patients$age)

## -----------------------------------------------------------------------------
# Utility
utilitymod <- create_StateVals(utility_tbl, n = transmod_coef_def$n, 
                               hesim_data = hesim_dat)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = transmod_coef_def$n,
                                method = "starting", hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = transmod_coef_def$n,
                               hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

## -----------------------------------------------------------------------------
econmod <- IndivCtstm$new(trans_model = transmod,
                          utility_model = utilitymod,
                          cost_models = costmods)

## -----------------------------------------------------------------------------
ptm <- proc.time()
econmod$sim_disease(max_t = 60, max_age = 120)
proc.time() - ptm
head(econmod$disprog_)

## -----------------------------------------------------------------------------
econmod$sim_stateprobs(t = 0:60)

## ---- echo = FALSE, fig.width = 7, fig.height = 4-----------------------------
mean_stateprobs <- econmod$stateprobs_[, .(prob_mean = mean(prob)),
                                       by = c("strategy_id", "state_id", "t")]
mean_stateprobs[, state_name := factor(state_id,
                                       levels = 1:nrow(tmat),
                                       labels = colnames(tmat))]

ggplot(mean_stateprobs, aes(x = t, y = prob_mean, col = factor(strategy_id))) +
  geom_line() + facet_wrap(~state_name, scales = "free_y") +
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "Strategy") +
  theme(legend.position = "bottom") +
  theme_bw()

## -----------------------------------------------------------------------------
econmod$sim_qalys(dr = .015)
econmod$sim_costs(dr = .06)

## -----------------------------------------------------------------------------
ce_sim <- econmod$summarize()
summary(ce_sim, labels = labs) %>%
  format()

## -----------------------------------------------------------------------------
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = 0.015, dr_costs = .06,
                     k = seq(0, 25000, 500))
icer(cea_pw_out, labels = labs) %>%
  format(digits_qalys = 3)

## ---- eval = FALSE, echo = FALSE----------------------------------------------
#  # Rather than use a piecewise exponential distribution, approximate
#  # mortality rates with a Weibull distribution
#  fit_from_pwexp <- function(n = 10000, rate, time = mr_times){
#    sim_pwexp <- rpwexp(10000, rate = rate, time = time)
#  
#    fit_wei <- flexsurvreg(formula = Surv(sim_pwexp) ~ 1, dist = "weibull")
#    sim_wei <- rweibull(n,
#                        shape = exp(fit_wei$res.t["shape", "est"]),
#                        scale = exp(fit_wei$res.t["scale", "est"]))
#  
#    res <- fit_wei
#    attr(res, "sim_dist") <- list(wei = summary(sim_wei),
#                                  pwexp = summary(sim_pwexp))
#    return(res)
#  }
#  fit_mort_wei <- fit_from_pwexp(rate = mr)
#  fit_mort2_wei <- fit_from_pwexp(rate = mr + .02) # .02 = mean of omrPTHR
#  
#  # Sample the parameters
#  transmod_coef2 <- define_rng({
#    mr <- multi_normal_rng(mu = fit_mort_wei$res.t[, "est"],
#                           Sigma = vcov(fit_mort_wei))
#    mr_omrPTHR <- multi_normal_rng(mu = fit_mort2_wei$res.t[, "est"],
#                                   Sigma = vcov(fit_mort2_wei))
#    list(
#      mr_shape = matrixv(mr$shape),
#      mr_scale = matrixv(mr$scale),
#      mr_omrPTHR_shape = matrixv(mr_omrPTHR$shape),
#      mr_omrPTHR_scale = matrixv(mr_omrPTHR$scale)
#    )
#  }, n = n_samples)
#  transmod_coef2 <- eval_rng(transmod_coef2)
#  
#  # Set parameters of transitions
#  transition4_params <- params_surv(
#    coefs = list(shape = transmod_coef2$mr_shape,
#                 scale = transmod_coef2$mr_scale),
#    dist = "weibull")
#  
#  transition6_params <- params_surv(
#    coefs = list(shape = transmod_coef2$mr_omrPTHR_shape,
#                 scale = transmod_coef2$mr_omrPTHR_scale),
#    dist = "weibull")
#  
#  # Simulation
#  ## Construct
#  econmod2 <- econmod$clone(deep = TRUE)
#  econmod2$trans_model$params[[4]] <- transition4_params
#  econmod2$trans_model$params[[6]] <- transition6_params
#  econmod2$trans_model$params[[8]]<- transition4_params
#  
#  ## Simulate
#  econmod2$sim_disease(max_t = 60, max_age = 120)
#  econmod2$sim_qalys(dr = .015)
#  econmod2$sim_costs(dr = .06)
#  ce_sim2 <- econmod2$summarize()
#  ce_sim2$qalys[, .(mean = mean(qalys)),
#                  by = c("strategy_id")]
#  ce_sim2$costs[category == "total",
#              .(mean = mean(costs)),
#                by = c("strategy_id")]

