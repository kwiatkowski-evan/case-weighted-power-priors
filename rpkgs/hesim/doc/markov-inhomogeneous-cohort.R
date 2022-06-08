## ---- out.width = "700px", echo = FALSE---------------------------------------
knitr::include_graphics("markov-inhomogeneous.png")

## ---- warning = FALSE, message = FALSE----------------------------------------
library("hesim")
library("data.table")

# Treatment strategies
strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("Standard prosthesis", "New prosthesis"))

# Patients
ages <- seq(55, 75, 5)
age_weights <- c(.05, .1, .4, .25, .20)
sex_weights <- c(.65, .35)
weights <- rep(age_weights, times = 2) * 
           rep(sex_weights, each = length(ages))
patients <- data.table(patient_id = 1:10,
                       grp_id = 1:10,
                       sex = rep(c("Female", "Male"), each = length(ages)),
                       age = rep(ages, times = 2),
                       patient_wt = weights)

# Health states
states <- data.table(state_id = 1:4,
                     state_name = c("Primary THR", "Successful primary",
                                    "Revision THR", "Successful revision"))

# Combined data
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
print(hesim_dat)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
print(labs)

## ---- echo = FALSE------------------------------------------------------------
library("kableExtra")
tpmat <- matrix(
  c(0, "C", 0, 0, "omrPTHR",
    0, "C", "rr", 0, "mr",
    0, 0, 0, "C", "omrRTHR + mr",
    0, 0, "rrr", "C", "mr",
    0, 0, 0, 0, 1),
  nrow = 5, ncol = 5, byrow = TRUE
)
colnames(tpmat) <- rownames(tpmat) <- names(labs$state_id)
knitr::kable(tpmat) %>%
  kable_styling()

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
print(mort_tbl)

## ---- echo = FALSE------------------------------------------------------------
# Coefficients
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

## -----------------------------------------------------------------------------
params <- list(
  # Transition probabilities
  ## Operative mortality following primary THR
  omrPTHR_shape1 = 2, 
  omrPTHR_shape2 = 98,
  
  ## Revision rate for prosthesis
  rr_coef = rr_coef,
  rr_vcov = rr_vcov,
    
  ## Mortality_rates
  mr = mort_tbl,
  
  ## Operative mortality following revision THR
  omrRTHR_shape1 = 2,
  omrRTHR_shape2 = 98,
  
  ## re-revision rate
  rrr_shape1 = 4,
  rrr_shape2 = 96,
  
  # Utility
  u_mean = c(PrimaryTHR = 0, SuccessP = .85, Revision = .30, SuccessR = .75),
  u_se = c(PrimaryTHR = 0, SuccessP = .03, Revision = .03, SuccessR = .04),
  
  # Costs
  c_med_mean = c(PrimaryTHR = 0, SuccessP = 0, Revision = 5294, SuccessR = 0),
  c_med_se = c(PrimaryTHR = 0, SuccessP = 0, Revision = 1487, SuccessR = 0),
  c_Standard = 394,
  c_NP1 = 579
)

## -----------------------------------------------------------------------------
rng_def <- define_rng({
  list( 
    omrPTHR = beta_rng(shape1 = omrPTHR_shape1, shape2 = omrPTHR_shape2),
    rr_coef = multi_normal_rng(mu = rr_coef, Sigma = rr_vcov),
    mr_male = fixed(mr$male, names = mr$age_lower),
    mr_female = fixed(mr$female, names = mr$age_lower),
    omrRTHR = beta_rng(shape1 = omrRTHR_shape1, shape2 = omrRTHR_shape2),
    rrr = beta_rng(shape1 = rrr_shape1, shape2 = rrr_shape2),
    u = beta_rng(mean = u_mean, sd = u_se),
    c_med = gamma_rng(mean = c_med_mean, sd = c_med_se),
    c_Standard = c_Standard,
    c_NP1 = c_NP1
  )
}, n = 500)

## -----------------------------------------------------------------------------
input_data <- expand(hesim_dat, by = c("strategies", "patients"))
head(input_data)

## -----------------------------------------------------------------------------
transitions_def <- define_tparams({
  # Regression for revision risk
  male <- ifelse(sex == "Female", 0, 1)
  np1 <- ifelse(strategy_name == "Standard prosthesis", 0, 1)
  scale <- exp(rr_coef$cons + rr_coef$age * age + rr_coef$male * male + 
                 rr_coef$np1 * np1)
  shape <- exp(rr_coef$lngamma)
  rr <- 1 - exp(scale * ((time - 1)^shape - time^shape))
    
  # Mortality rate
  age_new <- age + time
  mr <- mr_female[["35"]] * (sex == "Female" & age_new >= 35 & age_new < 45) +
        mr_female[["45"]] * (sex == "Female" & age_new >= 45 & age_new < 55) +
        mr_female[["55"]] * (sex == "Female" & age_new >= 55 & age_new < 65) +
        mr_female[["65"]] * (sex == "Female" & age_new >= 65 & age_new < 75) +
        mr_female[["75"]] * (sex == "Female" & age_new >= 75 & age_new < 85) +
        mr_female[["85"]] * (sex == "Female" & age_new >= 85) +
        
        mr_male[["35"]] * (sex == "Male" & age_new >= 35 & age_new < 45) +
        mr_male[["45"]] * (sex == "Male" & age_new >= 45 & age_new < 55) +
        mr_male[["55"]] * (sex == "Male" & age_new >= 55 & age_new < 65) +
        mr_male[["65"]] * (sex == "Male" & age_new >= 65 & age_new < 75) +
        mr_male[["75"]] * (sex == "Male" & age_new >= 75 & age_new < 85) +
        mr_male[["85"]] * (sex == "Male" & age_new >= 85)
  
  list(
    tpmatrix = tpmatrix(
      0, C, 0,   0, omrPTHR,
      0, C, rr,  0, mr,
      0, 0, 0,   C, omrRTHR + mr,
      0, 0, rrr, C, mr,
      0, 0, 0,   0, 1)
  )
}, times = 1:60)

statevals_def <- define_tparams({
  c_prosthesis <- ifelse(strategy_name == "Standard prosthesis",
                         c_Standard,
                         c_NP1)
  list(
    utility = u,
    costs = list(
      prosthesis = c_prosthesis,
      medical = c_med
    )
  )
})

## -----------------------------------------------------------------------------
mod_def <- define_model(tparams_def = list(transitions_def, 
                                           statevals_def),
                        rng_def = rng_def, 
                        params = params)

## ----echo = FALSE-------------------------------------------------------------
ptm <- proc.time()

## ----econmod------------------------------------------------------------------
cost_args <- list(
  prosthesis = list(method = "starting"),
  medical = list(method = "wlos")
)
econmod <- create_CohortDtstm(mod_def, input_data,
                              cost_args = cost_args)

## ----echo = FALSE-------------------------------------------------------------
run_time1 <- time <- proc.time() - ptm

## ----simStateprobs------------------------------------------------------------
econmod$sim_stateprobs(n_cycles = 60)

## ----simStateprobsPlot, echo = FALSE, fig.width = 7, fig.height = 4-----------
library("ggplot2")
theme_set(theme_bw())
stateprob_summary <- econmod$stateprobs_[patient_id == 2, 
                                          .(count_mean = mean(prob) * 1000,
                                            count_lower = 1000 * quantile(prob, .025),
                                            count_upper = 1000 * quantile(prob, .975)),
                                          by = c("strategy_id", "patient_id", "state_id", "t")]
set_labels(stateprob_summary, labels = labs, 
           new_names = c("strategy_name", "state_name"))
ggplot(stateprob_summary, aes(x = t, y = count_mean)) +
  geom_line(aes(col = strategy_name)) +
  geom_ribbon(aes(x = t, ymin = count_lower, ymax = count_upper,
                  fill = strategy_name), alpha = .3) +
  facet_wrap(~state_name, scale = "free_y") +
  xlab("Year") + ylab("Count") +
  scale_fill_discrete("Strategy") + scale_color_discrete("Strategy")

## ----simStateVals-------------------------------------------------------------
econmod$sim_qalys(dr = .015, integrate_method = "riemann_right")
econmod$sim_costs(dr = .06, integrate_method = "riemann_right")

## ----cea----------------------------------------------------------------------
ce_sim <- econmod$summarize(by_grp = TRUE)
wtp <- seq(0, 25000, 500)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = 0.015, dr_costs = .06,
                     k = wtp)

summary(ce_sim, labels = labs)[grp == 2] %>%
  format()

## ----icerSubgroup-------------------------------------------------------------
icer(cea_pw_out) %>%
  format(pivot_from = "outcome")

## ----icerOverall--------------------------------------------------------------
ce_sim <- econmod$summarize(by_grp = FALSE)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = .015, dr_costs = .06,
                     k = wtp)
icer(cea_pw_out, labels = labs) %>%
  format(digits_qalys = 3)

## ----echo = FALSE-------------------------------------------------------------
ptm <- proc.time()

## -----------------------------------------------------------------------------
params_rng <- eval_rng(rng_def, params = params)

## -----------------------------------------------------------------------------
# Create new variables
tpmatrix_data <- expand(hesim_dat, by = c("strategies", "patients"),
                        times = 1:60)
tpmatrix_data[, time_stop := ifelse(is.infinite(time_stop), 61, time_stop)]
tpmatrix_data[, male := ifelse(sex == "Male", 1, 0)]
tpmatrix_data[, np1 := ifelse(strategy_name == "New prosthesis", 1, 0)]
tpmatrix_data[, age_new := age + time_stop]
tpmatrix_data[, cons := 1]

## -----------------------------------------------------------------------------
x_cons <- as.matrix(tpmatrix_data$cons)
x_rr <- as.matrix(tpmatrix_data[, .(cons, age, male, np1)])

head(x_cons)
head(x_rr)

## ---- echo = FALSE------------------------------------------------------------
mort_tbl2 <- melt(data.table(mort_tbl), 
                  id.vars = c("age_lower", "age_upper"),
                  variable.name = "sex", value.name = "rate",
                  variable.factor = FALSE)
mort_tbl2[, sex := ifelse(sex == "female", "Female", "Male")]
mort_tbl2[, age_sex := paste0(sex, "_", age_lower, "_", age_upper)]
mort_tbl2 <- mort_tbl2[, .(age_lower, age_upper, sex, age_sex, rate)]
setorderv(mort_tbl2, c("sex", "age_lower"))

## -----------------------------------------------------------------------------
print(mort_tbl2)

## -----------------------------------------------------------------------------
create_age_sex_dummies <- function(data, mort_tbl) {
  x <- matrix(0, nrow = nrow(data), ncol = nrow(mort_tbl2))
  colnames(x) <- mort_tbl2$age_sex
  for (i in 1:nrow(mort_tbl2)) {
    x[, i] <- ifelse(data$sex == mort_tbl2[i]$sex & 
                      data$age_new >= mort_tbl2[i]$age_lower & 
                      data$age_new < mort_tbl2[i]$age_upper,
                      1, 0)
  }
  return(x)
}
x_age_sex <- create_age_sex_dummies(tpmatrix_data, mort_tbl2)
head(x_age_sex)

## -----------------------------------------------------------------------------
tpmat_id <- tpmatrix_id(tpmatrix_data, 500)
head(tpmat_id)

## -----------------------------------------------------------------------------
# "Transformed" parameters
time <- tpmat_id$time_stop
omrPTHR <- c(x_cons %*% params_rng$omrPTHR)
omrRTHR <- omrPTHR
rr_scale <- exp(c(x_rr %*% t(params_rng$rr_coef[, -1])))
rr_shape <- exp(c(x_cons %*% params_rng$rr_coef$lngamma))
rr <- 1 - exp(rr_scale * ((time - 1)^rr_shape - time^rr_shape))
mr <- c(x_age_sex %*% t(cbind(params_rng$mr_female, params_rng$mr_male)))
rrr <- c(x_cons %*% params_rng$rrr)

# Transition probability matrix
tpmat <- tpmatrix(
  0, C, 0,   0, omrPTHR,
  0, C, rr,  0, mr,
  0, 0, 0,   C, omrRTHR + mr,
  0, 0, rrr, C, mr,
  0, 0, 0,   0, 1
)
tpmat[1, ]

## -----------------------------------------------------------------------------
tparams_transprobs <- tparams_transprobs(tpmat, tpmat_id)

## ----echo = FALSE-------------------------------------------------------------
run_time2 <- time <- proc.time() - ptm

## -----------------------------------------------------------------------------
transmod <- CohortDtstmTrans$new(params = tparams_transprobs)
econmod$trans_model <- transmod 
econmod$sim_stateprobs(n_cycles = 60)
econmod$sim_qalys(dr = .015, integrate_method = "riemann_right")
econmod$sim_costs(dr = .06, integrate_method = "riemann_right")

ce_sim <- econmod$summarize(by_grp = TRUE)
summary(ce_sim, labels = labs)[grp == 2] %>%
  format()

## ---- hide = TRUE-------------------------------------------------------------
data.table(
  Method = c("define_model()", "Custom tpmatrix()"),
  rbind(
    run_time1[1:3],
    run_time2[1:3]
  )
) %>%
  kable() %>%
  kable_styling()

