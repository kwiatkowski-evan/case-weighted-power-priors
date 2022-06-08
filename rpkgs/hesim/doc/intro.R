## ---- include = FALSE---------------------------------------------------------
library("knitr")
opts_chunk$set(cache = TRUE)

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
library("kableExtra")
psm <- c("N-state partitioned survival model (PSM)", "`hesim::Psm`")
cdtstm <- c("Cohort discrete time state transition model (cDTSTM)", "`hesim::CohortDtstm`")
ictstm <- c("Individual-level continuous time state transition model (iCTSTM)", "`hesim::IndivCtstm`")
tbl <- rbind(psm, cdtstm, ictstm) 
colnames(tbl) <- c("Economic model", "R6 class")
knitr::kable(tbl, row.names = FALSE)  %>% # Pipe imported from kableExtra
  kableExtra::kable_styling()

## ---- out.width = "600px", echo = FALSE---------------------------------------
knitr::include_graphics("econ-eval-process-hesim.png")

## ----warning = FALSE, message = FALSE-----------------------------------------
library("hesim")
library("data.table")

# Treatment strategies
strategies <- data.table(strategy_id = c(1, 2),
                         strategy_name  = c("SOC", "New"))

# Patients
n_patients <- 1000
patients <- data.table(patient_id = 1:n_patients,
                       age = rnorm(n_patients, mean = 45, sd = 7),
                       female = rbinom(n_patients, size = 1, prob = .51))
patients[, grp_id := ifelse(female == 1, 1, 2)]
patients[, grp_name := ifelse(female == 1, "Female", "Male")]

# (Non-death) health states
states <- data.table(state_id = c(1, 2),
                     state_name = c("Stage 1", "Stage 2")) 

# Transitions
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c("Stage 1", "Stage 2", "Death")
transitions <- create_trans_dt(tmat)
transitions[, trans := factor(transition_id)]

# Combining
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states,
                        transitions = transitions)
print(hesim_dat)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
print(labs)

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
tbl <- rbind(
  c("`hesim::CohortDtstm`",
    "Custom", 
    "`hesim::tparams_transprobs`", 
    "`msm::msm`"),
  c("`hesim::CohortDtstm`", 
    "Multinomial logistic regressions",
    "`hesim::params_mlogit_list`",
    "`hesim::multinom_list`"),
  c("`hesim::Psm`", 
    "Independent survival models", 
    "`hesim::params_surv_list`", 
    "`hesim::flexsurvreg_list`"),
  c("`hesim::IndivCtstm`", 
    "Multi-state model (joint likelihood)", 
    "`hesim::params_surv`", 
    "`flexsurv::flexsurvreg`"),
  c("`hesim::IndivCtstm`", 
    "Multi-state model (transition-specific)", 
    "`hesim::params_surv_list`", 
    "`hesim::flexsurvreg_list`")  
)
colnames(tbl) <- c("Economic model (R6 class)", "Statistical model", 
                   "Parameter object", "Model object")
knitr::kable(tbl, row.names = FALSE) %>%
  kableExtra::kable_styling() %>%
  kableExtra::collapse_rows(columns = 1, valign = "top")

## ---- message = FALSE, warning = FALSE----------------------------------------
library("flexsurv")
mstate_data <- data.table(mstate3_exdata$transitions)
mstate_data[, trans := factor(trans)]
fit_wei <- flexsurv::flexsurvreg(Surv(years, status) ~ trans + 
                                                       factor(strategy_id):trans +
                                                       age:trans + 
                                                       female: trans +
                                                       shape(trans), 
                                 data = mstate_data, 
                                 dist = "weibull")

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
tbl <- rbind(
  c("Predicted means", 
    "`hesim::tparams_mean`", 
    "`hesim::stateval_tbl`"),
  c("Linear model", 
    "`hesim::params_lm`", 
    "`stats::lm`")
)
colnames(tbl) <- c("Statistical model", "Parameter object", "Model object")
knitr::kable(tbl, row.names = FALSE) %>%
  kableExtra::kable_styling() %>%
    kableExtra::collapse_rows(columns = 1, valign = "top")

## -----------------------------------------------------------------------------
# Utility
utility_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = mstate3_exdata$utility$mean,
             se = mstate3_exdata$utility$se),
  dist = "beta"
)

# Costs
drugcost_tbl <- stateval_tbl(
  data.table(strategy_id = strategies$strategy_id,
             est = mstate3_exdata$costs$drugs$costs),
  dist = "fixed"
)

medcost_tbl <- stateval_tbl(
  data.table(state_id = states$state_id,
             mean = mstate3_exdata$costs$medical$mean,
             se = mstate3_exdata$costs$medical$se),
  dist = "gamma"
)

## -----------------------------------------------------------------------------
print(utility_tbl)
print(drugcost_tbl)
print(medcost_tbl)

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
dtstm <- c("`hesim::CohortDtstm`", "`hesim::CohortDtstmTrans`",
         "`hesim::StateVals`", "`hesim::StateVals`")
psm <- c("`hesim::Psm`", "`hesim::PsmCurves`",
         "`hesim::StateVals`", "`hesim::StateVals`")
ictstm <- c("`hesim::IndivCtstm`", "`hesim::IndivCtstmTrans`",
         "`hesim::StateVals`", "`hesim::StateVals`")
tbl <- rbind(dtstm, psm, ictstm)
colnames(tbl) <- c("Economic model", "Disease model", "Utility model", "Cost model(s)")
knitr::kable(tbl, row.names = FALSE) %>%
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
n_samples <- 1000

## ----warning = FALSE, message = FALSE-----------------------------------------
transmod_data <- expand(hesim_dat, 
                        by = c("strategies", "patients", "transitions"))
head(transmod_data)

## -----------------------------------------------------------------------------
transmod <- create_IndivCtstmTrans(fit_wei, transmod_data,
                                   trans_mat = tmat, n = n_samples,
                                   uncertainty = "normal")
class(transmod)

## -----------------------------------------------------------------------------
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples, hesim_data = hesim_dat)
costmods <- list(drugs = drugcostmod,
                 medical = medcostmod)

## -----------------------------------------------------------------------------
ictstm <- IndivCtstm$new(trans_model = transmod,
                         utility_model = utilitymod,
                         cost_models = costmods)

## ----echo = FALSE, message = FALSE, warning = FALSE---------------------------
cdtstm_methods <- c("`hesim::CohortDtstm`", "$sim_stateprobs()", "$sim_qalys()", "$sim_costs()")
psm_methods <- c("`hesim::Psm`", "$sim_survival() and $sim_stateprobs()", "$sim_qalys()", "$sim_costs()")
ictstm_methods <- c("`hesim::IndivCtstm`", "$sim_disease() and $sim_stateprobs()", "$sim_qalys()", "$sim_costs()")
tbl <- rbind(cdtstm_methods, psm_methods, ictstm_methods)
colnames(tbl) <- c("Economic model (R6 class)", "Disease progression", "QALYs", "Costs")
knitr::kable(tbl, row.names = FALSE) %>%
  kableExtra::kable_styling()

## -----------------------------------------------------------------------------
ictstm$sim_disease()
head(ictstm$disprog_)

## -----------------------------------------------------------------------------
ictstm$sim_stateprobs(t = c(0:10))
head(ictstm$stateprobs_)

## -----------------------------------------------------------------------------
# QALYs
ictstm$sim_qalys(dr = .03)
head(ictstm$qalys_)

# Costs
ictstm$sim_costs(dr = .03)
head(ictstm$costs_)

## -----------------------------------------------------------------------------
ce <- ictstm$summarize(by_grp = FALSE)
print(ce)

## -----------------------------------------------------------------------------
cea_out <- cea(ce, dr_qalys = .03, dr_costs = .03)
cea_pw_out <- cea_pw(ce, dr_qalys = .03, dr_costs = .03, comparator = 1)

## ----ceac_plot, warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4----
library("ggplot2")
plot_ceac(cea_out, labels = labs) +
  theme_minimal()

