## ---- include = FALSE---------------------------------------------------------
ggplot2::theme_set(ggplot2::theme_bw())

## ---- out.width = "500px", echo = FALSE---------------------------------------
knitr::include_graphics("reversible-illness-death.png")

## -----------------------------------------------------------------------------
tmat <- rbind(c(NA, 1, 2),
              c(3, NA, 4),
              c(NA, NA, NA))
colnames(tmat) <- rownames(tmat) <- c("Healthy", "Sick", "Death")
print(tmat)

## ---- warning = FALSE, message = FALSE----------------------------------------
library("hesim")
library("data.table")
strategies <- data.table(strategy_id = c(1, 2),
                         strategy_name = c("SOC", "New"))
n_patients <- 1000
patients <- data.table(patient_id = 1:n_patients,
                       age = rnorm(n_patients, mean = 45, sd = 7),
                       female = rbinom(n_patients, size = 1, prob = .51))
states <- data.table(state_id = c(1, 2),
                     state_name = rownames(tmat)[1:2]) # Non-death health states
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients, 
                        states = states)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
labs$transition_id <- c("Healthy-> Sick" = 1, 
                        "Healthy -> Death" = 2,
                        "Sick -> Healthy" = 3,
                        "Sick -> Death" = 4)
print(labs)

## ---- warning = FALSE, message = FALSE----------------------------------------
library("flexsurv")
n_trans <- max(tmat, na.rm = TRUE) # Number of transitions
wei_fits_cr <- vector(length = n_trans, mode = "list") 
for (i in 1:length(wei_fits_cr)){
  wei_fits_cr[[i]] <- flexsurv::flexsurvreg(Surv(years, status) ~ factor(strategy_id), 
                                            data = mstate3_exdata$transitions, 
                                            subset = (trans == i) , 
                                            dist = "weibull") 
}
wei_fits_cr <- flexsurvreg_list(wei_fits_cr)

## -----------------------------------------------------------------------------
wei_fits_cf <- vector(length = n_trans, mode = "list") 
for (i in 1:length(wei_fits_cf)){
  wei_fits_cf[[i]] <- flexsurv::flexsurvreg(Surv(Tstart, Tstop, status) ~ factor(strategy_id), 
                                            data = mstate3_exdata$transitions, 
                                            subset = (trans == i) , 
                                            dist = "weibull") 
}
wei_fits_cf <- flexsurvreg_list(wei_fits_cf)

## -----------------------------------------------------------------------------
utility_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = mstate3_exdata$utility$mean,
                                       se = mstate3_exdata$utility$se),
                            dist = "beta")
head(utility_tbl)


## -----------------------------------------------------------------------------
drugcost_tbl <- stateval_tbl(data.table(strategy_id = strategies$strategy_id,
                                       est = mstate3_exdata$costs$drugs$costs),
                            dist = "fixed")
medcost_tbl <- stateval_tbl(data.table(state_id = states$state_id,
                                       mean = mstate3_exdata$costs$medical$mean,
                                       se = mstate3_exdata$costs$medical$se),
                            dist = "gamma")

## -----------------------------------------------------------------------------
n_samples <- 1000

## -----------------------------------------------------------------------------
transmod_data <- expand(hesim_dat, 
                        by = c("strategies", "patients"))
head(transmod_data)

## -----------------------------------------------------------------------------
transmod_cr <- create_IndivCtstmTrans(wei_fits_cr, transmod_data,
                                      trans_mat = tmat, n = n_samples,
                                      clock = "reset",
                                      start_age = patients$age)
transmod_cf <- create_IndivCtstmTrans(wei_fits_cf, transmod_data,
                                      trans_mat = tmat, n = n_samples,
                                      clock = "forward",
                                      start_age = patients$age)

## -----------------------------------------------------------------------------
# Predict hazard
transmod_data_pat1 <- transmod_data[patient_id == 1]
predict_haz <- function(fits, clock){
 transmod_cr_pat1 <- create_IndivCtstmTrans(fits, transmod_data_pat1,
                                            trans_mat = tmat, 
                                            clock = clock,
                                            uncertainty = "none")
  haz <- transmod_cr_pat1$hazard(t = seq(0, 20, 1))
  title_clock <- paste(toupper(substr(clock, 1, 1)), 
                       substr(clock, 2, nchar(clock)), sep="")
  haz[, clock := title_clock]
  return(haz[, ])
}

## ----predicted_hazards_plot, fig.width = 7, fig.height = 4--------------------
# Plot hazards
library("ggplot2")
haz <- rbind(predict_haz(wei_fits_cr, "reset"),
             predict_haz(wei_fits_cf, "forward"))
set_labels(haz, labels = labs, 
           new_names = c("strategy_name", "trans_name"))
ggplot(haz[t > 0], 
       aes(x = t, y = hazard, col = clock, linetype = strategy_name)) + 
  geom_line() + 
  facet_wrap(~trans_name) +
  xlab("Years") + ylab("Hazard") +
  scale_linetype_discrete(name = "Strategy") +
  scale_color_discrete(name = "Clock") 

## -----------------------------------------------------------------------------
# Utility
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

# Costs
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)
medcostmod <- create_StateVals(medcost_tbl, n = n_samples, hesim_data = hesim_dat)
costmods <- list(Drug = drugcostmod,
                 Medical = medcostmod)

## -----------------------------------------------------------------------------
econmod_cr <- IndivCtstm$new(trans_model = transmod_cr,
                             utility_model = utilitymod,
                             cost_models = costmods)
econmod_cf <- IndivCtstm$new(trans_model = transmod_cf,
                             utility_model = utilitymod,
                             cost_models = costmods)

## ---- cache = FALSE-----------------------------------------------------------
# "Clock reset"
econmod_cr$sim_disease()
head(econmod_cr$disprog_)

# "Clock forward"
econmod_cf$sim_disease()

## -----------------------------------------------------------------------------
econmod_cr$sim_stateprobs(t = seq(0, 20 , 1/12)) 

## ----stprobs_by_strategy_plot, fig.width = 7, fig.height = 4------------------
# Short function to create state probability "dataset" for plotting
summarize_stprobs <- function(stateprobs){
  x <- stateprobs[, .(prob_mean = mean(prob)),
                  by = c("strategy_id", "state_id", "t")]
  set_labels(x, labels = labs, new_names = c("strategy_name", "state_name"))
}

# Plot of state probabilities
stprobs_cr <- summarize_stprobs(econmod_cr$stateprobs_)
ggplot(stprobs_cr, aes(x = t, y = prob_mean, col = strategy_name)) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "Strategy") +
  theme(legend.position = "bottom") 

## ----stprobs_by_timescale_plot, fig.width = 7, fig.height = 4-----------------
econmod_cf$sim_stateprobs(t = seq(0, 20 , 1/12)) 
stprobs_cf <- summarize_stprobs(econmod_cf$stateprobs_)

# Compare "clock forward" and "clock reset" cases
stprobs <- rbind(data.table(stprobs_cf, clock = "Forward"),
                 data.table(stprobs_cr, clock = "Reset"))
ggplot(stprobs[strategy_id == 1], 
       aes(x = t, y = prob_mean, col = clock)) +
  geom_line() + facet_wrap(~state_name) + 
  xlab("Years") + ylab("Probability in health state") +
  scale_color_discrete(name = "Clock") +
  theme(legend.position = "bottom") 

## -----------------------------------------------------------------------------
econmod_cr$sim_qalys(dr = c(0,.03))
head(econmod_cr$qalys_)

## ----qalys_plot, fig.width = 6, fig.height = 4--------------------------------
qalys_summary <- econmod_cr$qalys_[, .(mean = mean(qalys)),
                                    by = c("strategy_id", "state_id", "dr")]
set_labels(qalys_summary, labels = labs,
           new_names = c("strategy_name", "state_name"))
ggplot(qalys_summary[dr == .03],
       aes(x = strategy_name, y = mean, fill = state_name)) + 
  geom_bar(stat = "identity") +
  scale_fill_discrete(name = "") +
  xlab("Strategy") + ylab("Mean QALYs") 

## -----------------------------------------------------------------------------
econmod_cr$sim_costs(dr = 0.03)
head(econmod_cr$costs_)

## -----------------------------------------------------------------------------
ce_sim <- econmod_cr$summarize()
format(summary(ce_sim, labels = labs))
cea_out <- cea(ce_sim, dr_qalys = .03, dr_costs = .03)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)

