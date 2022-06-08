## ---- include = FALSE---------------------------------------------------------
ggplot2::theme_set(ggplot2::theme_bw())

## ---- warning = FALSE, message = FALSE----------------------------------------
library("data.table")
library("hesim")
set.seed(101)
strategies <- data.table(strategy_id = c(1, 2, 3),
                         strategy_name = paste0("Strategy ", 1:3))
patients <- data.table(patient_id = seq(1, 3),
                       age = c(45, 50, 60),
                       female = c(0, 0, 1))
states <- data.frame(state_id =  seq(1, 3),
                     state_name = c("Progression free (1st line)",
                                    "Progression free (2nd line)",
                                    "Progressed (2nd line)"),
                        stringsAsFactors = FALSE)
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
print(hesim_dat)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
print(labs)

## ---- warning = FALSE, message = FALSE----------------------------------------
library("flexsurv")
surv_est_data <- psm4_exdata$survival
fit1 <- flexsurv::flexsurvreg(Surv(endpoint1_time, endpoint1_status) ~ age + female + 
                                   factor(strategy_id),
                              data = surv_est_data,
                              dist = "weibull")
fit2 <- flexsurv::flexsurvreg(Surv(endpoint2_time, endpoint2_status) ~ age + female +
                                   factor(strategy_id),
                              data = surv_est_data,
                              dist = "weibull")
fit3 <- flexsurv::flexsurvreg(Surv(endpoint3_time, endpoint3_status) ~ age + female +
                                   factor(strategy_id),
                              data = surv_est_data,
                              dist = "weibull")
psfit_wei <- flexsurvreg_list(fit1, fit2, fit3)

## ---- warning = FALSE, message = FALSE----------------------------------------
utility_tbl <- stateval_tbl(tbl = data.frame(state_id = states$state_id,
                                             min = psm4_exdata$utility$lower,
                                             max = psm4_exdata$utility$upper),
                            dist = "unif")

## ---- warning = FALSE, message = FALSE----------------------------------------
drugcost_tbl <- stateval_tbl(tbl = data.frame(strategy_id = strategies$strategy_id,
                                           est = psm4_exdata$costs$drugs$costs),
                          dist = "fixed")

## ---- warning = FALSE, message = FALSE----------------------------------------
medcosts_fit <- lm(costs ~ female + state_name, data = psm4_exdata$costs$medical)

## ---- warning = FALSE, message = FALSE----------------------------------------
n_samples <- 100

## ---- warning = FALSE, message = FALSE----------------------------------------
surv_input_data <- expand(hesim_dat, by = c("strategies", "patients"))
survmods <- create_PsmCurves(psfit_wei, input_data = surv_input_data, n = n_samples,
                               uncertainty = "bootstrap", est_data = surv_est_data)

## ---- warning = FALSE, message = FALSE----------------------------------------
utilitymod <- create_StateVals(utility_tbl, n = n_samples, hesim_data = hesim_dat)

## ---- warning = FALSE, message = FALSE----------------------------------------
drugcostmod <- create_StateVals(drugcost_tbl, n = n_samples, hesim_data = hesim_dat)

## ---- warning = FALSE, message = FALSE----------------------------------------
medcost_data <- expand(hesim_dat, by = c("strategies", "patients", "states"))
medcostmod <- create_StateVals(medcosts_fit, input_data = medcost_data, 
                                 n = n_samples)

## ---- warning = FALSE, message = FALSE----------------------------------------
psm <- Psm$new(survival_models = survmods,
               utility_model = utilitymod,
               cost_models = list(medical = medcostmod, 
                                  drug = drugcostmod))

## ----survcurves, warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4----
# Simulate
times <- seq(0, 10, by = .1)
psm$sim_survival(t = times)
head(psm$survival_)

# Plot
labs$curve <- c("Curve 1" = 1, "Curve 2" = 2, "Curve 3" = 3)
autoplot(psm$survival_, labels = labs,
         ci = TRUE, ci_style = "ribbon")

## ----stateprobs, warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4----
# Simulate
psm$sim_stateprobs()
head(psm$stateprobs_)

# Plot
library("ggplot2")
stateprobs <- psm$stateprobs_[sample == 30 & patient_id == 1]
stateprobs[, state := factor(state_id, 
                             levels = rev(unique(state_id)))]
stateprobs[, strategy := factor(strategy_id, labels = names(labs$strategy_id))]
ggplot(stateprobs[strategy_id %in% c(1, 2)],
            aes(x = t, y = prob, fill = state, group = state)) + 
  geom_area(alpha = .65) + facet_wrap(~strategy) +
  xlab("Years") + ylab("Proportion in state") +
  scale_fill_manual(name = "Health state",
                    values = c("gray92", "green4", "orange", "purple"),
                    labels = c("Death", 
                               rev(hesim_dat$states$state_name))) +
  guides(fill = guide_legend(reverse = TRUE,
                             nrow = 2, byrow = TRUE)) + 
  theme(legend.position = "bottom")

## ---- warning = FALSE, message = FALSE----------------------------------------
# Costs
psm$sim_costs(dr = .03)
head(psm$costs_)

# QALYs
psm$sim_qalys(dr = .03)
head(psm$qalys_)

## -----------------------------------------------------------------------------
ce_sim <- psm$summarize()
cea_out <- cea(ce_sim, dr_qalys = .03, dr_costs = .03)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = .03, dr_costs = .03)

