## ---- include = FALSE---------------------------------------------------------
library("ggplot2")
theme_set(theme_bw())

## ---- out.width = "700px", echo = FALSE---------------------------------------
knitr::include_graphics("markov-cohort.png")

## ---- warning = FALSE, message = FALSE----------------------------------------
library("hesim")
library("data.table")
strategies <- data.table(strategy_id = 1:2,
                         strategy_name = c("Monotherapy", "Combination therapy"))
patients <- data.table(patient_id = 1)
states <- data.table(state_id = 1:3,
                     state_name = c("State A", "State B", "State C"))
hesim_dat <- hesim_data(strategies = strategies,
                        patients = patients,
                        states = states)
print(hesim_dat)

## -----------------------------------------------------------------------------
labs <- get_labels(hesim_dat)
print(labs)

## -----------------------------------------------------------------------------
trans_mono <- matrix(c(1251, 350, 116, 17,
                       0, 731, 512, 15,
                       0, 0, 1312, 437,
                       0, 0, 0, 469),
                      ncol = 4, nrow = 4, byrow = TRUE)
colnames(trans_mono) <- rownames(trans_mono) <-  c("A", "B", "C", "D")
print(trans_mono)

## -----------------------------------------------------------------------------
params <- list(
  alpha_mono = trans_mono, 
  lrr_mean = log(.509), 
  lrr_lower <- log(.365),
  lrr_upper = log(.710),
  c_dmed_mean = c(A = 1701, B = 1774, C = 6948),
  c_cmed_mean = c(A = 1055, B = 1278, C = 2059),
  c_zido = 2278,
  c_lam = 2086.50,
  u = 1
)

## -----------------------------------------------------------------------------
rng_def <- define_rng({
  lrr_se <- (lrr_upper - lrr_lower)/(2 * qnorm(.975)) # Local object 
                                                      # not returned
  list( # Parameters to return
    p_mono = dirichlet_rng(alpha_mono),
    rr_comb = lognormal_rng(lrr_mean, lrr_se),
    c_zido = c_zido,
    c_lam = c_lam,
    c_dmed = gamma_rng(mean = c_dmed_mean, sd = c_dmed_mean),
    c_cmed = gamma_rng(mean = c_cmed_mean, sd = c_cmed_mean),
    u = u
  )
}, n = 1000)

## -----------------------------------------------------------------------------
input_data <- expand(hesim_dat, by = c("strategies", "patients"))
head(input_data)

## -----------------------------------------------------------------------------
tparams_def <- define_tparams({
  ## The treatment effect (relative risk) is transformed so that it varies by 
  ## strategies and only applies for the first 2 years (Monotherapy is 
  ## the reference strategy). 
  rr <- ifelse(strategy_name == "Monotherapy" | time > 2, 1, rr_comb)
  
  list(
    tpmatrix = tpmatrix(
      C, p_mono$A_B * rr, p_mono$A_C * rr, p_mono$A_D * rr,
      0, C,               p_mono$B_C * rr, p_mono$B_D * rr,
      0, 0,               C,               p_mono$C_D * rr,
      0, 0,               0,               1
    ),
    utility = u,
    costs = list(
        drug = ifelse(strategy_name == "Monotherapy" | time > 2,
                      c_zido, c_zido + c_lam),
        community_medical = c_cmed,
        direct_medical = c_dmed
    )
  )
}, times = c(2, Inf))

## -----------------------------------------------------------------------------
mod_def <- define_model(tparams_def = tparams_def, 
                        rng_def = rng_def, 
                        params = params)

## ----econmod------------------------------------------------------------------
econmod <- create_CohortDtstm(mod_def, input_data)

## ----simStateprobs------------------------------------------------------------
econmod$sim_stateprobs(n_cycles = 20)
head(econmod$stateprobs_)

## ----simStateprobsPlot, fig.width = 7, fig.height = 6-------------------------
autoplot(econmod$stateprobs_, labels = labs,
         ci = TRUE, ci_style = "ribbon")

## ----simQALYs-----------------------------------------------------------------
econmod$sim_qalys(dr = 0, integrate_method = "riemann_right")
head(econmod$qalys_)

## ----simCosts-----------------------------------------------------------------
econmod$sim_costs(dr = 0.06, integrate_method = "riemann_right")
head(econmod$costs_)

## ----cea----------------------------------------------------------------------
ce_sim <- econmod$summarize()
wtp <- seq(0, 25000, 500)
cea_pw_out <- cea_pw(ce_sim, comparator = 1, dr_qalys = 0, dr_costs = .06,
                     k = wtp)

## ----icer---------------------------------------------------------------------
format(icer(cea_pw_out))

## ----ceac, fig.width = 6, fig.height = 4--------------------------------------
plot_ceac(cea_pw_out, labels = labs)

