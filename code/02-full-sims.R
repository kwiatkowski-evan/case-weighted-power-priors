rm(list = ls())

# for (idx_scale in seq(624, 854, by = 23)){

ptm <- proc.time()

if (Sys.getenv("USER") %in% c("kwiatkoe", "evankwiatkowski")) {
  setwd("~/Documents/GitHub/case-adaptive-pp/large-sample-share/code")
  Packages <- c("ggplot2", "survival", "ks", "survminer", "MASS", "hesim", "rjags", "psborrow", "cmdstanr","pracma")
  lapply(Packages, library, character.only = TRUE)
  idx_scale <- 601
} else { # longleaf
  # install.packages("survival", lib = "../rpkgs", repos='http://cran.us.r-project.org')
  library(ks, lib.loc = "../rpkgs/")
  library(survival, lib.loc = "../rpkgs/")
  library(psborrow, lib.loc = "../rpkgs/")
  library(deSolve, lib.loc = "../rpkgs/")
  library(mstate, lib.loc = "../rpkgs/")
  library(muhaz, lib.loc = "../rpkgs/")
  library(flexsurv, lib.loc = "../rpkgs/")
  library(hesim, lib.loc = "../rpkgs/")
  library(MASS)
  library(dplyr)
  library(rjags)
  library(cmdstanr)
  library(pracma)
  args <- commandArgs(trailingOnly = TRUE)  # sequence from batch file
  idx_scale  <- as.numeric(args[1]);
}

source("01-config.R")

X3beta <- c(rep(0, trt_n + ctrl_n),
            (X3_vec[idx_scale] == 0) * rep(c(log(c(1/((1:8) * 2), (1:8) * 2)), rep(0, length = 34)), each = 2) +
              (X3_vec[idx_scale] %in% c(1, 2)) * rep(1, h_n)) * scale_list[idx_scale]

outer_table <- matrix(NA, nrow = length(scale_list), ncol = length(names))
colnames(outer_table) <- names
outermost <- matrix(NA, nrow = idx_samp, ncol = length(names))
colnames(outermost) <- names
outer_a0 <- matrix(NA, nrow = (n.seg - 1) * h_n, ncol = idx_samp)

formula2 <- formula(paste("nu ~ treat +", paste(cov.names, collapse=" + "), " + interval + offset(log(expo)) - 1"))
formula1 <- formula(paste("nu ~ treat +", paste(cov.names, collapse=" + "), " + offset(log(expo))"))
formula2_no_trt <- formula(paste("nu ~ ", paste(cov.names, collapse=" + "), " + interval + offset(log(expo)) - 1"))
formula1_no_trt <- formula(paste("nu ~ ", paste(cov.names, collapse=" + "), " + offset(log(expo))"))

if (T1E_vec[idx_scale] == 1){ # (if T1E = 1 then treatment effect set to zero)
  mu2["treat"] <- 0
}

# Reformat data to create rate matrices used to generate piecewise exponential data. ratemat for events, ratemat2 for censoring
X.long <- X.all[rep(seq_len(nrow(X.all)), each = length(cut.time)), ]
for (j in 1:length(cut.time)){
  X.long[, paste0("interval", j)] <- rep(c(rep(0, j - 1), 1, rep(0, length(cut.time) - j)), nrow(X.all))
}
X.long <- as.matrix(X.long)
ratemat <- matrix(X.long[, c("treat", cov.names, interval.names)] %*% 
                    mu2[c("treat", cov.names, interval.names)], 
                  ncol = length(cut.time), byrow = T)
ratemat2 <- matrix(rep(theta2, nrow(X.all)), ncol = length(cut.time), byrow = T)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
for(idx_outer in 1:idx_samp){ # number of simulations (i.e. RCT and external data generated)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  set.seed(23 * as.integer(mod(idx_scale, 800) * idx_samp + idx_outer))  #  05-19-2020
  trt.idx  <- sample(which(X.all$treat == 1), trt_n)
  ctrl.idx <- sample(which(X.all$treat == 0 & X.all$external == 0), ctrl_n)
  # h.idx    <- sample(which(X.all$external == 1), h_n)
  # h.ecog.idx <- is.na(X.all[h.idx, "ecog" ]) # length h.idx
  h.idx <- 1:h_n
  h.ecog.idx <- h.idx
  
  t.sim <- 
    rpwexp(n = trt_n + ctrl_n + h_n,
           time = cut.time,
           rate = as.matrix(
             exp(ratemat[c(trt.idx, ctrl.idx, h.idx), ] + 
                   (X3_vec[idx_scale] %in% c(0, 1)) * (matrix(c(rep(0, trt_n + ctrl_n), rep(1, h_n)) * X3beta, ncol = 1) %*% matrix(c(1, 1), ncol = length(cut.time), nrow = 1)) + 
                   (X3_vec[idx_scale] %in% c(2)) * (matrix(c(rep(0, trt_n + ctrl_n), rep(1, h_n)) * X3beta, ncol = 1) %*% matrix(c(0, 1), ncol = length(cut.time), nrow = 1))
             ))) 
  # MAJOR ERROR 9/6/21
  # ALWAYS CHECK 9/14/21
  
  
  if (cen_vec[idx_scale] == 0){ # low censoring
    # y.sim <- pmin(t.sim, tail(sort(t.sim))[1])
    c.sim <- rpwexp(n = trt_n + ctrl_n + h_n,
                    rate = as.matrix(exp(ratemat2[c(trt.idx, ctrl.idx, h.idx), ] * 1.4)),  # low censoring
                    time = cut.time)
    y.sim <- c(t.sim[1:(trt_n + ctrl_n)], pmin(t.sim, c.sim)[(trt_n + ctrl_n + 1):(trt_n + ctrl_n + h_n)])
    nu.sim <- as.numeric(y.sim == t.sim)
    
    if (sum(nu.sim) == (trt_n + ctrl_n + h_n)){ ## redo if no external controls are censored
      c.sim <- tail(sort(t.sim[(trt_n + ctrl_n + 1):(trt_n + ctrl_n + h_n)]))[1]
      y.sim <- c(t.sim[1:(trt_n + ctrl_n)], pmin(t.sim, c.sim)[(trt_n + ctrl_n + 1):(trt_n + ctrl_n + h_n)])
      nu.sim <- as.numeric(y.sim == t.sim) # equals 1 if observed 
    }
  }
  if (cen_vec[idx_scale] == 1){ # high censoring 
    c.sim <- rpwexp(n = trt_n + ctrl_n + h_n,
                    rate = as.matrix(exp(ratemat2[c(trt.idx, ctrl.idx, h.idx), ] * 0.9)), # high censoring 
                    time = cut.time)
    y.sim <- c(t.sim[1:(trt_n + ctrl_n)], pmin(t.sim, c.sim)[(trt_n + ctrl_n + 1):(trt_n + ctrl_n + h_n)])
    nu.sim <- as.numeric(y.sim == t.sim)
  }
  
  # Simulated design matrix
  sim.data           <- data.frame(rbind(cbind(t = y.sim, nu = nu.sim,
                                               X.all[c(trt.idx, ctrl.idx, h.idx), 
                                                     c("treat", "external", cov.names, "ecog")])))
  # Add patient ID to sim.data
  sim.data$id <- seq.int(nrow(sim.data))
  # Reformat data using survSplit()
  sim.split          <- survSplit(sim.data, cut = cut.time, end = "t", start = "t0", 
                                  event = "nu", episode = "interval")
  sim.split$interval <- factor(sim.split$interval - 1)
  sim.split$expo     <- sim.split$t - sim.split$t0
  
  ########################################################################################################################################################################################################
  # Add iptw
  ps_formula <- formula(paste("1 - external ~", paste(cov.names, collapse=" + ")))
  ps_fits <- glm(ps_formula,
                 family = binomial(),
                 data = sim.data)
  ps <- predict(ps_fits, type = "response")
  iptw_att <- ps / (1 - ps)
  sim.data <- sim.data %>% dplyr::mutate(iptw_att = ifelse(external == 1, iptw_att, 1))
  
  ########################################################################################################################################################################################################
  
  # Fit model for historical events
  if (n.seg >= 3){
    sim.poisson <- glm(formula2_no_trt, 
                       family = poisson(link = "log"), 
                       data = sim.split[sim.split$external == 1, ])
  } else {
    sim.poisson <- glm(formula1_no_trt, 
                       family = poisson(link = "log"), 
                       data = sim.split[sim.split$external == 1, ])
    names(sim.poisson$coefficients)[names(sim.poisson$coefficients) == "(Intercept)"] <- "interval1"
  }
  lambda <- mvrnorm(1, sim.poisson$coefficients[c(cov.names, interval.names)], 
                    vcov(sim.poisson)[c(cov.names, interval.names), 
                                      c(cov.names, interval.names)])
  # Fit model for historical censoring
  if (n.seg >= 3){
    sim.poisson <- glm(1 - nu ~ interval + offset(log(expo)) - 1,
                       family = poisson(link = "log"),
                       data = sim.split[sim.split$external == 1, ])
  } else {
    sim.poisson <- glm(1 - nu ~ offset(log(expo)),
                       family = poisson(link = "log"),
                       sim.split[sim.split$external == 1, ])
    names(sim.poisson$coefficients)[names(sim.poisson$coefficients) == "(Intercept)"] <- "interval1"
  }
  theta <- mvrnorm(1, sim.poisson$coefficients[interval.names],
                   vcov(sim.poisson)[interval.names, interval.names])
  
  # Fit model for RCT events
  if (n.seg >= 3){
    sim.poisson.lambda <- glm(formula2,
                              family = poisson(link = "log"), 
                              data = sim.split[sim.split$external == 0, ])
  } else {
    sim.poisson.lambda <- glm(formula1,
                              family = poisson(link = "log"), 
                              data = sim.split[sim.split$external == 0, ])
    names(sim.poisson.lambda$coefficients)[names(sim.poisson.lambda$coefficients) == "(Intercept)"] <- "interval1"
  }
  # Fit model for historical censoring
  if (n.seg >= 3){
    sim.poisson.theta   <- glm(1 - nu ~ interval + offset(log(expo)) - 1,
                               family = poisson(link = "log"),
                               data = sim.split[sim.split$external == 1, ])
  } else {
    sim.poisson.theta   <- glm(1 - nu ~ offset(log(expo)),
                               family = poisson(link = "log"),
                               data = sim.split[sim.split$external == 1, ])
    names(sim.poisson.theta$coefficients)[names(sim.poisson.theta$coefficients) == "(Intercept)"] <- "interval1"
  }
  
  ########################################################################################################################################################################################################
  
  innermost  <- matrix(NA, nrow = idx_n, ncol = length(names))
  colnames(innermost) <- names
  a0.outer <- matrix(NA, nrow = nrow(sim.split[sim.split$external == 1, ]), ncol = idx_n)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  for(idx in 1:idx_n){ # repetitions per design (i.e. number of imputations per exposure time interval in external controls)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    
    sim.split.idx <- sim.split # 8/3/21
    
    # Add imputation for multiple intervals in historical controls
    if (n.seg >= 3){ 
      for (i in which(sim.split.idx$external == 1)){
        for(j in 1:length(inner.t)){
          if (sim.split.idx$nu[i] == 0 & sim.split.idx$interval[i] == j & sim.split.idx$expo[i] == diff[j]){
            temp.t <- rexp(1, rate = exp(as.matrix(sim.split.idx[i, c(cov.names)]) %*%
                                           as.matrix(lambda[c(cov.names)]) + 
                                           lambda[paste0("interval", j)]))
            temp.c  <- rexp(1, rate = exp(theta[paste0("interval", j)]))
            temp.y  <- pmin(temp.t, temp.c)
            temp.nu <- as.numeric(temp.y == temp.t)
            sim.split.idx$expo[i] <- diff[j] + temp.y
            sim.split.idx$nu[i]   <- temp.nu
          }
        }
      }
    }
    
    # subset to historical only
    sim.split.h <- sim.split.idx[sim.split.idx$external == 1, ]
    
    box.p.log <- matrix(NA, ncol = idx_a0_n, nrow = nrow(sim.split.h))
    # repetitions per a0 (i.e. number of compatibility scores computed per exposure time interval)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    for (idx_a0 in 1:idx_a0_n){
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
      lambda <- mvrnorm(1, sim.poisson.lambda$coefficients[c(cov.names, interval.names)], 
                        vcov(sim.poisson.lambda)[c(cov.names, interval.names), 
                                                 c(cov.names, interval.names)])
      theta <- mvrnorm(1, sim.poisson.theta$coefficients[interval.names],
                       vcov(sim.poisson.theta)[interval.names, interval.names])
      # debugging code - prevent theta from being positive
      for (i in 1:length(theta)){
        if (theta[i] > 0) theta[i] <- min(theta)
      }
      
      # generate y.h.pred incorporating variability in estimating lambda and assess compatibility, use interval + expo
      for(i in 1:nrow(sim.split.h)){
        j <- c(sim.split.h[i, "interval"])
        y.h.pred <- rexp(1E4, rate = exp(as.matrix(sim.split.h[i, c(cov.names)]) %*%
                                           as.matrix(lambda[c(cov.names)]) +
                                           lambda[paste0("interval", j)]) +
                           exp(theta[paste0("interval", j)]))
        fhat_outer         <- kde(log(y.h.pred))
        fhat_predict_outer <- predict(fhat_outer, x = log(y.h.pred))
        box.p.log[i, idx_a0] <- mean(fhat_predict_outer <= predict(fhat_outer, x = log(sim.split.h$expo[i])))
      }
      # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    } # end for (idx_a0 in 1:idx_a0_n)
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
    ## Compute case-specific weights
    a0 <- rowMeans(box.p.log)
    a0.outer[, idx] <- a0
    sim.split.h$a0 <- a0
    # Compute average weights a0 overall and by ecog score
    innermost[idx, "mean_a0"] <- mean(a0)
    innermost[idx, c("a0_mean_q0.1", "a0_mean_q0.25", "a0_mean_q0.5", "a0_mean_q0.75", "a0_mean_q0.9")] <- quantile(a0, probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    innermost[idx, "mean_a0_r"] <- mean(rotate3(a0, p = 1, q = 0))
    innermost[idx, c("a0_mean_r_q0.1", "a0_mean_r_q0.25", "a0_mean_r_q0.5", "a0_mean_r_q0.75", "a0_mean_r_q0.9")] <- quantile(rotate3(a0, p = 1, q = 0), probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    innermost[idx, "a0_wt_ecog0"] <- mean(sim.split.h[sim.split.h[,"ecog"] %in% 0, "a0"])
    innermost[idx, "a0_wt_ecog1"] <- mean(sim.split.h[sim.split.h[,"ecog"] %in% 1, "a0"])
    innermost[idx, "a0_wt_ecogNA"] <- mean(sim.split.h[is.na(sim.split.h[,"ecog"]), "a0"])
    innermost[idx, "a0_wt_obs"] <- mean(sim.split.h[sim.split.h$nu == 1, "a0"])
    innermost[idx, "a0_wt_cen"] <- mean(sim.split.h[sim.split.h$nu == 0, "a0"])
    # Compute average weights a0 by interval, and by ecog and interval combination
    for (i1 in c(0, 1)){
      for (j1 in 1:(n.seg - 1)){
        innermost[idx, paste0("a0_wts_", i1, "_", j1)] <- 
          mean(sim.split.h[sim.split.h[,"ecog"] %in% i1 & sim.split.h[,"interval"] == j1, "a0"])
      }
    }
    for (j1 in 1:(n.seg - 1)){
      innermost[idx, paste0("a0_wts_", 2, "_", j1)] <-  mean(sim.split.h[is.na(sim.split.h[,"ecog"]) & sim.split.h[,"interval"] == j1, "a0"])
      innermost[idx, paste0("a0_wt_int", j1)] <- mean(sim.split.h[sim.split.h[,"interval"] == j1, "a0"])
      innermost[idx, paste0("a0_wt_int", j1, c("_q0.1", "_q0.25", "_q0.5", "_q0.75", "_q0.9"))] <- quantile(sim.split.h[sim.split.h[,"interval"] == j1, "a0"], probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
    }
    
    # Final analysis for adaptive weights a0_v1
    wts <- c(rep(1, times = sum(sim.split.idx$external == 0)), a0)
    if (n.seg >= 3){
      sim.poisson <- glm(formula2,
                         family = poisson(link = "log"),
                         data = sim.split.idx,
                         weights = wts)
    } else {
      sim.poisson <- glm(formula1,
                         family = poisson(link = "log"),
                         data = sim.split.idx,
                         weights = wts)
    }
    innermost[idx, paste0(c("HR_", "lower_", "upper_"), "a0_v1")] <-
      c(exp(coef(sim.poisson))["treat"],
        exp(coef(sim.poisson)["treat"])*exp(-1.96*summary(sim.poisson)$coefficients["treat", 2]),
        exp(coef(sim.poisson)["treat"])*exp(1.96*summary(sim.poisson)$coefficients["treat", 2]))
    
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  } # ends for(idx in 1:idx_n) # repetitions per design (i.e. number of imputations per exposure time interval in external controls)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  #  For each idx_scale and idx_sample, summarize data from idx_n repeats
  innermost <- data.frame(innermost)
  innermost$bias_a0 <- (innermost$HR_a0 - exp(mu2["treat"])) # BIAS
  # Compute MSE, HR, credible interval width, and power for "a0_v1"
  innermost[[paste0("MSE_", "a0_v1")]]      <- 
    (innermost[[paste0("HR_", "a0_v1")]] - exp(mu2["treat"])) ^ 2
  innermost[[paste0("CI_width_", "a0_v1")]] <- 
    (innermost[[paste0("upper_", "a0_v1")]] - innermost[[paste0("lower_", "a0_v1")]])
  innermost[[paste0("CI_", "a0_v1")]]       <- 
    (innermost[[paste0("lower_", "a0_v1")]] < exp(mu2["treat"]) & 
       exp(mu2["treat"]) < innermost[[paste0("upper_", "a0_v1")]])
  innermost[[paste0("power_", "a0_v1")]]    <- 
    ((innermost[[paste0("lower_", "a0_v1")]] < 1 & innermost[[paste0("upper_", "a0_v1")]] < 1) 
     # | (innermost[[paste0("lower_", "a0_v1")]] > 1 & innermost[[paste0("upper_", "a0_v1")]] > 1)
    )
  
  outermost[idx_outer, ] <- colMeans(innermost)
  
  # Final analysis using fixed power prior weights & adaptive weights
  for (i in 1:length(ppList)){
    if (ppList[i] == "a0"){
      wts <- c(rep(1, times = sum(sim.split$external == 0)), rowMeans(a0.outer))
    } else if (substring(ppList[i], 1, 1) == "p"){
      wts <- c(rep(1, times = sum(sim.split$external == 0)), 
               rotate3(x = rowMeans(a0.outer), 
                       p = as.numeric(substring(ppList[i], 2, 2)), 
                       q = as.numeric(substring(ppList[i], 5, 5))) *
                 g(x = 2 * mean(a0.outer), 
                   q = as.numeric(substring(ppList[i], 9, 13))))
    # } else if (substring(ppList[i], 1, 1) == "q"){
    #   wts <- c(rep(1, times = sum(sim.split$external == 0)), 
    #            rotate2(rowMeans(a0.outer), as.numeric(substring(ppList[i], 2, 5)), 0) * 
    #              g(2 * mean(a0.outer), as.numeric(substring(ppList[i], 9, 13))))
    # } else if (substring(ppList[i], 1, 1) == "e"){
    #   wts <- c(rep(1, times = sum(sim.split$external == 0)), 
    #            rotate(rowMeans(a0.outer), 0, as.numeric(substring(ppList[i], 2, 5))) * 
    #              g(2 * mean(a0.outer), as.numeric(substring(ppList[i], 9, 13))))
    } else {
      wts <- c(rep(1, times = sum(sim.split$external == 0)), rep(as.numeric(ppList[i]), sum(sim.split$external == 1)))
    }
    if (n.seg >= 3){
      sim.poisson <- glm(formula2,
                         family = poisson(link = "log"),
                         data = sim.split,
                         weights = wts)
    } else {
      sim.poisson <- glm(formula1,
                         family = poisson(link = "log"),
                         data = sim.split,
                         weights = wts)
    }
    outermost[idx_outer, paste0(c("HR_", "lower_", "upper_"), ppList[i])] <-
      c(exp(coef(sim.poisson))["treat"],
        exp(coef(sim.poisson)["treat"])*exp(-1.96*summary(sim.poisson)$coefficients["treat", 2]),
        exp(coef(sim.poisson)["treat"])*exp(1.96*summary(sim.poisson)$coefficients["treat", 2]))
    if (ppList[i] == "a0"){
      outermost[idx_outer, "a0_logHR"] <- coef(sim.poisson)["treat"]
      outermost[idx_outer, "SE_a0"] <- summary(sim.poisson)$coefficients["treat", 2]
    }
  }
  
  # # Final analysis using frailty model
  # f_formula <- formula(paste("Surv(t, nu) ~ treat + frailty(external) +",
  #                            paste(cov.names, collapse=" + ")))
  # coxf <- tryCatch(coxph(f_formula, data = sim.data),
  #                  error = function(e) coxph(f_formula, data = sim.data))
  # outermost[idx_outer, c("HR_f", "lower_f", "upper_f", paste0(cov.names, "_f"))] <-
  #   c(exp(coef(coxf))["treat"],
  #     exp(coef(coxf)["treat"])*exp(-1.96*summary(coxf)$coefficients["treat", "se(coef)"]),
  #     exp(coef(coxf)["treat"])*exp(1.96*summary(coxf)$coefficients["treat", "se(coef)"]),
  #     summary(coxf)$coefficients[cov.names, "coef"])
  
  # Final analysis using cox model
  cox_formula <- formula(paste("Surv(t, nu) ~ treat + ", paste(cov.names, collapse=" + ")))
  cox <- coxph(cox_formula, data = sim.data[sim.data$external == 0, ])
  outermost[idx_outer, c("HR_cox", "lower_cox", "upper_cox", paste0(cov.names, "_cox"))] <-
    c(exp(coef(cox))["treat"],
      exp(coef(cox)["treat"])*exp(-1.96*summary(cox)$coefficients["treat", "se(coef)"]),
      exp(coef(cox)["treat"])*exp(1.96*summary(cox)$coefficients["treat", "se(coef)"]),
      summary(cox)$coefficients[cov.names, "coef"])
  
  # Final analysis using iptw
  iptw_formula <-  formula(paste("Surv(t, nu) ~ treat + cluster(id) + ", paste(cov.names, collapse=" + ")))
  iptw <- coxph(iptw_formula, 
                data = sim.data,
                weights = iptw_att)
  outermost[idx_outer, c("HR_iptw", "lower_iptw", "upper_iptw", paste0(cov.names, "_iptw"))] <-
    c(exp(coef(iptw))["treat"],
      exp(coef(iptw)["treat"])*exp(-1.96*summary(iptw)$coefficients["treat", "robust se"]),
      exp(coef(iptw)["treat"])*exp(1.96*summary(iptw)$coefficients["treat", "robust se"]),
      summary(iptw)$coefficients[cov.names, "coef"])
  
  # Final analysis using commensurate priors
  if (run_com == TRUE){ 
    
    survival.dat <- list(
      N = nrow(sim.split),
      nu = sim.split$nu,
      trt = sim.split$treat,
      X1 = sim.split$age,
      X2 = sim.split$sex2,
      a0 = c(rep(1, times = sum(sim.split$external == 0)), a0), # 8/28/21
      offset = log(sim.split$expo),
      int_d_1 = as.numeric(sim.split$interval == 1),
      int_d_2 = as.numeric(sim.split$interval == 2),
      E = as.numeric(sim.split$external == 0)
    )
    
    skip_to_next <- FALSE
    file <- "com_cauchy.stan"
    mod <- cmdstan_model(file, pedantic = TRUE)
    tryCatch(fit_com_cauchy <- mod$sample(
      data = survival.dat,
      seed = as.integer(idx_scale*idx_samp + idx_outer),
      chains = 2,
      parallel_chains = 2,
      refresh = 0,
      iter_warmup = 500,
      iter_sampling = 2500
    ), error = function(e) {skip_to_next <<- TRUE})
    
    if (skip_to_next == FALSE){
      res_com_cauchy <- fit_com_cauchy$summary("gamma", "mean", sig = ~ mean(. <= 0))
      outermost[idx_outer, c("HR_com_cauchy")] <- as.numeric(exp(res_com_cauchy["mean"]))
      outermost[idx_outer, c("power_com_cauchy")] <- as.numeric(res_com_cauchy["sig"] > 0.975)
      outermost[idx_outer, c("MSE_com_cauchy")] <- (outermost[idx_outer, c("HR_com_cauchy")] - exp(mu2["treat"])) ^ 2
    }
    
    skip_to_next <- FALSE
    file <- "com_a0.stan"
    mod <- cmdstan_model(file, pedantic = TRUE)
    tryCatch(fit_com_a0 <- mod$sample(
      data = survival.dat,
      seed = as.integer(idx_scale*idx_samp + idx_outer),
      chains = 2,
      parallel_chains = 2,
      refresh = 0,
      iter_warmup = 500,
      iter_sampling = 2500
    ), error = function(e) {skip_to_next <<- TRUE})
    
    if (skip_to_next == FALSE){
      res_com_a0 <- fit_com_a0$summary("gamma", "mean", sig = ~ mean(. <= 0))
      outermost[idx_outer, c("HR_com_a0")] <- as.numeric(exp(res_com_a0["mean"]))
      outermost[idx_outer, c("power_com_a0")] <- as.numeric(res_com_a0["sig"] > 0.975)
      outermost[idx_outer, c("MSE_com_a0")] <- (outermost[idx_outer, c("HR_com_a0")] - exp(mu2["treat"])) ^ 2
    }
    
  } # end commensurate prior
  
  ########################################################################################################################################################################################################
  
  # Compute MSE, HR, credible interval width, and power for each analysis method
  for (i in 1:length(nList)){
    outermost[, paste0("MSE_", nList[i])]      <- 
      (outermost[, paste0("HR_", nList[i])] - exp(mu2["treat"])) ^ 2
    outermost[, paste0("CI_width_", nList[i])] <- 
      (outermost[, paste0("upper_", nList[i])] - outermost[, paste0("lower_", nList[i])])
    outermost[, paste0("CI_", nList[i])]       <- 
      (outermost[, paste0("lower_", nList[i])] < exp(mu2["treat"]) & 
         exp(mu2["treat"]) < outermost[, paste0("upper_", nList[i])])
    outermost[, paste0("power_", nList[i])]    <- 
      ((outermost[, paste0("lower_", nList[i])] < 1 & outermost[, paste0("upper_", nList[i])] < 1) 
       # | (outermost[, paste0("lower_", nList[i])] > 1 & outermost[, paste0("upper_", nList[i])] > 1)
      )
  }
  
  outermost[idx_outer, "nu_trt"]  <- sum(sim.data$nu[sim.data$treat == 1 & sim.data$external == 0])
  outermost[idx_outer, "nu_ctrl"] <- sum(sim.data$nu[sim.data$treat == 0 & sim.data$external == 0])
  outermost[idx_outer, "nu_h"]    <- sum(sim.data$nu[sim.data$treat == 0 & sim.data$external == 1])
  outermost[idx_outer, "T1E"]     <- T1E_vec[idx_scale] 
  
  outer_a0[1:length(a0), idx_outer] <- a0
  print(idx_outer)
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
} # ends for(idx_outer in 1:idx_samp)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

outermost[, "scale"]   <- scale_list[idx_scale]  # scale
if (save_a0 == TRUE){ 
  write.csv(outer_a0, paste0("../output - a0/",  formatC(idx_scale, width=4, flag="0"), ".csv"))
  # write.csv(outer_a0, paste0("../output-a0/",  formatC(idx_scale, width=4, flag="0"), ".csv"))
}
write.csv(c(colMeans(outermost)), paste0("../output/",  formatC(idx_scale, width=4, flag="0"), ".csv"))

proc.time() - ptm

# }