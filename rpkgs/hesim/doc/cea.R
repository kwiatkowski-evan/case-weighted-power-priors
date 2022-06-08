## ---- include = FALSE---------------------------------------------------------
library("knitr")
opts_chunk$set(fig.width = 8, fig.height = 4.5)

## ----ce_output, warning = FALSE, message = FALSE------------------------------
set.seed(131)
n_samples <- 1000

# cost
c <- vector(mode = "list", length = 6)
names(c) <- c("Strategy 1, Grp 1", "Strategy 1, Grp 2", "Strategy 2, Grp 1",
              "Strategy 2, Grp 2", "Strategy 3, Grp 1", "Strategy 3, Grp 2")
c[[1]] <- rlnorm(n_samples, 2, .1)
c[[2]] <- rlnorm(n_samples, 2, .1)
c[[3]] <- rlnorm(n_samples, 11, .15)
c[[4]] <- rlnorm(n_samples, 11, .15)
c[[5]] <- rlnorm(n_samples, 11, .15)
c[[6]] <- rlnorm(n_samples, 11, .15)

# effectiveness
e <- c
e[[1]] <- rnorm(n_samples, 8, .2)
e[[2]] <- rnorm(n_samples, 8, .2)
e[[3]] <- rnorm(n_samples, 10, .8)
e[[4]] <- rnorm(n_samples, 10.5, .8)
e[[5]] <- rnorm(n_samples, 8.5, .6)
e[[6]] <- rnorm(n_samples, 11, .6)

# cost and effectiveness by strategy and simulation
library("data.table")
ce <- data.table(sample = rep(seq(n_samples), length(e)),
                 strategy = rep(paste0("Strategy ", seq(1, 3)), 
                                each = n_samples * 2),
                 grp = rep(rep(c("Group 1", "Group 2"),
                               each = n_samples), 3),
                 cost = do.call("c", c), qalys = do.call("c", e))
head(ce)

## ----cea, warning = FALSE, message = FALSE------------------------------------
library("hesim")
ktop <- 200000
cea_out <-  cea(ce, k = seq(0, ktop, 500), sample = "sample", strategy = "strategy",
                grp = "grp", e = "qalys", c = "cost")

## ----cea_pw-------------------------------------------------------------------
cea_pw_out <-  cea_pw(ce,  k = seq(0, ktop, 500), comparator = "Strategy 1",
                      sample = "sample", strategy = "strategy", grp = "grp",
                       e = "qalys", c = "cost")

## -----------------------------------------------------------------------------
library("magrittr") # Use pipes
icer(cea_pw_out, k = 50000) %>%
  format()

## -----------------------------------------------------------------------------
head(cea_pw_out$delta)

## ----ceplane_plot-------------------------------------------------------------
library("ggplot2")
theme_set(theme_bw())
plot_ceplane(cea_pw_out, k = 50000)

## ----mce_example_setup, echo = -1, warning = FALSE, message = FALSE-----------
library("knitr")
ce <- ce[, nmb := 50000 * qalys - cost]
random_rows <- sample(1:n_samples, 10)
nmb <- dcast(ce[sample %in% random_rows & grp == "Group 2"], 
                sample ~ strategy, value.var = "nmb")
setnames(nmb, colnames(nmb), c("sample", "nmb1", "nmb2", "nmb3"))
nmb <- nmb[, maxj := apply(nmb[, .(nmb1, nmb2, nmb3)], 1, which.max)]
nmb <- nmb[, maxj := factor(maxj, levels = c(1, 2, 3))]

## ----mce_example, echo = -1---------------------------------------------------
kable(nmb, digits = 0, format = "html")
mce <- prop.table(table(nmb$maxj))
print(mce)

## ----mce_plot, warning = FALSE, message = FALSE-------------------------------
plot_ceac(cea_out)

## ----ceac_plot----------------------------------------------------------------
plot_ceac(cea_pw_out)

## ----ceaf_plot----------------------------------------------------------------
plot_ceaf(cea_out)

## ----evpi_example_a-----------------------------------------------------------
# Expected net monetary benefit
enmb <- ce[, .(enmb = mean(nmb)), by = c("strategy", "grp")]
enmb <- dcast(enmb, strategy ~ grp, value.var = "enmb")
strategymax_g2 <- which.max(enmb[[3]]) # Optimal strategy group 2

# Net monetary benefit (with perfect and current information)
nmb <- nmb[, nmbpi := apply(nmb[, .(nmb1, nmb2, nmb3)], 1, max)]
nmb <- nmb[, nmbci := nmb[[strategymax_g2 + 1]]]
kable(nmb, digits = 0, format = "html")

## ----evpi_example_b-----------------------------------------------------------
enmbpi <- mean(nmb$nmbpi)
enmbci <- mean(nmb$nmbci)
print(enmbpi)
print(enmbci)
print(enmbpi - enmbci)

## ----evpi_plot----------------------------------------------------------------
plot_evpi(cea_out)

## ----totevpi, fig.width = 6---------------------------------------------------
w_dt <- data.table(grp = paste0("Group ", seq(1, 2)), w = c(0.25, .75))
evpi <- cea_out$evpi
evpi <- merge(evpi, w_dt, by = "grp")
totevpi <- evpi[,lapply(.SD, weighted.mean, w = w),
                by = "k", .SDcols = c("evpi")]
ggplot(totevpi, aes(x = k, y = evpi)) +
  geom_line() + xlab("Willingness to pay") +
  ylab("Total EVPI") +
  scale_x_continuous(breaks = seq(0, ktop, 100000), 
                     label = scales::dollar) +
  scale_y_continuous(label = scales::dollar) +
  theme(legend.position = "bottom") 

## ----totenmb------------------------------------------------------------------
# Compute total expected NMB with one-size fits all treatment
ce <- merge(ce, w_dt, by = "grp")
totenmb <- ce[, .(totenmb = weighted.mean(nmb, w = w)), by = c("strategy")]
totenmb_max <- max(totenmb$totenmb)

## ----ptotenmb-----------------------------------------------------------------
# Compute total expected NMB with individualized treatment
itotenmb_grp_max <- apply(as.matrix(enmb[, -1]), 2, max)
itotenmb_max <- sum(itotenmb_grp_max * w_dt$w)

## ----evic2--------------------------------------------------------------------
# Compute EVIC
totnmb_scenarios <- c(itotenmb_max, totenmb_max)
names(totnmb_scenarios) <- c("Individualized total expected NMB",
                              "One-size fits all total expected NMB")
evic <- totnmb_scenarios[1] - totnmb_scenarios[2]
names(evic) <- "EVIC"
print(evic)
print(evic/150000)

