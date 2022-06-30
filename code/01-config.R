set.seed(101)

# Data are generated to mimic the external control (RCT) and randomized
# randomized controlled trial (RCT) data sets analyzed, without including 
# individual-level data.
external <- c(rep(1, 547),        # EC
              rep(0, 850))        # RCT
treat    <- c(rep(0, 547 + 425),  # control
              rep(1, 425))        # treated
sex2     <- c(rep(0, 243),        # female EC
              rep(1, 304),        # male EC
              rep(0, 330),        # female RCT
              rep(1, 520))        # male RCT 
age <- c(rnorm(243, 66.18, 9.63), # female EC
         rnorm(304, 67.79, 8.50), # male EC
         rnorm(330, 62.64, 9.68), # female RCT
         rnorm(520, 63.52, 9.09)) # male RCT
X.all <- data.frame(external, treat, sex2, age)

# Fitted parameters from proportional hazards model for events 
# (see equations (1) and (2)).
mu2 <- c(-0.316176599, 0.009002796, 0.231969682, -6.832722957, -6.768120447)
names(mu2) <- c(treat, age, sex2, interval1, interval2)

# Fitted parameters from proportional hazards model for censoring 
# (see equation (5)).
theta2 <- c(-5.493429, -6.533711)
names(theta2) <- c("interval1", "interval2")

cut.time <- c(0, 207)
inner.t  <- cut.time[-1]
diff     <- 207 # need to find old code
n.seg    <- 3

interval.names      <- c("interval1", "interval2")
#########################################################################################################

idx_samp   <- 1 # number of simulations (i.e. RCT and external data generated)
idx_n      <- 10  # repetitions per design (i.e. number of imputations per exposure time interval in external controls) # only relevant if n.seg > 2 
idx_a0_n   <- 1   # repetitions per a0 (i.e. number of compatibility scores computed per exposure time interval)
trt_n      <- 200 # number of treated subjects to use
ctrl_n     <- 100 # number of control subjects to use
h_n        <- 100 # number of external (historical) subjects to use
run_com    <- FALSE  # run commensurate prior (takes time)
save_a0    <- FALSE # save weights a0

#########################################################################################################

# scale_list <- rep(rep(rep(c(seq(-log(3), log(3), length = 9), 0), each = 12), times = 10), times = 2)
# T1E_vec    <- rep(rep(rep(c(1, 1, 1, 0, 0), each = 120), times = 2), times = 2)
# X3_vec     <- rep(rep(c(0, 1), each = 600), times = 2)
# cen_vec    <- rep(c(0, 1), each = 1200)

# scale_list <- rep(0, 2400)
# T1E_vec    <- rep(1, 2400)
# X3_vec     <- rep(0, 2400)
# cen_vec    <- rep(0, 2400)

# scale_list <- rep(rep(c(seq(-log(3), log(3), length = 7), 0), each = 50), times = 2)
# T1E_vec    <- rep(c(0, 1), each = 400)
# X3_vec     <- rep(0, 800)
# cen_vec    <- rep(0, 800)

scale_source <- seq(-log(3), log(3), length = 25)
scale_source[5]  <- log(1/2)
scale_source[10] <- log(3/4)
scale_source[16] <- log(4/3)
scale_source[21] <- log(2)
# exp(scale_source)

scale_list <- rep(c(rep(scale_source, each = 8)), times = 12)
T1E_vec    <- rep(rep(c(0, 1), each = 200), times = 6)
X3_vec     <- rep(rep(c(0, 1, 2), each = 400), times = 2)
cen_vec    <- rep(c(0, 1), each = 1200) # 9-7-21

table(scale_list, T1E_vec, X3_vec, cen_vec)
#########################################################################################################

cov.names  <- c("age", "sex2")
a <- c("HR", "MSE", "lower", "upper", "CI_width", "CI", "power", "se",
       cov.names)
b <- c("a0", "0", "0.25", "0.5", "0.75", "1", "f", "iptw", "cox")
c <- c("com_cauchy", "com_a0") #, "com_gamma", "com_no_ext", "com_full_ext")
d <- c("a0_v1")

# eList <- paste0("e", format(seq(0.15, 0.4, by = 0.05), nsmall = 2))
# qList <- paste0("q", format(c(1, 1.25, 1.5, 2, 3), nsmall = 2))

# pqList <- c("p0_q1",
#             "p1_q4",
#             "p1_q2",
#             "p1_q1",
#             "p2_q1",
#             "p4_q1",
#             "p1_q0")
pqList <- c("p1_q0", "p3_q0", "p5_q0", "p7_q0")
tList <- paste("t", format(c(0.950, 0.955, 0.960, 0.965, 0.970, 0.975,
                             0.990, 0.991, 0.993, 0.995, 0.997, 0.999, 0), nsmall = 3), sep = "_")

# eList <- paste(expand.grid(eList, tList)$Var1, expand.grid(eList, tList)$Var2, sep = "_")
# qList <- paste(expand.grid(qList, tList)$Var1, expand.grid(qList, tList)$Var2, sep = "_")
pqList <- paste(expand.grid(pqList, tList)$Var1, expand.grid(pqList, tList)$Var2, sep = "_")

# nList  <- c(b, eList, qList)
# ppList <- c(nList[1:6], eList, qList)
nList  <- c(b, pqList)
ppList <- c(nList[1:6], pqList)

names <- c(paste(expand.grid(a,b)$Var1, expand.grid(a,b)$Var2, sep = "_"),
           paste(expand.grid(a,c)$Var1, expand.grid(a,c)$Var2, sep = "_"),
           paste(expand.grid(a,d)$Var1, expand.grid(a,d)$Var2, sep = "_"),
           "SE_a0", "a0_logHR", "bias_a0", "scale", "mean_a0", 
           "a0_wt_ecog0", "a0_wt_ecog1", "a0_wt_ecogNA", 
           c(paste0("a0_wt_int", 1:(n.seg - 1))),
           c(paste(expand.grid(c("a0_mean", paste0("a0_wt_int", 1:(n.seg - 1))), c("q0.1", "q0.25", "q0.5", "q0.75", "q0.9"))$Var1,
                   expand.grid(c("a0_mean", paste0("a0_wt_int", 1:(n.seg - 1))), c("q0.1", "q0.25", "q0.5", "q0.75", "q0.9"))$Var2, sep = "_")),
           "mean_a0_r", "a0_mean_r_q0.1", "a0_mean_r_q0.25", "a0_mean_r_q0.5", "a0_mean_r_q0.75", "a0_mean_r_q0.9",
           c(paste("a0_wts", expand.grid(c(0,1,2),1:(n.seg - 1))$Var1, 
                   expand.grid(c(0,1,2),1:(n.seg - 1))$Var2, sep = "_")),
           "nu_trt", "nu_ctrl", "nu_h", "T1E",
           "a0_wt_obs", "a0_wt_cen",
           # paste(expand.grid(a,eList)$Var1, expand.grid(a,eList)$Var2, sep = "_"),
           # paste(expand.grid(a,qList)$Var1, expand.grid(a,qList)$Var2, sep = "_")
           paste(expand.grid(a,pqList)$Var1, expand.grid(a,pqList)$Var2, sep = "_"))

########################################################################################################################################################################################################

# rotate <- function(x, q, e){
#   e*x + 0.5 * (1 - e) + q*0
# }
# 
# rotate2 <- function(x, q, e){
#   (((2 * (x - 0.5)) ^ 3 + 1) / 2 + (q - 1 ) * 0.5) / q
# }

rotate3 <- function(x, p, q){
   ((2 * (x - 0.5)) ^ p + 1) / 2
}

g <- function(x, q){
  1 / (1 + exp(-(200*(x - q))))
}
