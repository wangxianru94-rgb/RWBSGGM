#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript run_time_rwbs.R <n> <p>")
}

n <- as.numeric(args[1])
p <- as.numeric(args[2])
alpha0 <- 0.25 * p

cat("Running RWBS with n =", n, ", p =", p, ", alpha0 =", alpha0, "\n")

time_total_start <- proc.time()[3]

if (!requireNamespace("RWBSGGM", quietly = TRUE)) {
  stop("Package 'RWBSGGM' is not installed. Please install the package first with R CMD INSTALL.")
}

library(RWBSGGM)
library(MASS)

tau0 <- 0.5
delta <- 0.1 * n
delta1 <- 0.1
m1 <- delta
kmax <- 10
c1 <- 1 * sqrt(1 / delta1)
c2 <- 2
M <- 30
c0 <- c1
item <- 1
s_values <- seq(-20, 20, length.out = 2 * 1000 + 1)

graph <- chain.graph2(p, alpha0)
precision1 <- graph[[1]]
precision2 <- graph[[2]]

X1 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision1))
X2 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision2))
X3 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision1))
X4 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision2))
x <- rbind(X1, X2, X3, X4)
Y <- rep(0, p)

# ================== Step 1: initial RWBS + cache construction ==================
time_step1_start <- proc.time()[3]

ran.intervals <- gen.intervals(n, M)
NN <- n
vp.tau <- 0.5 * c2 * p * log(p) / n
record <- cusum.record(ran.intervals, x, delta, c1, c2)
vp.estimate <- vpcusum(NN, ran.intervals, delta, vp.tau, record)

cat("vp =\n")
print(vp.estimate)

alpha <- c(1, vp.estimate, n)
tau.length <- length(alpha)
cp.number <- tau.length - 2
tau.location.initial <- c(rep(0, kmax - length(alpha)), alpha)

# 把 betahat.all 和 cache 都放进 Step 1
if (cp.number == 0) {
  betahat.all <- NULL
  prefix_total <- NULL
} else {
  betahat.all <- estimate_graphs_on_segments(alpha, c1, x)
  prefix_total <- build_refinement_cache(alpha, betahat.all, x)
}

time_step1_end <- proc.time()[3]
step1_time <- time_step1_end - time_step1_start
cat("Step 1 time (initial RWBS + graph caching):", step1_time, "seconds\n")

# ================== Step 2: pure refinement ==================
time_step2_start <- proc.time()[3]

if (cp.number == 0) {
  tau.location.refine <- c(rep(0, kmax - length(alpha)), alpha)
} else {
  tau.location.refine <- refinement_reuse(alpha, prefix_total, kmax = kmax)
}

cpt.detect.res <- numeric()
cpt.detect.res <- c(rep(0, 2 * kmax - length(cpt.detect.res)), cpt.detect.res)

time_step2_end <- proc.time()[3]
step2_time <- time_step2_end - time_step2_start
cat("Step 2 time (pure refinement):", step2_time, "seconds\n")

res <- list()
res[[1]] <- cp.number
res[[2]] <- tau.location.initial
res[[3]] <- tau.location.refine
res[[4]] <- cpt.detect.res

time_total_end <- proc.time()[3]
total_time <- time_total_end - time_total_start
cat("Total time for both steps:", total_time, "seconds\n")

time_summary <- data.frame(
  n = n,
  p = p,
  alpha0 = alpha0,
  step1_time = step1_time,
  step2_time = step2_time,
  total_time = total_time
)

write.csv(
  time_summary,
  sprintf("time_rwbs_n%d_p%d_alpha%.0f_summary.csv", n, p, alpha0),
  row.names = FALSE
)

cat(
  "Finished running, time summary saved to ",
  sprintf("time_rwbs_n%d_p%d_alpha%d_summary.csv", n, p, alpha0),
  "\n",
  sep = ""
)