############################################################
# Main script for reproducing simulation results in the paper
#
# Usage:
#   Rscript run_RWBS.R <n> <p> <alpha0>
#
# Example:
#   Rscript run_RWBS.R 600 20 4
#
# This script:
#   1. generates synthetic data under the model setting;
#   2. runs RWBS estimation, refinement, and confidence interval construction;
#   3. repeats the experiment B times in parallel;
#   4. writes all outputs to CSV files.
############################################################

##############################
# 1. Read command-line arguments
##############################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_RWBS.R <n> <p> <alpha0>")
}

n <- as.numeric(args[1])
p <- as.numeric(args[2])
alpha0 <- as.numeric(args[3])

##############################
# 2. Load required packages
##############################
if (!requireNamespace("RWBSGGM", quietly = TRUE)) {
  stop("Package 'RWBSGGM' is not installed. Please install the package first with R CMD INSTALL.")
}

library(RWBSGGM)
library(foreach)
library(doParallel)
library(MASS)

##############################
# 5. Simulation settings
##############################
B <- 100                  # number of replications
n_cores <- 20            # number of parallel workers

delta <- 0.1 * n         # minimal spacing parameter
delta1 <- 0.1
kmax <- 10               # upper bound for number of change points
c1 <- 0.75 * sqrt(1 / delta1)
c2 <- 0.75
M <- 30                  # number of random WBS intervals

# True change-point locations used for Hausdorff distance evaluation
true_cp <- c(1, floor(0.25 * n), floor(0.5 * n), floor(0.75 * n), n)

##############################
# 6. Helper functions
##############################

# Convert a vector of estimated change points to a standardized format:
# - remove missing values
# - keep positive locations only
# - sort and enforce boundary points 1 and n
standardize_rwbs_cp <- function(cp_vec, n) {
  cp_vec <- as.numeric(cp_vec)
  cp_vec <- cp_vec[!is.na(cp_vec)]
  cp_vec <- cp_vec[cp_vec > 0]
  cp_vec <- unique(cp_vec)
  cp_vec <- sort(cp_vec)
  
  if (length(cp_vec) == 0 || cp_vec[1] != 1) {
    cp_vec <- c(1, cp_vec)
  }
  if (tail(cp_vec, 1) != n) {
    cp_vec <- c(cp_vec, n)
  }
  
  cp_vec
}

# Flatten a confidence interval matrix into a fixed-length numeric vector.
# Output length is 2*kmax so that all replications can be saved into a table.
flatten_confidence <- function(confidence, kmax) {
  out <- rep(NA_real_, 2 * kmax)
  
  if (is.null(confidence) || all(is.na(confidence))) {
    return(out)
  }
  
  confidence <- as.matrix(confidence)
  if (nrow(confidence) == 0) {
    return(out)
  }
  
  vals <- as.vector(t(confidence))  # row-wise: l1,u1,l2,u2,...
  keep_len <- min(length(vals), 2 * kmax)
  out[seq_len(keep_len)] <- vals[seq_len(keep_len)]
  out
}

# Convert the 3D array betahat.all into a long data frame for CSV export.
# Each row corresponds to one node within one segment.
flatten_betahat <- function(betahat.all) {
  if (is.null(betahat.all) || all(is.na(betahat.all))) {
    return(data.frame())
  }
  
  dims <- dim(betahat.all)
  if (length(dims) != 3) {
    return(data.frame())
  }
  
  p_dim <- dims[1]
  pm1_dim <- dims[2]
  seg_dim <- dims[3]
  
  rows <- vector("list", seg_dim)
  
  for (s in seq_len(seg_dim)) {
    mat_s <- betahat.all[, , s, drop = FALSE]
    mat_s <- matrix(mat_s, nrow = p_dim, ncol = pm1_dim)
    
    df_s <- as.data.frame(mat_s)
    colnames(df_s) <- paste0("beta", seq_len(pm1_dim))
    df_s$node <- seq_len(p_dim)
    df_s$segment <- s
    rows[[s]] <- df_s
  }
  
  do.call(rbind, rows)
}

# Convert a list of change-point vectors into a rectangular table.
make_cp_table <- function(cp_list, B) {
  cp_lengths <- sapply(cp_list, length)
  max_cp <- max(cp_lengths, na.rm = TRUE)
  cp_matrix <- matrix(NA_real_, nrow = B, ncol = max_cp)
  
  for (b in seq_len(B)) {
    cp_b <- cp_list[[b]]
    cp_matrix[b, seq_len(length(cp_b))] <- cp_b
  }
  
  colnames(cp_matrix) <- paste0("cp", seq_len(max_cp))
  data.frame(replication = seq_len(B), cp_matrix)
}

# Safe extraction helper for change-point vectors
safe_cp <- function(x) {
  if (is.null(x) || all(is.na(x))) {
    return(c(NA_real_))
  }
  as.numeric(x)
}

##############################
# 7. Set up parallel workers
##############################
cl <- makeCluster(n_cores)

invisible(clusterEvalQ(cl, {
  library(MASS)
  library(RWBSGGM)
}))

registerDoParallel(cl)

cat("Running RWBS with n =", n, ", p =", p, ", alpha0 =", alpha0, "\n")
cat("Parallel cores used:", n_cores, "\n")

##############################
# 8. Main simulation loop
##############################
res_list <- foreach(
  b = 1:B,
  .packages = c("MASS", "RWBSGGM"),
  .export = c(
    "n", "p", "alpha0", "delta", "delta1", "kmax",
    "c1", "c2", "M", "true_cp",
    "standardize_rwbs_cp",
    "flatten_confidence",
    "flatten_betahat"
  ),
  .errorhandling = "pass"
) %dopar% {
  
  # Use a different seed for each replication
  set.seed(1000 + b)
  
  ###################################
  # Generate data under Model Setting
  ###################################
  graph <- chain.graph2(p, alpha0)
  precision1 <- graph[[1]]
  precision2 <- graph[[2]]
  
  X1 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision1))
  X2 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision2))
  X3 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision1))
  X4 <- mvrnorm(floor(0.25 * n), rep(0, p), solve(precision2))
  x <- rbind(X1, X2, X3, X4)
  
  ###################################
  # Run unified RWBS procedure
  ###################################
  fit_res <- tryCatch(
    {
      rwbs_fit(
        x = x,
        M = M,
        delta = delta,
        delta1 = delta1,
        kmax = kmax,
        c1 = c1,
        c2 = c2,
        vp.tau = 0.2 * c2 * p * log(p) / n,
        use_fast_refinement = FALSE,
        do_inference = TRUE,
        argmax_quantile_fun = argmax.quantile
      )
    },
    error = function(e) e
  )
  
  ###################################
  # If one replication fails, keep NA outputs
  ###################################
  if (inherits(fit_res, "error")) {
    return(list(
      init_cp = NA_real_,
      refine_cp = NA_real_,
      init_khat = NA_integer_,
      refine_khat = NA_integer_,
      init_hausdorff = NA_real_,
      refine_hausdorff = NA_real_,
      cpt.detect.res = rep(NA_real_, 2 * kmax),
      confidence = rep(NA_real_, 2 * kmax),
      betahat.all = NULL,
      error_msg = conditionMessage(fit_res)
    ))
  }
  
  ###################################
  # Standardize estimated change points
  ###################################
  init_cp_std <- standardize_rwbs_cp(fit_res$alpha.initial, n)
  refine_cp_std <- standardize_rwbs_cp(fit_res$alpha.refine, n)
  
  ###################################
  # Compute Hausdorff distances
  ###################################
  init_hd <- hausdorff.distance(
    v1 = init_cp_std,
    v2 = true_cp,
    n = n
  ) / n
  
  refine_hd <- hausdorff.distance(
    v1 = refine_cp_std,
    v2 = true_cp,
    n = n
  ) / n
  
  ###################################
  # Flatten confidence intervals for export
  ###################################
  conf_vec <- flatten_confidence(fit_res$confidence, kmax = kmax)
  
  ###################################
  # Collect outputs for this replication
  ###################################
  list(
    init_cp = init_cp_std,
    refine_cp = refine_cp_std,
    init_khat = length(init_cp_std) - 2,
    refine_khat = length(refine_cp_std) - 2,
    init_hausdorff = init_hd,
    refine_hausdorff = refine_hd,
    cpt.detect.res = fit_res$cpt.detect.res,
    confidence = conf_vec,
    betahat.all = fit_res$betahat.all,
    error_msg = NA_character_
  )
}

stopCluster(cl)

##############################
# 9. Organize outputs
##############################
init_cp_results <- lapply(res_list, function(z) {
  if (inherits(z, "error")) return(c(NA_real_))
  safe_cp(z$init_cp)
})

refine_cp_results <- lapply(res_list, function(z) {
  if (inherits(z, "error")) return(c(NA_real_))
  safe_cp(z$refine_cp)
})

init_khat <- sapply(res_list, function(z) {
  if (inherits(z, "error")) return(NA_integer_)
  z$init_khat
})

refine_khat <- sapply(res_list, function(z) {
  if (inherits(z, "error")) return(NA_integer_)
  z$refine_khat
})

init_hausdorff <- sapply(res_list, function(z) {
  if (inherits(z, "error")) return(NA_real_)
  z$init_hausdorff
})

refine_hausdorff <- sapply(res_list, function(z) {
  if (inherits(z, "error")) return(NA_real_)
  z$refine_hausdorff
})

cpt_detect_res_mat <- do.call(rbind, lapply(res_list, function(z) {
  if (inherits(z, "error")) return(rep(NA_real_, 2 * kmax))
  z$cpt.detect.res
}))

confidence_mat <- do.call(rbind, lapply(res_list, function(z) {
  if (inherits(z, "error")) return(rep(NA_real_, 2 * kmax))
  z$confidence
}))

error_msg <- lapply(res_list, function(z) {
  if (inherits(z, "error")) return(as.character(z))
  z$error_msg
})

##############################
# 10. Convert outputs to tables
##############################
init_cp_table <- make_cp_table(init_cp_results, B)
refine_cp_table <- make_cp_table(refine_cp_results, B)

init_khat_table <- data.frame(
  replication = seq_len(B),
  k_hat = init_khat
)

refine_khat_table <- data.frame(
  replication = seq_len(B),
  k_hat = refine_khat
)

init_hausdorff_table <- data.frame(
  replication = seq_len(B),
  hausdorff = init_hausdorff
)

refine_hausdorff_table <- data.frame(
  replication = seq_len(B),
  hausdorff = refine_hausdorff
)

cpt_detect_res_table <- data.frame(
  replication = seq_len(B),
  cpt_detect_res_mat
)
colnames(cpt_detect_res_table) <- c(
  "replication",
  paste0("detect_ci_", seq_len(2 * kmax))
)

confidence_table <- data.frame(
  replication = seq_len(B),
  confidence_mat
)
colnames(confidence_table) <- c(
  "replication",
  paste0("conf_", seq_len(2 * kmax))
)

error_table <- data.frame(
  replication = seq_len(B),
  error_msg = unlist(error_msg)
)

##############################
# 11. Organize betahat.all
##############################
betahat_tables <- lapply(seq_len(B), function(b) {
  z <- res_list[[b]]
  if (inherits(z, "error") || is.null(z$betahat.all)) {
    return(NULL)
  }
  
  df_b <- flatten_betahat(z$betahat.all)
  if (nrow(df_b) == 0) {
    return(NULL)
  }
  
  df_b$replication <- b
  df_b
})

betahat_table <- do.call(rbind, betahat_tables)

##############################
# 12. Write outputs to CSV
##############################
write.csv(
  init_cp_table,
  sprintf("RWBS_n%d_p%d_alpha%d_initial_cp_estimates.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  refine_cp_table,
  sprintf("RWBS_n%d_p%d_alpha%d_refine_cp_estimates.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  init_khat_table,
  sprintf("RWBS_n%d_p%d_alpha%d_initial_khat.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  refine_khat_table,
  sprintf("RWBS_n%d_p%d_alpha%d_refine_khat.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  init_hausdorff_table,
  sprintf("RWBS_n%d_p%d_alpha%d_initial_hausdorff.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  refine_hausdorff_table,
  sprintf("RWBS_n%d_p%d_alpha%d_refine_hausdorff.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  cpt_detect_res_table,
  sprintf("RWBS_n%d_p%d_alpha%d_cpt_detect_res.csv", n, p, alpha0),
  row.names = FALSE
)

write.csv(
  confidence_table,
  sprintf("RWBS_n%d_p%d_alpha%d_confidence.csv", n, p, alpha0),
  row.names = FALSE
)

if (!is.null(betahat_table) && nrow(betahat_table) > 0) {
  write.csv(
    betahat_table,
    sprintf("RWBS_n%d_p%d_alpha%d_betahat_all.csv", n, p, alpha0),
    row.names = FALSE
  )
}

write.csv(
  error_table,
  sprintf("RWBS_n%d_p%d_alpha%d_error.csv", n, p, alpha0),
  row.names = FALSE
)

cat("Finished n =", n, ", p =", p, ", alpha0 =", alpha0, "\n")