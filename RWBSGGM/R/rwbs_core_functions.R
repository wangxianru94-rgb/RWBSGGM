
############################################################
# RWBS Core Functions
#
# This file contains the full implementation of the RWBS procedure,
# including:
#   - Initial estimation via WBS
#   - Refinement via dynamic programming
#   - Confidence interval construction
#
# Main entry:
#   rwbs_fit()
#
############################################################
# ==========================================================
# Step 0: random WBS intervals
# ==========================================================
gen.intervals <- function(N, M) {
  result.list <- matrix(nrow = M, ncol = 2)
  for (i in 1:M) {
    result.list[i, ] <- c(sort(round(runif(2, min = 0, max = N))))
  }
  return(result.list)
}

# ==========================================================
# Step 1: local gain / one-change objective
# ==========================================================
one.change.group.lasso <- function(t, s, e, x, c1, c2) {
  n <- nrow(x)
  p <- ncol(x)
  
  fit11 <- rep(0, p)
  fit21 <- rep(0, p)
  fit31 <- rep(0, p)
  
  for (k in 1:p) {
    Y <- x[, k]
    X <- x[, -k, drop = FALSE]
    
    fit2 <- glmnet(X[s:t, , drop = FALSE], Y[s:t],
                   family = "gaussian", intercept = FALSE, standardize = FALSE)
    fit.value <- predict(fit2, newx = X[s:t, , drop = FALSE],
                         s = c1 * sqrt(log(p) / (t - s + 1)))
    fit21[k] <- sum((Y[s:t] - fit.value)^2) / n
  }
  
  for (k in 1:p) {
    Y <- x[, k]
    X <- x[, -k, drop = FALSE]
    
    fit3 <- glmnet(X[(t + 1):e, , drop = FALSE], Y[(t + 1):e],
                   family = "gaussian", intercept = FALSE, standardize = FALSE)
    fit.value <- predict(fit3, newx = X[(t + 1):e, , drop = FALSE],
                         s = c1 * sqrt(log(p) / (e - t)))
    fit31[k] <- sum((Y[(t + 1):e] - fit.value)^2) / n
  }
  
  for (k in 1:p) {
    Y <- x[, k]
    X <- x[, -k, drop = FALSE]
    
    fit1 <- glmnet(X[s:e, , drop = FALSE], Y[s:e],
                   family = "gaussian", intercept = FALSE, standardize = FALSE)
    fit.value <- predict(fit1, newx = X[s:e, , drop = FALSE],
                         s = c1 * sqrt(log(p) / (e - s + 1)))
    fit11[k] <- sum((Y[s:e] - fit.value)^2) / n
  }
  
  res <- sum(fit11 - fit21 - fit31) +
    c2 * (log(p) / (e - s)) -
    c2 * (log(p) / (e - t)) -
    c2 * (log(p) / (t - s))
  
  return(list(res = res))
}

# ==========================================================
# Step 1: intervalwise search
# ==========================================================
vpcusum.one.interval <- function(s, e, x, delta, c1, c2) {
  projected.data.max <- rep(0, e)
  for (t in (s + delta):(e - delta)) {
    projected.data.max[t] <- one.change.group.lasso(t, s, e, x, c1, c2)$res
  }
  projected.data.max.max <- max(projected.data.max)
  estimate <- which.max(projected.data.max)
  return(list(value = projected.data.max.max, estimate = estimate))
}

cusum.record <- function(ran.intervals, x, delta, c1, c2) {
  M <- nrow(ran.intervals)
  record <- matrix(0, nrow = M, ncol = 2)
  
  for (m in 1:M) {
    s <- ran.intervals[m, 1]
    e <- ran.intervals[m, 2]
    
    if (e - s > 2 * delta) {
      out.one.interval <- vpcusum.one.interval(s, e, x, delta, c1, c2)
      record[m, 1] <- out.one.interval$value
      record[m, 2] <- out.one.interval$estimate
    }
  }
  return(record)
}

vpcusum <- function(NN, ran.intervals, delta, vp.tau, record) {
  start.set <- c(0, NN)
  output.set <- c()
  
  while (length(start.set) > 1) {
    s <- start.set[1]
    e <- start.set[2]
    
    record.temp <- matrix(
      record[which(ran.intervals[, 1] > s & ran.intervals[, 2] < e + 1), ],
      ncol = 2
    )
    
    if (nrow(record.temp) == 0) {
      start.set <- start.set[-1]
    } else {
      out.temp <- record.temp[which.max(record.temp[, 1]), ]
      if (out.temp[1] > vp.tau) {
        start.set <- sort(c(start.set, out.temp[2]))
        output.set <- c(output.set, out.temp[2])
      } else {
        start.set <- start.set[-1]
      }
    }
  }
  return(sort(output.set))
}

# ==========================================================
# Step 2: original refinement
# ==========================================================
refinement <- function(tau.length, alpha, c1, x) {
  n <- nrow(x)
  p <- ncol(x)
  m1 <- ceiling(0.12 * n)
  kmax <- 10
  
  cp.number <- tau.length - 2
  segment.number <- tau.length - 1
  betahat.all <- array(0, dim = c(p, p - 1, segment.number))
  
  for (i in 1:segment.number) {
    for (k in 1:p) {
      Y <- x[, k]
      X <- x[, -k, drop = FALSE]
      fit4 <- glmnet(
        X[alpha[i]:alpha[i + 1], , drop = FALSE],
        Y[alpha[i]:alpha[i + 1]],
        family = "gaussian",
        intercept = TRUE,
        standardize = FALSE
      )
      betahat.all[k, , i] <-
        as.vector(coef(fit4, s = c1 * sqrt(log(p) / (alpha[i + 1] - alpha[i] + 1))))[-1]
    }
  }
  
  f <- matrix(10000, nrow = segment.number, ncol = n)
  positionhat <- matrix(0, nrow = kmax, ncol = n)
  
  fit1 <- rep(0, p)
  for (j in m1:n) {
    for (k in 1:p) {
      Y <- x[, k]
      X <- x[, -k, drop = FALSE]
      fit1[k] <- sum((Y[1:j] - X[1:j, , drop = FALSE] %*% betahat.all[k, , 1])^2) / n
    }
    f[1, j] <- sum(fit1)
  }
  
  fit2 <- rep(0, p)
  for (i in 2:segment.number) {
    for (j in (i * m1):n) {
      fu <- rep(10000, n)
      lf <- rep(10000, n)
      for (z in ((i - 1) * m1):(j - m1)) {
        if (((i - 1) * m1) > (j - m1)) break
        for (k in 1:p) {
          Y <- x[, k]
          X <- x[, -k, drop = FALSE]
          fit2[k] <- sum((Y[(z + 1):j] - X[(z + 1):j, , drop = FALSE] %*% betahat.all[k, , i])^2) / n
        }
        lf[z] <- sum(fit2)
        fu[z] <- f[i - 1, z] + lf[z]
      }
      f[i, j] <- min(fu)
      positionhat[i, j] <- which.min(fu)
    }
  }
  
  lhat <- c(rep(0, cp.number), n)
  for (i in cp.number:1) {
    lhat[i] <- positionhat[i + 1, lhat[i + 1]]
  }
  
  lhat <- c(rep(0, kmax - tau.length), 1, lhat)
  return(lhat)
}

# ==========================================================
# Step 2: faster refinement using cached graph estimates
# ==========================================================
estimate_graphs_on_segments <- function(alpha, c1, x) {
  p <- ncol(x)
  segment.number <- length(alpha) - 1
  betahat.all <- array(0, dim = c(p, p - 1, segment.number))
  
  for (i in seq_len(segment.number)) {
    s <- alpha[i]
    e <- alpha[i + 1]
    lambda_i <- c1 * sqrt(log(p) / (e - s + 1))
    Xi <- x[s:e, , drop = FALSE]
    
    for (k in seq_len(p)) {
      Y <- Xi[, k]
      X <- Xi[, -k, drop = FALSE]
      
      fit <- glmnet(
        X, Y,
        family = "gaussian",
        intercept = FALSE,
        standardize = FALSE
      )
      
      betahat.all[k, , i] <- as.vector(coef(fit, s = lambda_i))[-1]
    }
  }
  
  betahat.all
}

build_refinement_cache <- function(alpha, betahat.all, x) {
  n <- nrow(x)
  p <- ncol(x)
  segment.number <- length(alpha) - 1
  prefix_total <- matrix(0, nrow = segment.number, ncol = n + 1)
  
  for (i in seq_len(segment.number)) {
    resid_sum_t <- rep(0, n)
    
    for (k in seq_len(p)) {
      Y_all <- x[, k]
      X_all <- x[, -k, drop = FALSE]
      beta_ik <- betahat.all[k, , i]
      r <- as.vector(Y_all - X_all %*% beta_ik)
      resid_sum_t <- resid_sum_t + (r^2) / n
    }
    
    prefix_total[i, 2:(n + 1)] <- cumsum(resid_sum_t)
  }
  
  prefix_total
}

refinement_reuse <- function(alpha, prefix_total, kmax = 10) {
  n <- ncol(prefix_total) - 1
  tau.length <- length(alpha)
  cp.number <- tau.length - 2
  segment.number <- tau.length - 1
  m1 <- ceiling(0.12 * n)
  
  interval_cost <- function(seg_id, a, b) {
    prefix_total[seg_id, b + 1] - prefix_total[seg_id, a]
  }
  
  f <- matrix(Inf, nrow = segment.number, ncol = n)
  positionhat <- matrix(0, nrow = kmax, ncol = n)
  
  for (j in m1:n) {
    f[1, j] <- interval_cost(1, 1, j)
  }
  
  if (segment.number >= 2) {
    for (i in 2:segment.number) {
      for (j in (i * m1):n) {
        z_candidates <- ((i - 1) * m1):(j - m1)
        cost_vec <- f[i - 1, z_candidates] +
          sapply(z_candidates, function(z) interval_cost(i, z + 1, j))
        best_idx <- which.min(cost_vec)
        f[i, j] <- cost_vec[best_idx]
        positionhat[i, j] <- z_candidates[best_idx]
      }
    }
  }
  
  lhat <- c(rep(0, cp.number), n)
  if (cp.number >= 1) {
    for (i in cp.number:1) {
      lhat[i] <- positionhat[i + 1, lhat[i + 1]]
    }
  }
  
  c(rep(0, kmax - tau.length), 1, lhat)
}

# ==========================================================
# utilities for post-refinement graph estimation
# ==========================================================
build_segment_windows <- function(alpha) {
  segment.number <- length(alpha) - 1
  Deltahat <- min(diff(alpha))
  s.j.all <- rep(0, 2 * segment.number)
  
  for (i in 1:(2 * segment.number)) {
    if (i %% 2 != 0) s.j.all[i] <- alpha[ceiling(i / 2)] + ceiling(Deltahat / 10)
    if (i %% 2 == 0) s.j.all[i] <- alpha[ceiling(i / 2) + 1] - ceiling(Deltahat / 10)
  }
  
  for (i in 1:segment.number) {
    if (s.j.all[2 * i] - s.j.all[2 * i - 1] < 2) {
      s.j.all[2 * i - 1] <- s.j.all[2 * i - 1] - 1
      s.j.all[2 * i] <- s.j.all[2 * i] + 1
    }
  }
  
  s.j.all
}

estimate_segment_graphs_cv <- function(x, alpha) {
  p <- ncol(x)
  tau.length <- length(alpha)
  segment.number <- tau.length - 1
  cp.number <- tau.length - 2
  
  betahat.all <- array(0, dim = c(p, p - 1, segment.number))
  deltahat.all <- array(0, dim = c(p, p - 1, cp.number))
  s.j.all <- build_segment_windows(alpha)
  
  for (k in 1:p) {
    Y <- x[, k]
    X <- x[, -k, drop = FALSE]
    
    for (i in 1:segment.number) {
      X.true <- X[s.j.all[2 * i - 1]:s.j.all[2 * i], , drop = FALSE]
      Y.true <- Y[s.j.all[2 * i - 1]:s.j.all[2 * i]]
      lambda_candidate <- seq(2^0, 2^2, 0.05) *
        sqrt(log(p) / (s.j.all[2 * i] - s.j.all[2 * i - 1] + 1))
      
      cv_fit <- cv.glmnet(
        X.true, Y.true,
        family = "gaussian",
        intercept = TRUE,
        standardize = FALSE,
        lambda = lambda_candidate,
        nfolds = 3
      )
      
      best_lambda <- cv_fit$lambda.min
      fit_best <- glmnet(X.true, Y.true, alpha = 1, lambda = best_lambda)
      betahat.all[k, , i] <- coef(fit_best)[-1]
    }
    
    for (i in 1:cp.number) {
      deltahat.all[k, , i] <- betahat.all[k, , i + 1] - betahat.all[k, , i]
    }
  }
  
  list(
    betahat.all = betahat.all,
    deltahat.all = deltahat.all,
    s.j.all = s.j.all
  )
}

compute_sigma_components <- function(x, alpha, betahat.all, deltahat.all, s.j.all) {
  n <- nrow(x)
  p <- ncol(x)
  segment.number <- length(alpha) - 1
  cp.number <- length(alpha) - 2
  
  Deltahat.all <- rep(0, cp.number)
  for (i in 1:cp.number) {
    Deltahat.all[i] <- norm(deltahat.all[, , i], "F")
  }
  pshihat.all <- Deltahat.all / sqrt(p)
  
  epsilonhat.all <- matrix(0, nrow = p, ncol = n)
  for (k in 1:p) {
    Y <- x[, k]
    X <- x[, -k, drop = FALSE]
    for (i in 1:segment.number) {
      epsilonhat.all[k, alpha[i]:alpha[i + 1]] <-
        Y[alpha[i]:alpha[i + 1]] -
        X[alpha[i]:alpha[i + 1], , drop = FALSE] %*% betahat.all[k, , i]
    }
  }
  
  Sigmahat.all <- array(0, dim = c(p, p, segment.number))
  for (i in 1:segment.number) {
    Sigmahat.all[, , i] <- cov(x[s.j.all[2 * i - 1]:s.j.all[2 * i], , drop = FALSE])
  }
  
  sigmastarhat.all <- matrix(0, nrow = cp.number, ncol = 2)
  sigma.all1.all <- matrix(0, nrow = cp.number, ncol = p)
  sigma.all2.all <- matrix(0, nrow = cp.number, ncol = p)
  
  for (i in 1:cp.number) {
    for (k in 1:p) {
      sigma.all1.all[i, k] <-
        deltahat.all[k, , i] %*% Sigmahat.all[-k, -k, i + 1] %*% deltahat.all[k, , i] /
        Deltahat.all[i]^2
      sigma.all2.all[i, k] <-
        deltahat.all[k, , i] %*% Sigmahat.all[-k, -k, i] %*% deltahat.all[k, , i] /
        Deltahat.all[i]^2
    }
    sigmastarhat.all[i, ] <- c(
      sqrt(sum(sigma.all1.all[i, ])),
      sqrt(sum(sigma.all2.all[i, ]))
    )
  }
  
  R.all <- matrix(0, nrow = n, ncol = cp.number)
  for (j in 1:cp.number) {
    for (i in 1:n) {
      tmp <- 0
      for (k in 1:p) {
        tmp <- tmp + epsilonhat.all[k, i] * x[i, -k] %*% deltahat.all[k, , j]
      }
      R.all[i, j] <- tmp / sqrt(p)
    }
  }
  
  Rbar.all <- matrix(0, nrow = cp.number, ncol = 2)
  for (i in 1:cp.number) {
    Rbar.all[i, 1] <- mean(R.all[s.j.all[2 * i + 1]:s.j.all[2 * i + 2], i])
    Rbar.all[i, 2] <- mean(R.all[s.j.all[2 * i - 1]:s.j.all[2 * i], i])
  }
  
  sigmahat.all <- matrix(0, nrow = cp.number, ncol = 2)
  for (i in 1:cp.number) {
    sigmahat.all[i, 1] <-
      sqrt(mean((R.all[s.j.all[2 * i + 1]:s.j.all[2 * i + 2], i] - Rbar.all[i, 1])^2) /
             Deltahat.all[i]^2)
    sigmahat.all[i, 2] <-
      sqrt(mean((R.all[s.j.all[2 * i - 1]:s.j.all[2 * i], i] - Rbar.all[i, 2])^2) /
             Deltahat.all[i]^2)
  }
  
  list(
    sigmahat.all = sigmahat.all,
    sigmastarhat.all = sigmastarhat.all,
    pshihat.all = pshihat.all,
    epsilonhat.all = epsilonhat.all,
    Sigmahat.all = Sigmahat.all,
    R.all = R.all
  )
}

construct_confidence_intervals <- function(alpha, sigmahat.all, sigmastarhat.all,
                                           pshihat.all, argmax_quantile_fun) {
  cp.number <- length(alpha) - 2
  
  confidence <- matrix(nrow = cp.number, ncol = 2)
  quantile <- matrix(nrow = cp.number, ncol = 2)
  
  for (i in 1:cp.number) {
    quantile[i, ] <- argmax_quantile_fun(
      sigma1 = sigmahat.all[i, 1],
      sigma1_star = sigmastarhat.all[i, 1],
      sigma2 = sigmahat.all[i, 2],
      sigma2_star = sigmastarhat.all[i, 2]
    )
    
    confidence[i, ] <- c(
      alpha[i + 1] - as.numeric(quantile[i, 2]) / pshihat.all[i]^2,
      alpha[i + 1] - as.numeric(quantile[i, 1]) / pshihat.all[i]^2
    )
  }
  
  list(confidence = confidence, quantile = quantile)
}

evaluate_detection_and_ci <- function(alpha, confidence, n, kmax = 10) {
  cpt.true <- c(1, 0.25 * n, 0.5 * n, 0.75 * n, n)
  cpt.true.actual <- cpt.true[-c(1, length(cpt.true))]
  
  cpt.est.refine <- alpha
  cpt.est.refine.actual <- alpha[-c(1, length(alpha))]
  
  detect.distance <- floor(0.25 * min(diff(cpt.true)))
  cpt.detect.res <- numeric()
  
  if (length(cpt.est.refine) > 2) {
    for (i in 1:length(cpt.true.actual)) {
      abs_diff <- abs(cpt.est.refine.actual - cpt.true.actual[i])
      closest.index <- which.min(abs_diff)
      if (abs_diff[closest.index] <= detect.distance) {
        cpt.detect.res[(i - 1) * 2 + 1] <- 1
        if ((cpt.true.actual[i] >= confidence[closest.index, 1]) &&
            (cpt.true.actual[i] <= confidence[closest.index, 2])) {
          cpt.detect.res[(i - 1) * 2 + 2] <- 1
        } else {
          cpt.detect.res[(i - 1) * 2 + 2] <- 0
        }
      } else {
        cpt.detect.res[(i - 1) * 2 + 1] <- 0
        cpt.detect.res[(i - 1) * 2 + 2] <- 0
      }
    }
  } else {
    cpt.detect.res <- rep(0, 2 * length(cpt.true.actual))
  }
  
  cpt.detect.res <- c(rep(0, 2 * kmax - length(cpt.detect.res)), cpt.detect.res)
  cpt.detect.res
}

# ==========================================================
# unified main function
# ==========================================================
rwbs_fit <- function(x, M = 30, delta = ceiling(0.1 * nrow(x)),
                     delta1 = 0.1, kmax = 10,
                     c1 = 0.75 * sqrt(1 / delta1),
                     c2 = 0.75,
                     vp.tau = NULL,
                     use_fast_refinement = FALSE,
                     do_inference = TRUE,
                     argmax_quantile_fun = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  
  if (is.null(vp.tau)) {
    vp.tau <- 0.2 * c2 * p * log(p) / n
  }
  
  ran.intervals <- gen.intervals(n, M)
  record <- cusum.record(ran.intervals, x, delta, c1, c2)
  vp.estimate <- vpcusum(n, ran.intervals, delta, vp.tau, record)
  alpha.initial <- c(1, vp.estimate, n)
  
  tau.length <- length(alpha.initial)
  cp.number <- tau.length - 2
  tau.location.initial <- c(rep(0, kmax - length(alpha.initial)), alpha.initial)
  
  if (cp.number == 0) {
    alpha.refine <- alpha.initial
    tau.location.refine <- c(rep(0, kmax - length(alpha.refine)), alpha.refine)
  } else {
    if (use_fast_refinement) {
      betahat.cache <- estimate_graphs_on_segments(alpha.initial, c1, x)
      prefix_total <- build_refinement_cache(alpha.initial, betahat.cache, x)
      tau.location.refine <- refinement_reuse(alpha.initial, prefix_total, kmax = kmax)
    } else {
      tau.location.refine <- refinement(tau.length, alpha.initial, c1, x)
    }
    alpha.refine <- tau.location.refine[tau.location.refine != 0]
  }
  
  out <- list(
    alpha.initial = alpha.initial,
    tau.location.initial = tau.location.initial,
    alpha.refine = alpha.refine,
    tau.location.refine = tau.location.refine
  )
  
  if (do_inference) {
    if (is.null(argmax_quantile_fun)) {
      stop("Please provide argmax_quantile_fun when do_inference = TRUE.")
    }
    
    graph.est <- estimate_segment_graphs_cv(x, alpha.refine)
    sigma.comp <- compute_sigma_components(
      x = x,
      alpha = alpha.refine,
      betahat.all = graph.est$betahat.all,
      deltahat.all = graph.est$deltahat.all,
      s.j.all = graph.est$s.j.all
    )
    
    ci.out <- construct_confidence_intervals(
      alpha = alpha.refine,
      sigmahat.all = sigma.comp$sigmahat.all,
      sigmastarhat.all = sigma.comp$sigmastarhat.all,
      pshihat.all = sigma.comp$pshihat.all,
      argmax_quantile_fun = argmax_quantile_fun
    )
    
    detect.res <- evaluate_detection_and_ci(
      alpha = alpha.refine,
      confidence = ci.out$confidence,
      n = n,
      kmax = kmax
    )
    
    out$betahat.all <- graph.est$betahat.all
    out$deltahat.all <- graph.est$deltahat.all
    out$s.j.all <- graph.est$s.j.all
    out$sigmahat.all <- sigma.comp$sigmahat.all
    out$sigmastarhat.all <- sigma.comp$sigmastarhat.all
    out$pshihat.all <- sigma.comp$pshihat.all
    out$quantile <- ci.out$quantile
    out$confidence <- ci.out$confidence
    out$cpt.detect.res <- detect.res
  }
  
  return(out)
}