# Minimal installation and usage example for RWBSGGM

library(RWBSGGM)

set.seed(123)
p <- 10
n <- 80

G <- chain.graph2(p, alpha = 2)

X1 <- MASS::mvrnorm(n / 2, rep(0, p), solve(G$graph1))
X2 <- MASS::mvrnorm(n / 2, rep(0, p), solve(G$graph2))
X <- rbind(X1, X2)

fit <- rwbs_fit(
  x = X,
  M = 5,
  delta = ceiling(0.1 * n),
  delta1 = 0.1,
  kmax = 5,
  do_inference = FALSE
)

print(fit$alpha.initial)
print(fit$alpha.refine)
