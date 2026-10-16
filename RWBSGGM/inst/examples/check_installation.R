# Basic installation check for RWBSGGM

library(RWBSGGM)

stopifnot(is.function(rwbs_fit))
stopifnot(is.function(argmax.quantile))
stopifnot(is.function(chain.graph2))
stopifnot(is.function(hausdorff.distance))

cat("RWBSGGM installation check passed.\n")
