\name{argmax.quantile}
\alias{argmax.quantile}
\title{Quantiles of a simulated Brownian argmax distribution}
\description{
Simulates the Brownian-motion-based argmax distribution used for
change-point inference and returns its empirical 5\% and 95\% quantiles.
}
\usage{
argmax.quantile(sigma1, sigma1_star, sigma2, sigma2_star,
                s_values = seq(-20, 20, length.out = 2 * 1000 + 1),
                nsim = 2000)
}
\arguments{
  \item{sigma1}{A scale parameter for the positive side of the argmax process.}
  \item{sigma1_star}{A drift/scale parameter for the positive side of the argmax process.}
  \item{sigma2}{A scale parameter for the negative side of the argmax process.}
  \item{sigma2_star}{A drift/scale parameter for the negative side of the argmax process.}
  \item{s_values}{The grid of values on which the simulated process is evaluated.}
  \item{nsim}{The number of Monte Carlo simulations.}
}
\value{
A numeric vector containing the empirical 5\% and 95\% quantiles of the
simulated argmax locations.
}
\details{
This function can be supplied to \code{\link{rwbs_fit}} through the
\code{argmax_quantile_fun} argument when confidence intervals are requested.
}
