\name{rwbs_fit}
\alias{rwbs_fit}
\title{RWBS estimation, refinement, and confidence intervals}
\description{
Implements the RWBS procedure for multiple change-point estimation in
high-dimensional Gaussian graphical models. The function first obtains
initial change-point estimates using randomly sampled intervals and then
refines the estimated locations. When requested, it additionally estimates
segment-specific graphical models and constructs confidence intervals for
the estimated change-point locations.
}
\usage{
rwbs_fit(x, M = 30, delta = ceiling(0.1 * nrow(x)), delta1 = 0.1,
         kmax = 10, c1 = 0.75 * sqrt(1/delta1), c2 = 0.75,
         vp.tau = NULL, use_fast_refinement = FALSE,
         do_inference = TRUE, argmax_quantile_fun = NULL)
}
\arguments{
  \item{x}{A numeric matrix with observations in rows and variables in columns.}
  \item{M}{The number of randomly generated intervals used in the initial RWBS search.}
  \item{delta}{The minimum interval length used in the initial change-point search.}
  \item{delta1}{A tuning parameter used in the default value of \code{c1}.}
  \item{kmax}{The maximum number of change points represented in the padded location output.}
  \item{c1}{The tuning constant used in the nodewise graphical-model estimation step.}
  \item{c2}{The tuning constant used in the segmentation criterion.}
  \item{vp.tau}{The threshold for the initial change-point selection. If \code{NULL}, a data-dependent default is used.}
  \item{use_fast_refinement}{Logical; if \code{TRUE}, use the cached implementation of the refinement step.}
  \item{do_inference}{Logical; if \code{TRUE}, construct confidence intervals in addition to estimating change points.}
  \item{argmax_quantile_fun}{A function used to obtain the quantiles of the limiting argmax distribution. It is required when \code{do_inference = TRUE}; \code{\link{argmax.quantile}} can be supplied.}
}
\value{
A list containing the initial and refined change-point estimates. The main
components are \code{alpha.initial}, \code{tau.location.initial},
\code{alpha.refine}, and \code{tau.location.refine}. When
\code{do_inference = TRUE}, the returned list additionally contains
segment-wise graphical-model estimates and the quantities used to construct
confidence intervals, including \code{confidence}, \code{quantile}, and
\code{cpt.detect.res}.
}
\details{
The returned \code{alpha.initial} and \code{alpha.refine} vectors include the
sample boundaries 1 and \code{n}. The \code{tau.location.initial} and
\code{tau.location.refine} objects are zero-padded representations with
length controlled by \code{kmax}. For inference, the function requires an
argmax-quantile function; the package function \code{argmax.quantile} provides
the corresponding simulation-based approximation used in the numerical
implementation.
}
\seealso{
\code{\link{argmax.quantile}}
}
