\name{standardize_cp_format}
\alias{standardize_cp_format}
\title{Standardize change-point location format}
\description{
Converts a vector of change-point locations to the boundary-inclusive
format used by the package.
}
\usage{
standardize_cp_format(cp_vec, n, has_boundary = FALSE)
}
\arguments{
  \item{cp_vec}{A numeric or integer vector of change-point locations.}
  \item{n}{The sample size.}
  \item{has_boundary}{Logical; if \code{TRUE}, the input is treated as already boundary-inclusive and is sorted, deduplicated, and truncated to the sample range.}
}
\value{
If \code{has_boundary = FALSE}, the returned vector contains the internal
change points together with the two sample boundaries, 1 and \code{n}. If
\code{has_boundary = TRUE}, the supplied locations are returned after
basic boundary and uniqueness processing.
}
