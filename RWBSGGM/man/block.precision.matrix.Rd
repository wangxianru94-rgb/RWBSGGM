\name{block.precision.matrix}
\alias{block.precision.matrix}
\title{Generate a block-structured precision matrix}
\description{
Generates a block-structured precision matrix used in the numerical
experiments.
}
\usage{
block.precision.matrix(dd, alpha)
}
\arguments{
  \item{dd}{The dimension of the matrix.}
  \item{alpha}{A parameter controlling the subset of locations that is randomly permuted before the block structure is formed.}
}
\value{
A \code{dd} by \code{dd} block-structured precision matrix.
}
