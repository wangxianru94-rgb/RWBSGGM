\name{block.precision.matrix2}
\alias{block.precision.matrix2}
\title{Generate an alternative block-structured precision matrix}
\description{
Generates an alternative block-structured precision matrix used in the
numerical experiments.
}
\usage{
block.precision.matrix2(dd, alpha)
}
\arguments{
  \item{dd}{The dimension of the matrix.}
  \item{alpha}{A parameter controlling the locations that are randomly permuted before the block structure is formed.}
}
\value{
A \code{dd} by \code{dd} block-structured precision matrix.
}
