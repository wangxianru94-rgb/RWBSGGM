\name{band.matrix}
\alias{band.matrix}
\title{Generate a band-type covariance matrix}
\description{
Generates a square band-type covariance matrix whose entries decay
geometrically with the distance from the diagonal.
}
\usage{
band.matrix(dd)
}
\arguments{
  \item{dd}{The dimension of the matrix.}
}
\value{
A \code{dd} by \code{dd} symmetric matrix with diagonal entries equal to 1
and off-diagonal entries determined by the decay factor 0.8.
}
