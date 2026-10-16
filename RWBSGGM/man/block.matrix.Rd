\name{block.matrix}
\alias{block.matrix}
\title{Generate a block-structured matrix}
\description{
Generates a block-structured matrix used in the simulation settings.
Variables are grouped in blocks of size five, with a random diagonal.
}
\usage{
block.matrix(dd)
}
\arguments{
  \item{dd}{The dimension of the matrix.}
}
\value{
A \code{dd} by \code{dd} block-structured matrix.
}
