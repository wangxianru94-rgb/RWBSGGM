\name{hausdorff.distance}
\alias{hausdorff.distance}
\title{Hausdorff distance between two change-point sets}
\description{
Computes the Hausdorff distance between two vectors of change-point
locations.
}
\usage{
hausdorff.distance(v1, v2, n)
}
\arguments{
  \item{v1}{A vector of estimated or candidate change-point locations.}
  \item{v2}{A vector of true or reference change-point locations.}
  \item{n}{The sample size, used as the returned distance when either input set is empty.}
}
\value{
A nonnegative numeric value giving the Hausdorff distance between the two
input sets.
}
