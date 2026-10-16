\name{chain.graph}
\alias{chain.graph}
\title{Generate two chain-type precision matrices}
\description{
Generates two standardized chain-type precision matrices for simulation
studies, with a controlled difference between the two graph structures.
}
\usage{
chain.graph(p, alpha)
}
\arguments{
  \item{p}{The number of nodes.}
  \item{alpha}{A parameter controlling the portion of locations permuted between the two graph structures.}
}
\value{
A list with components \code{graph1} and \code{graph2}, containing the two
generated standardized precision matrices.
}
