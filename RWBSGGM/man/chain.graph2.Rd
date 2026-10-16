\name{chain.graph2}
\alias{chain.graph2}
\title{Generate an alternative pair of chain-type precision matrices}
\description{
Generates an alternative pair of standardized chain-type precision matrices
for simulation studies.
}
\usage{
chain.graph2(p, alpha)
}
\arguments{
  \item{p}{The number of nodes.}
  \item{alpha}{A parameter controlling the locations permuted between the two graph structures.}
}
\value{
A list with components \code{graph1} and \code{graph2}, containing the two
generated standardized precision matrices.
}
