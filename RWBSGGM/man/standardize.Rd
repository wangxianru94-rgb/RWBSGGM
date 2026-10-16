\name{standardize}
\alias{standardize}
\title{Column-standardize a data matrix}
\description{
Centers and scales each column of a data matrix using its sample mean and
sample standard deviation.
}
\usage{
standardize(X)
}
\arguments{
  \item{X}{A numeric data matrix with observations in rows and variables in columns.}
}
\value{
A matrix with each column centered by its sample mean and scaled by its
sample standard deviation.
}
