% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindBestTransform.R
\name{FindBestTransform}
\alias{FindBestTransform}
\title{Find the Power Transformation Parameter.}
\usage{
FindBestTransform(x, grid.length = 50)
}
\arguments{
\item{x}{an n-by-p data frame or matrix of count data. Samples should be in the rows.}

\item{grid.length}{how many distinct points of alpha should be searched within the interval [0, 1]? Default is 50.}
}
\value{
the value of alpha to be used within the power transformation.
}
\description{
Use this function to find a constant value of alpha to be used for transforming count data. The power transformation parameter \code{alpha}, which approximately fits transformed data to the Poisson log-linear model, is selected using a grid search within the interval [0, 1].
}
\note{
This function is copied from \code{PoiClaClu} package and modified to control the total number of grid search.
}
\examples{
set.seed(2128)
counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
                            allZero.rm = FALSE, tag.samples = TRUE)

x <- counts$x
FindBestTransform(x)

}
\author{
Dincer Goksuluk
}
