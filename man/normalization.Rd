% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz_methods.R
\docType{methods}
\name{normalization}
\alias{normalization}
\alias{normalization,nblda-method}
\alias{normalization,nblda_trained-method}
\title{Accessors for the 'type' slot.}
\usage{
\S4method{normalization}{nblda}(object)

\S4method{normalization}{nblda_trained}(object)
}
\arguments{
\item{object}{an \code{nblda} or \code{nblda_trained} object.}
}
\description{
This slot stores the name of normalization method. Normalization is defined using \code{type} argument in \code{\link{trainNBLDA}} function.
}
\examples{
set.seed(2128)
counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
                            allZero.rm = FALSE, tag.samples = TRUE)
x <- t(counts$x + 1)
y <- counts$y
xte <- t(counts$xte + 1)
ctrl <- nbldaControl(folds = 2, repeats = 2)

fit <- trainNBLDA(x = x, y = y, type = "mle", tuneLength = 10,
                  metric = "accuracy", train.control = ctrl)

normalization(fit)

}
\seealso{
\code{\link{trainNBLDA}}
}
