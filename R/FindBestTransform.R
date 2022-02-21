#' @title Find the Power Transformation Parameter.
#'
#' @description Use this function to find a constant value of alpha to be used for transforming count data. The power transformation parameter \code{alpha}, which approximately fits transformed data to the Poisson log-linear model, is selected using a grid search within the interval [0, 1].
#'
#' @param x an n-by-p data frame or matrix of count data. Samples should be in the rows.
#' @param grid.length how many distinct points of alpha should be searched within the interval [0, 1]? Default is 50.
#'
#' @return the value of alpha to be used within the power transformation.
#'
#' @author Dincer Goksuluk
#'
#' @note This function is copied from \code{PoiClaClu} package and modified to control the total number of grid search.
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#'
#' x <- counts$x
#' FindBestTransform(x)
#'
#' @name FindBestTransform
#' @rdname FindBestTransform
#'
#' @export
FindBestTransform <- function(x, grid.length = 50){
  alphas <- seq(.01, 1, len = grid.length)
  gof <- rep(NA, length(alphas))
  for (alpha in alphas){
    gof[alphas == alpha] <- GoodnessOfFit(x^alpha, type = "mle")
  }
  return(alphas[which.min(abs(gof - (nrow(x) - 1)*(ncol(x) - 1)))])
}
