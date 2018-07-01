#' @title Estimate Shrinked Overdispersions
#'
#' @description Use this function to shrink initial estimates of overdispersions towards a target value.
#'
#' @param obs a numeric vector. Initial dispersion estimates for each feature.
#' @param shrinkTarget a numeric value. initial dispersion estimates are shrinked towards
#' this value. If NULL, target value is estimated from initial dispersion estimates. See notes.
#' @param delta a numeric value. This is the weight that is used for shrinkage algorithm. If 0,
#' no shrinkage is performed on intiial values. If equals 1, initial values are forced to shrinked
#' to target value. If NULL, weight are automatically estimated from initial disperson estimates.
#'
#' @note This function is modified using source code from \code{\link[sSeq]{getAdjustDisp}}.
#'
#' @return a list with elements of initial and adjusted (shrinked) dispersion estimates, shrinkage target and
#' weight that is used to shrink towards target value. See related paper for detailed information on shrinkage
#' algorithm (Yu et. al., 2013).
#' \item{initial}{initial estimates for dispersions estimated from method-of-momnets.}
#' \item{adj}{shrinked dispersion estimates.}
#' \item{cmp}{mean and variance of initial estimates.}
#' \item{delta}{a weight used for shrinkage estimates. See Yu et. al. (2013) for details.}
#' \item{target}{shrinkage target for initial dispersion estimates.}
#'
#' @author Dincer Goksuluk
#'
#' @references Yu, D., Huber, W., & Vitek, O. (2013). Shrinkage estimation of dispersion in Negative Binomial models
#' for RNA-seq experiments with small sample size. Bioinformatics, 29(10), 1275-1282.
#'
#' @seealso \code{\link[sSeq]{getT}}, \code{\link[sSeq]{getAdjustDisp}}
#'
#' @examples
#' set.seed(2128)
#' initial <- runif(10, 0, 4)
#'
#' getShrinkedDispersions(initial, 0)  # shrink towards 0.
#' getShrinkedDispersions(initial, 0, delta = 1)  # force to shrink 0.
#'
#' @name getShrinkedDispersions
#' @rdname getShrinkedDispersions
#'
#' @importFrom stats var
#'
#' @export
getShrinkedDispersions <- function(obs, shrinkTarget = NULL, delta = NULL){
  # This function returns the shrinked dispersion estimates.
  # Args:
  #   obs: a numeric vector. dispersion estimates for each feature.
  #   shrinkTarket: a numeric value. initial dispersion estimates are shrinked towards this value. If NULL, target value
  #                 is estimated from initial dispersion estimates.
  #   delta: a numeric value. This is the weight that is used for shrinkage algorithm. If 0, no shrinkage is
  #          performed on intiial values. If equals 1, initial values are forced to shrinked to target value.
  #          If NULL, weight are automatically estimated from initial disperson estimates.

  if (!is.null(delta)){
    if (delta < 0 | delta > 1){
      stop(warning("'delta' should be within [0, 1]. Calculation stopped..."))
    }
  }

  obs[is.na(obs)] <- 0
  upBound <- shrinkTarget

  # if (verbose) {
  #   print(paste("Shrink initial dispersion estimates toward ", shrinkTarget, ".", sep = ""))
  # }

  S.mat <- var(obs, na.rm = T)
  cmp <- data.frame(mean = mean(obs, na.rm = T), sigmasq.part = S.mat)
  mean.mat <- rep(upBound, length(obs))
  dif.mat <- obs - mean.mat
  dif2.mat <- sum(dif.mat^2)

  if (is.null(delta)){
    delta <- ((length(obs) - 2) * S.mat/(dif2.mat))   ## 1 - delta
    delta <- pmax(0, delta)
  }

  disp.shrinked <- delta * mean.mat + (1 - delta) * obs

  return(list(initial = obs, adj = disp.shrinked, cmp = cmp, delta = delta, target = shrinkTarget))
}
