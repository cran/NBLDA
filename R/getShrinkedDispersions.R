#' @title Estimate Shrinked Overdispersions
#'
#' @description Use this function to shrink initial estimates of overdispersions towards a target value.
#'
#' @param obs a numeric vector. Initial dispersion estimates for each feature.
#' @param shrinkTarget a numeric value. Initial dispersion estimates are shrunk towards this value. If NULL, target value is estimated from the initial dispersion estimates. See notes.
#' @param delta a numeric value. This is the weight that is used within the shrinkage algorithm. If 0, no shrinkage is performed on initial values. If equals 1, initial values are forced to be shrunken to the target value. If NULL, weights are automatically estimated from the initial dispersion estimates.
#'
#' @note This function is modified from the source codes of \code{\link[sSeq]{getAdjustDisp}} function in the \bold{sSeq} Bioconductor package.
#'
#' @return a list with the elements of initial and adjusted (shrunken) dispersion estimates, shrinkage target, and
#' weights that are used to shrink towards the target value. See the related paper for detailed information on shrinkage
#' algorithm (Yu et al., 2013).
#' \item{initial}{initial dispersion estimates using the method-of-moments.}
#' \item{adj}{shrunken dispersion estimates.}
#' \item{cmp}{the means and variances of initial estimates.}
#' \item{delta}{a weight used for shrinkage estimates. See Yu et al. (2013) for details.}
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

  return(
    list(
      initial = obs,
      adj = disp.shrinked,
      cmp = cmp,
      delta = delta,
      target = shrinkTarget
    )
  )
}
