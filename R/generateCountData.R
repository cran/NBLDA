#' @title Generate Count Data
#'
#' @description This function can be used to generate counts, e.g RNA-Sequencing data, for both classification and clustering purposes.
#'
#' @usage
#' generateCountData(n, p, K, param, sdsignal = 1, DE = 0.3, allZero.rm = TRUE,
#'   tag.samples = FALSE)
#'
#' @param n number of samples.
#' @param p number of variables/features.
#' @param K number of classes.
#' @param param overdispersion parameter. This parameter is matched with the arguement \code{size} in \code{\link[stats]{rnbinom}} function.
#' Hence, Negative Binomial distribution aproximates to Poisson distribution as \code{param} increases.
#' @param sdsignal a nonzero numeric value. As \code{sdsignal} increases, the observed counts greatly differs among K classes.
#' @param DE a numeric value within the interval [0, 1]. This is the proportion of total number of variables that is significantly
#' different among K classes. The remaining part is assumed to be having no contribution to discrimination function.
#' @param allZero.rm a logical. If TRUE, columns having all zero cells are dropped.
#' @param tag.samples a logical. If TRUE, rownames are automatically generated. A tag for each sample such as "S1", "S2", etc.
#'
#' @return
#' \item{\code{x, xte}}{count data matrix for training and test set.}
#' \item{\code{y, yte}}{class labels for training and test set.}
#' \item{\code{truesf, truesfte}}{true size factors for training and test set. See Witten (2011) for more information on estimating size factors.}
#'
#' @author Dincer Goksuluk
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#' head(counts$x)
#'
#' @name generateCountData
#' @rdname generateCountData
#'
#' @importFrom stats rexp rnorm rnbinom runif sd var
#'
#' @export
generateCountData <- function(n, p, K, param, sdsignal = 1, DE = 0.3, allZero.rm = TRUE, tag.samples = FALSE) {
  # Args:
  #   n: number of samples
  #   p: number of features
  #   K: number of classes
  #   sdsignal: The extent to which the classes are different. If this equals zero then there
  #             are no class differences and if this is large then the classes are very different.
  #             This parameter is used to siulate dkj parameters.
  #   DE: probability of a gene being differentially expressed between classes.
  #   allZero.rm: a logical. If TRUE, features with only zero counts will be removed. Otherwise, it will be kept
  #               in the model as a presence of excess zeros. Might be suitable for zero-inflated scenarios.

  if (n < 4 * K){
    stop("We require n to be at least 4*K.")
  }
  if (DE < 0 || DE > 1){
    stop("DE should be a value within [0, 1].")
  }
  q0 <- rexp(p, rate = 1/25)   # Gene Total
  isDE <- runif(p) < DE
  classk <- matrix(NA, nrow = K, ncol = p)
  for (k in 1:K) {
    lfc <- rnorm(p, sd = sdsignal)   # dkj parametreleri.
    classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)   ## Hangi genlerin DE olduğu bilgisi. Bunları da generate edilen veri ile döndürebiliriz.
    ##
  }
  truesf <- runif(n) * 2 + 0.2
  truesfte <- runif(n) * 2 + 0.2
  conds <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  condste <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  x <- xte <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    for (k in 1:K) {
      if (conds[i] == k)
        x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ], size = param)
      if (condste[i] == k)
        xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ], size = param)
    }
  }
  if (allZero.rm){
    rm <- apply(x, 2, sum) == 0   # Removes all zero columns.
  } else {
    rm <- logical(ncol(x))
  }

  # x <- as.data.frame(x)
  # xte <- as.data.frame(xte)
  colnames(x) <- colnames(xte) <- paste("G", 1:p, sep = "")

  if (tag.samples){
    rownames(x) <- rownames(xte) <- paste("S", 1:n, sep = "")
  }

  return(structure(list(x = t(x[, !rm]), xte = t(xte[, !rm]), y = conds, yte = condste,
                        truesf = truesf, truesfte = truesfte), class = "count.data"))
}
