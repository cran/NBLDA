
#' @title Extract predictions from NBLDA model
#'
#' @description  This function predicts the class labels of test data for a given model.
#'
#' @param object a \code{nblda} object returned from \code{\link{trainNBLDA}}.
#' @param test.data a data frame or matrix whose class labels to be predicted.
#' @param return what should be returned? Predicted class labels or eveything?
#' @param \dots further arguments to be passed to or from methods.
#'
#' @return It is possible to return only predicted class labels or a list with elements which are used within prediction
#' process. These arguements are as follows:
#' \item{xte}{count data for test set.}
#' \item{nste}{normalized count data for test set.}
#' \item{ds}{estimates of offset parameter for each variable. See notes.}
#' \item{discriminant}{discriminant scores of each subject.}
#' \item{prior}{prior probabilities for each class.}
#' \item{ytehat}{predicted class labels for test set.}
#' \item{alpha}{power transformation parameter. If no transformation is requested, it returns NULL.}
#' \item{type}{normalization method.}
#' \item{dispersions}{dispersion estimates of each variable.}
#'
#' @note \code{d_kj} is simply used to re-parameterize the Negative Binomial mean as s_i*g_j*d_kj where s_i is the size
#' factor for subject i, g_j is the total count of variable j and d_kj is the offset parameter for variable j at class k.
#'
#' @author Dincer Goksuluk
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#' x <- t(counts$x + 1)
#' y <- counts$y
#' xte <- t(counts$xte + 1)
#' ctrl <- nbldaControl(folds = 2, repeats = 2)
#'
#' fit <- trainNBLDA(x = x, y = y, type = "mle", tuneLength = 10,
#'                   metric = "accuracy", train.control = ctrl)
#'
#' predict(fit, xte)
#'
#' @name predict
#' @rdname predict
#'
#' @importFrom stats predict
#' @method predict nblda
#'
#' @export
predict.nblda <- function(object, test.data, return = c("predictions", "everything"), ...){
  ## Args:
  ##  object: an object returned form trainNBLDA.
  ##  test.data: a data.frame or matrix for test set. Samples in the rows and features in the columns.
  ##  return: should return everything or only predictions?

  return <- match.arg(return)

  fitted.model <- object@result@finalModel
  xte <- test.data
  type <- fitted.model$type

  if (fitted.model$transform){
    alpha <- fitted.model$alpha
    if (alpha <= 0 || alpha > 1){
      stop("alpha must be between 0 and 1")
    }
    xte <- xte^alpha
  } else {
    alpha <- NULL
  }

  ### prior trained model'den alınacak.
  prior <- fitted.model$prior

  # null.out <- NullModel(x, type = type)  ### trained modelden alınacak
  # ns <- null.out$n
  nste <- NullModelTest(fitted.model$trainParameters, xte = xte)$nste

  uniq <- sort(unique(fitted.model$y))
  # ds <- GetD(ns, x, y, rho, beta)   ### trained modelden alınacak
  ds <- fitted.model$ds
  phihat <- fitted.model$dispersions$adj
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  for (k in 1:length(uniq)){
    dstar = ds[k, ]
    part2 = 1 + t(nste) * dstar * phihat
    part1 = dstar/part2

    part3 <- inv.phihat * log(part2)
    if (!is.null(id.inf)){
      part3.limits <- t(nste[ ,id.inf]) * dstar[id.inf]
      part3[id.inf, ] <- part3.limits  # Replace zero dispersed genes with limit values.
    }

    discriminant[ ,k] <- rowSums(xte * t(log(part1))) - colSums(part3) + log(prior[k])
  }

  pred <- uniq[apply(discriminant, 1, which.max)]

  if (return == "predictions"){
    return(pred)
  } else if (return == "everything"){
    list(xte = xte, nste = nste, ds = ds, discriminant = discriminant, prior = prior,
         ytehat = pred, alpha = alpha, type = type, dispersions = phihat)
  }
}
