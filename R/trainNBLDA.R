nblda <- function(x, y, xte = NULL, rhos = 0, beta = 1, type = c("mle", "deseq", "quantile", "tmm"),
                  prior = NULL, transform = FALSE, alpha = NULL, truephi = NULL, phi.epsilon = 0.05,
                  target = NULL, normalize.target = TRUE, delta = NULL, return.selected.features = FALSE, ...){

  # xte = NULL; rhos = best.rho; beta = tc$beta; type = type; prior = tc$prior;
  # transform = FALSE; alpha = alpha; return.selected.features = tc$return.selected.features;
  # truephi = tc$truephi; target = tc$target; normalize.target = tc$normalize.target;
  # delta = tc$delta; phi.epsilon = tc$phi.epsilon

  type <- match.arg(type)

  if (!is.null(target)){
    if (target < 0){
      warning("'target' should be non-negative. It is set to 0.")
      target <- 0
    }
  }

  if (is.null(phi.epsilon)){
    phi.epsilon <- 0
  }

  if (is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }

  if (is.null(rhos)){
    warning("'rhos' cannot be NULL. Setting rho = 0.")
  }

  # if (length(rhos) > 1){
  #   warning("Multiple values of rho is given. Only the first element of rho is used.")
  # }

  if (transform){
    if (is.null(alpha)){
      ### Train ve test seti buradaki "alpha" değerine göre transform ediliyor.
      ### Alpha değeri MLE altında size factor hesabı kullanılarak yapılıyor.
      ### TO DO ####
      ### Further Research: Farklı size factor kestirimleri altında power transformations'lar kullanılabilir.
      alpha <- FindBestTransform(x)
    }

    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }

  if (length(unique(y)) < 2){
    stop("Response variable must have at least two classes.")
  }

  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y)))
  }

  null.out <- NullModel(x, type = type)
  ns <- null.out$n
  nste <- NullModelTest(null.out, xte)$nste

  uniq <- sort(unique(y))

  ## Dispersion estimates.
  if (!is.null(truephi)){
    if (length(truephi) >= 2 & (length(truephi) != ncol(ns))){
      truephi <- truephi[1]
      disperhat <- rep(truephi, ncol(ns))
      warning("The length of \"truephi\" should be the same as number of features. Only the first element is used and replicated for all features.")

    } else if (length(truephi) == 1){
      disperhat <- rep(truephi, ncol(ns))

    } else if (length(truephi) >= 2 & (length(truephi) == ncol(ns))){
      disperhat <- truephi
    }
    disperhat.list = list(initial = disperhat, adj = disperhat, cmp = NULL, delta = NULL, target = NULL)

  } else {
    ####### TO DO #######
    # tt değerinin 0'a eşit olması durumu ile analiz yapılıyor ama normalde üstteki değere shrink ediliyor. Yu et al çalışmasındaki gibi.
    # Bunun nedeni araştırılacak. 0 Değerini almak problem yaratır mı diye incelenecek.
    if (is.null(target)){
      if (normalize.target){
        tt <- getT(t(x), sizeFactors = null.out$rawsizestr, verbose = FALSE, propForSigma = c(0, 1))$target
      } else {
        tt <- getT(t(x), sizeFactors = rep(1, nrow(x)), verbose = FALSE, propForSigma = c(0, 1))$target
      }

    } else {
      tt <- target
    }

    ### Moment estimation of gene-wise dispersions.
    rM = rowMeans(t(x))
    rV = rowVars(t(x))

    disp = (rV - rM)/rM^2   ## Dispersion estimates using method-of-moments.
    # disp0 = numeric()

    ## Negative dispersions are set to 0.
    disp0 <- sapply(disp, function(x){
      max(0, x)
    })

    # disperhat = getAdjustDisp(disp0, shrinkTarget = tt, verbose = FALSE)$adj
    disperhat.list = getShrinkedDispersions(obs = disp0, shrinkTarget = tt, delta = delta)
    disperhat <- disperhat.list$adj
    disperhat <- pmax(rep(0, length(disperhat)), disperhat) ## Negative dispersions are set to 0.
  }

  phihat <- as.numeric(disperhat)

  ### Shrink dispersion estimates to 0 using phi.epsilon.
  if (any(phihat <= phi.epsilon)){
    phihat[phihat <= phi.epsilon] <- 0
    disperhat <- disperhat.list$adj <- phihat
  }

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  save <- lapply(rhos, function(u){
    ds <- GetD(ns = ns, x = x, y = y, rho = u, beta = beta)
    discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

    # NBLDA
    for (k in 1:length(uniq)){
      dstar = ds[k, ]
      part2 = 1 + t(nste) * dstar * phihat
      part1 = dstar/part2

      part3 <- inv.phihat * log(part2)
      if (!is.null(id.inf)){
        part3.limits <- t(nste[ ,id.inf])*dstar[id.inf]
        part3[id.inf, ] <- part3.limits  # Replace zero dispersed genes with limit values.
      }

      discriminant[ ,k] <- rowSums(xte * t(log(part1))) - colSums(part3) + log(prior[k])
    }

    # Feature selection
    selectedFeaturesIdx.ds <- apply(ds, 2, function(x){
      all(x != 1)
    })

    selectedFeaturesIdx.disperhat <- disperhat != 0
    selectedFeaturesIdx <- which(selectedFeaturesIdx.ds | selectedFeaturesIdx.disperhat)
    selectedFeaturesNames <- colnames(x)[selectedFeaturesIdx]

    ## If all features are selected, return NULL.
    if (length(selectedFeaturesNames) == ncol(x)){
      selectedFeaturesIdx <- selectedFeaturesNames <- NULL
    }

    save.i <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                   alpha = alpha, rho = u, x = x, y = y, xte = xte, type = type, prior = prior, dispersions = disperhat.list,
                   transform = transform, trainParameters = null.out,
                   selectedFeatures = list(idx = selectedFeaturesIdx, names = selectedFeaturesNames))
    save.i
  })

  # class(save) <- "poicla"
  return(save)
}



#' @title Train Model over Different Tuning Parameters
#'
#' @description This function fits Negative Binomial classifier using various model parameters and finds the best model parameter
#' using the resampling based performance measures.
#'
#' @usage
#' trainNBLDA(x, y, type = c("mle", "deseq", "quantile", "tmm"),
#'   tuneLength = 10, metric = c("accuracy", "error"), train.control = nbldaControl(), ...)
#'
#' @param x a n-by-p data frame or matrix. Samples should be in the rows and variables in the columns. Used to train the classifier.
#' @param y a vector of length n. Each element corresponds to a class label of a sample. Integer and/or factor types are allowed.
#' @param type a character string indicating the type of normalization method within the NBLDA model. See details.
#' @param tuneLength a positive integer. This is the total number of levels to be used while tuning the model parameter(s).
#' @param metric which criteria should be used while determining the best parameter? overall accuracu or avarage number of misclassified samples?
#' @param train.control a list with control parameters to be used in NBLDA model. See \link{nbldaControl} for details.
#' @param ... further arguments. Deprecated.
#'
#' @return an \code{nblda} object with following slots:
#' \item{input}{an \code{nblda_input} object including the raw count data and response variable. See \code{\linkS4class{nblda_input}} for details.}
#' \item{result}{an \code{nblda_trained} object including the results from cross-validated and final models. See \code{\linkS4class{nblda_trained}} for details.}
#' \item{call}{a call expression.}
#'
#' @details NBLDA is proposed to classify count data from any field, e.g. economics, social sciences, genomics, etc.
#' In RNA-Seq studies, for example, normalization is used to adjust between-sample differences for downstream analysis.
#' \code{type} is used to define normalization method. Available options are "mle", "deseq", "quantile" and "tmm".
#' Since "deseq", "quantile" and "tmm" methods are originally proposed as robust methods to be used in RNA-Sequencing studies, one should
#' carefully define normalization types. In greater details, "deseq" estimates the size factors by dividing each sample by the geometric means
#' of the transcript counts (Anders and Huber, 2010). "tmm" trims the lower and upper side of the data by log fold changes to
#' minimize the log-fold changes between the samples and by absolute intensity (Robinson and Oshlack, 2010). "quantile" is
#' quantile normalization approach of Bullard et al (2010). "mle" (less robust) divides total counts of each sample to the grand total
#' counts (Witten, 2010). See related papers for mathematical backgrounds.
#'
#' @author Dincer Goksuluk
#'
#' @references Witten, DM (2011). Classification and clustering of sequencing data using a Poisson model.
#' Ann. Appl. Stat. 5(4), 2493--2518. doi:10.1214/11-AOAS493.
#'
#' Dong, K., Zhao, H., Tong, T., & Wan, X. (2016). NBLDA: negative binomial linear discriminant analysis for RNA-Seq data.
#' BMC Bioinformatics, 17(1), 369. http://doi.org/10.1186/s12859-016-1208-1.
#'
#' Anders S. Huber W. (2010). Differential expression analysis for sequence count data. Genome Biology, 11:R106
#'
#' Witten D. et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis
#' in a collection of cervical tumours and matched controls. BMC Biology, 8:58
#'
#' Robinson MD, Oshlack A (2010). A scaling normalization method for differential expression analysis of RNA-Seq data.
#' Genome Biology, 11:R25, doi:10.1186/gb-2010-11-3-r25
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
#' fit
#' nbldaTrained(fit)  # Cross-validated model summary.
#'
#' @name trainNBLDA
#' @rdname trainNBLDA
#'
#' @importFrom methods new
#'
#' @export
trainNBLDA <- function(x, y, type = c("mle", "deseq", "quantile", "tmm"), tuneLength = 10, metric = c("accuracy", "error"),
                       train.control = nbldaControl(), ...){
  # type = "deseq"
  # tuneLength = 10
  # metric = "accuracy"
  # train.control = ctrl2

  call <- match.call()
  type <- match.arg(type)
  metric <- match.arg(metric)
  y <- as.factor(y)
  tc <- train.control

  if (tc$transform){
    if (is.null(tc$alpha)){
      alpha <- FindBestTransform(x)
    }
    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    tc$alpha <- alpha
  } else {
    alpha <- NULL
  }

  # 1. Define sets of model parameter values to evaluate
  if (is.null(tc$rhos)){
    ns <- NullModel(x, type = type)$n
    uniq <- sort(unique(y))
    maxrho <- rep(NA, length(uniq))

    for (k in 1:length(uniq)){
      a <- colSums(x[y == uniq[k],]) + tc$beta
      b <- colSums(ns[y == uniq[k],]) + tc$beta
      maxrho[k] <- max(abs(a/b - 1) * sqrt(b), na.rm = TRUE)
    }

    tc$rhos <- rhos <- seq(0, max(maxrho, na.rm = TRUE) * (2/3), len = tuneLength)
  } else {
    rhos <- tc$rhos
  }

  # 2. For each parameter set do
  if (is.null(tc$foldIdx)){
    # Each element of the returned list is the indices of test samples in this fold.
    allFolds <- lapply(1:tc$repeats, function(x){
      tmp <- balanced.folds(y, nfolds = tc$folds)
      names(tmp) <- paste("Fold.", 1:tc$folds, sep = "")
      tmp
    })
    names(allFolds) <- paste("Repeat.", 1:tc$repeats, sep = "")

    tc$foldIdx <- allFolds
    repeats <- tc$repeats
    folds <- tc$folds

  } else {
    allFolds <- tc$foldIdx
    repeats <- length(allFolds)
    folds <- length(allFolds[[1]])

    tc$repeats <- repeats
    tc$folds <- folds
  }

  tuningRes <- tuningResults <- NULL
  tuneGrid <- expand.grid(rho = rhos, epsilon = tc$phi.epsilon)

  # for each fold within each repat, fit a plda model.
  # tuningResults <- foreach(iii = 1:nrow(tuneGrid), .combine = "c", .inorder = TRUE, .verbose = FALSE,
  #                          .export = c("tuneGrid", "allFolds")) %dopar% {
  #
  #
  # }

  foldRes <- lapply(allFolds, function(ii){
    errors <- acc <- nonzeros <- matrix(NA, nrow = length(rhos), ncol = length(ii))

    for (jj in 1:length(ii)){
      idx <- ii[[jj]]
      # out <- plda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = beta,
      #             type = "quantile", prior = prior, transform = FALSE)

      plist <- tc
      plist$transform <- FALSE
      plist <- c(list(x = x[-idx, ], y = y[-idx], xte = x[idx, ], type = type), plist)

      out <- do.call("nblda", plist)
      # out <- nblda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = tc$beta, type = type,
      #              prior = tc$prior, transform = FALSE, alpha = alpha, truephi = tc$true.dispersions,
      #              target = tc$target, normalize.target = tc$normalize.target, delta = tc$delta,
      #              return.selected.features = tc$return.selected.features)

      err.i <- nonzero.i <- acc.i <- NULL
      for (i in 1:length(out)){
        acc.i <- c(acc.i, sum(out[[i]][["ytehat"]] == y[idx]) / length(y[idx]))
        err.i <- c(err.i, sum(out[[i]][["ytehat"]] != y[idx]))
        nonzero.i <- c(nonzero.i, sum(colSums(out[[i]][["ds"]] != 1) != 0 | out[[i]]$dispersions$adj != 0))
      }

      acc[ ,jj] <- acc.i
      errors[ ,jj] <- err.i
      nonzeros[ ,jj] <- nonzero.i
    }

    list(errors = rowMeans(errors), accuracies = rowMeans(acc), nonzeros = rowMeans(nonzeros))
    # list(errors = errors, nonzeros = nonzeros)
  })

  err.sum <- nonzero.sum <- acc.sum <- 0
  for (r in 1:length(foldRes)){
    acc.sum <- acc.sum + foldRes[[r]]$accuracies
    err.sum <- err.sum + foldRes[[r]]$errors
    nonzero.sum <- nonzero.sum + foldRes[[r]]$nonzeros
  }

  tuningRes <- cbind(rhos, round(acc.sum/repeats, 10), round(err.sum/repeats, 10), round(nonzero.sum/repeats, 10))
  colnames(tuningRes) <- c("rho", "accuracy", "errors", "nonzero")

  if (metric == "accuracy"){
    best.idx <- max(which(tuningRes[, "accuracy"] == max(tuningRes[ ,"accuracy"])))
    best.rho <- tuningRes[, "rho"][best.idx]
    best.metric <- tuningRes[, "accuracy"][best.idx]
    best.nonzero <- tuningRes[, "nonzero"][best.idx]
  } else if (metric == "error"){
    best.idx <- max(which(tuningRes[, "errors"] == min(tuningRes[ ,"errors"])))
    best.rho <- tuningRes[, "rho"][best.idx]
    best.metric <- tuningRes[, "errors"][best.idx]
    best.nonzero <- tuningRes[, "nonzero"][best.idx]
  }

  # Set transform to FALSE since data is already transformed above.
  final.model <- nblda(x = x, y = y, xte = NULL, rhos = best.rho, beta = tc$beta, type = type, prior = tc$prior,
                       transform = FALSE, alpha = alpha, return.selected.features = tc$return.selected.features,
                       truephi = tc$truephi, target = tc$target, normalize.target = tc$normalize.target,
                       delta = tc$delta, phi.epsilon = tc$phi.epsilon)

  final.model[[1]]$transform <- tc$transform
  res <- list(tuning.results = tuningRes, best.rho = best.rho, final.model = final.model[[1]])

  ## Input object
  input.list <- new("nblda_input", x = x, y = y)  # Input data.

  # Trained model object
  result.list <- new("nblda_trained",
                     crossValidated = list(tuning.results = tuningRes, best.rho = best.rho,
                                           best.metric = best.metric, best.nonzero = best.nonzero, metric = metric),
                     finalModel = final.model[[1]],
                     control = tc)

  ## nblda object.
  new("nblda", input = input.list,
      result = result.list,
      call = call)
}
