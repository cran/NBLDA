
# NOTE: If setMethod is included in the same file with the corresponding function 'predict.nblda', then aliases
# are not required to be defined.

###### Show Methods  ########

#' @title Show Method for the S4 classes in NBLDA Package
#'
#' @description Pretty print the objects in S4 classes on R console.
#'
#' @param object an object of class \code{nblda, nblda_trained} and \code{nblda_input} to be printed.
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
#' show(fit)
#' show(inputs(fit))
#' show(nbldaTrained(fit))
#'
#' @name show
#' @rdname show
#'
#' @importFrom methods show
#' @method show nblda
show.nblda <- function(object){
  cat("\n", sep = " ")
  cat("  An object of class ", "\"", class(object), "\"","\n", sep = "")
  cat("  Model Description: Negative Binomial Linear Discriminant Analysis (NBLDA)", "\n\n", sep = "")

  tune.res <- object@result@crossValidated$tuning.results
  idx <- which(tune.res[ ,"rho"] == object@result@crossValidated$best.rho)
  acc <- tune.res[idx, "accuracy"]
  spars <- (ncol(object@input@x) - tune.res[idx, "nonzero"]) / ncol(object@input@x)

  cat("     Normalization : ", normalization(object), "\n\n")
  cat("       Accuracy(%) : ", sprintf("%6s", formatC(100*acc, digits = 2, format = "f")), "\n")
  cat("       Sparsity(%) : ", sprintf("%6s", formatC(100*spars, digits = 2, format = "f")), "\n\n")

  ### Notes:
  cat("-- NOTES -- \n")
  cat(" (1) Cross-validated accuracy is given. Run 'nbldaTrained()' for details. \n")
  cat(" (2) Sparsity is the percentage of removed variables from NBLDA model. \n")
  cat("     As sparsity increases, less variable included in the model. \n\n")
}

#' @rdname show
#'
#' @export
setMethod(f = "show", signature = signature(object = "nblda"), definition = show.nblda)


#' @rdname show
#'
#' @method show nblda_trained
show.nblda_trained <- function(object){
  n <- ceiling(length(object@finalModel$y))
  p <- ceiling(ncol(object@finalModel$x))
  nc <- ceiling(length(unique(object@finalModel$y)))
  ndigits <- length(unlist(strsplit(as.character(max(n, p, nc)), "")))
  class.names <- levels(object@finalModel$y)
  cNames <- paste(sapply(class.names, function(x)paste("\'", x, "\'", sep = "")), collapse = ", ")

  # PART 1. Model description and summary of sample/feature sizes. Class labels info.
  cat("\n", "Negative binomial linear discriminant analysis (NBLDA)", sep = "", "\n\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " samples", sep = ""), n), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " predictors", sep = ""), p), "\n")
  cat(sprintf(paste("%", ndigits + 1, "d", " classes: ", cNames, "  (Reference category: '", class.names[1], "')", sep = ""), nc), "\n\n")

  normalizationInfo <- if (normalization(object) == "none"){
    "Normalization is NOT performed."
  } else if (normalization(object) == "deseq"){
    "DESeq median ratio."
  } else if (normalization(object) == "mle"){
    "Maximum likelihood."
  } else if (normalization(object) == "quantile"){
    "Quantile."
  } else {
    "Trimmed-mean of M values."
  }

  # PART 2. Normalization and Resampling (Cross validation) info
  cat("Normalization: ", normalizationInfo, "\n", sep = "")
  cat("Power transformation is ", ifelse(object@control$transform, "performed.", "NOT performed."), "\n", sep = "")
  cat("Resampling: Cross-Validated (", object@control$folds, " folds, repeated ", object@control$repeats, " times)", sep = "", "\n")

  foldIdx <- object@control$foldIdx
  foldSampleSize <- n - unlist(lapply(unlist(foldIdx, recursive = FALSE), length))

  foldSampleSizeText <- if (length(foldSampleSize) > 5){
    paste(c(foldSampleSize[1:5], "..."), collapse = ", ", sep = "")
  } else {
    paste(c(foldSampleSize), collapse = ", ", sep = "")
  }

  cat("Summary of sample sizes: ", foldSampleSizeText, "\n", sep = "")

  selectedFeaturesInfo <- object@finalModel$selectedFeatures
  featureSelectionText <- if (is.null(selectedFeaturesInfo$names)){
    "All features are selected."
  } else {
    paste(length(object@finalModel$selectedFeatures$idx), " out of ", p, " features are selected.", sep = "")
  }
  cat("Summary of selected features: ", featureSelectionText, "\n")

  # PART 3. Tuning parameter info.
  # PLDA and PLDA2 methods have tuning parameter.
  tuningResults <- object@crossValidated$tuning.results

  cat("\n")
  cat(sprintf("%10s", "rho"), " ", sprintf("%8s", "Avg.Error"), " ", sprintf("%14s", "Avg.NonZeroFeat."), " ", sprintf("%8s", "Accuracy"),"\n")
  for (i in 1:nrow(tuningResults)){
    cat(sprintf("%10.5f", tuningResults[i, 1]), " ", sprintf("%9.2f", tuningResults[i, 3]), " ", sprintf("%15.2f", tuningResults[i, 4]), " ", sprintf("%9.4f", tuningResults[i, 2]), "\n")
  }

  # PART 4. Final notes on optimum model parameters.
  cat("\n")
  cat("The optimum model is obtained when rho = ",  sprintf("%6.5f", object@crossValidated$best.rho), " with an overall performance of","\n", sep = "")
  cat(ifelse(object@crossValidated$metric == "accuracy", "Accuracy = ", "Error = "),  sprintf("%6.4f", object@crossValidated$best.metric), " over folds. On the average ",
      object@crossValidated$best.nonzero, " out of ", p, " features was used", "\n", sep = "")
  cat("in the classifier." ,"\n\n", sep = "")

  cat("NOTES: The optimum model is selected using tuning parameter 'rho' which achieves \n")
  cat(" the lowest classification error (Avg.Error) or the highest Accuracy. The ", "\n")
  cat(" classification error is given as the average number of misclassified samples", "\n")
  cat(" over cross-validated folds. Similarly, the 'Avg.NonZeroFeat.' is the average ", "\n")
  cat(" number of non-zero features (the selected variables in the classification task)  ", "\n")
  cat(" over cross-validated folds. As the number of non-zero features decreases, the ", "\n")
  cat(" model becomes more sparse.", "\n\n")
}

#' @rdname show
#'
#' @export
setMethod(f = "show", signature = signature(object = "nblda_trained"), definition = show.nblda_trained)


#' @rdname show
#'
#' @method show nblda_input
show.nblda_input <- function(object){

  x <- object@x
  y <- object@y

  cat("x: Raw counts. A numeric ", class(x), ".", sep = "", "\n")
  cat(" dims(rows, columns): ", nrow(x), " by ", ncol(x), sep = "", "\n")

  rNames <- rownames(x)
  cNames <- colnames(x)

  rNamesText <- if (!is.null(rNames)){
    if (length(rNames) > 3){
      paste(paste(rNames[1:3], collapse = ", "), ", ...", sep = "")
    } else {
      paste(rNames, collapse = ", ")
    }
  }

  cNamesText <- if (!is.null(cNames)){
    if (length(cNames) > 3){
      paste(paste(cNames[1:3], collapse = ", "), ", ...", sep = "")
    } else {
      paste(cNames, collapse = ", ")
    }
  }

  cat(" rownames(", nrow(x), "): ", rNamesText, sep = "", "\n")
  cat(" colnames(", ncol(x), "): ", cNamesText, sep = "", "\n\n")
  cat("y: Class labels. A ", class(y), " vector.", sep = "", "\n")

  yText <- if (!is.null(y)){
    if (length(y) > 3){
      paste(paste(y[1:3], collapse = ", "), ", ...", sep = "")
    } else {
      paste(y, collapse = ", ")
    }
  }
  cat(" classes(", length(y), "): ", yText, sep = "", "\n")

  yLevels <- levels(y)

  yLevelsText <- if (length(yLevels) > 3){
    paste(paste(yLevels[1:3], collapse = ", "), ", ...", sep = "")
  } else {
    paste(yLevels, collapse = ", ")
  }
  cat(" levels(", length(levels(y)), "): ", yLevelsText, sep = "", "\n\n")

}

#' @rdname show
#'
#' @export
setMethod(f = "show", signature = signature(object = "nblda_input"), definition = show.nblda_input)

###### Getter Functions #########

##### type (normalization) #######
#' @title Accessors for the 'type' slot.
#'
#' @description This slot stores the name of normalization method. Normalization is defined using \code{type} argument in \code{\link{trainNBLDA}} function.
#'
#' @docType methods
#' @name normalization
#' @rdname normalization
#' @aliases normalization,nblda-method
#'
#' @param object an \code{nblda} or \code{nblda_trained} object.
#'
#' @seealso \code{\link{trainNBLDA}}
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
#' normalization(fit)
#'
#'@export
setMethod(f = "normalization", signature = signature(object = "nblda"), function(object){
  object@result@finalModel$type
})


#' @rdname normalization
#' @aliases normalization,nblda_trained-method
#'
#' @export
setMethod(f = "normalization", signature = signature(object = "nblda_trained"), function(object){
  object@finalModel$type
})



###### control ######
#' @title Accessors for the 'control' slot.
#'
#' @description This slot stores control parameters for training NBLDA model.
#'
#' @docType methods
#' @name control
#' @rdname control
#' @aliases control,nblda-method
#' @param object an \code{nblda} or \code{nblda_trained} object.
#'
#' @seealso \code{\link{trainNBLDA}}
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
#' control(fit)
#'
#'@export
setMethod(f = "control", signature = signature(object = "nblda"), function(object){
  object@result@control
})


#' @rdname control
#' @aliases control,nblda_trained-method
setMethod(f = "control", signature = signature(object = "nblda_trained"), function(object){
  object@control
})



###### nbldaTrained ######
#' @title Accessors for the 'crossValidated' slot.
#'
#' @description This slot stores the results for cross-validated model, e.g tuning results, optimum model parameters etc.
#'
#' @docType methods
#' @name nbldaTrained
#' @rdname nbldaTrained
#' @aliases nbldaTrained,nblda-method
#' @param object an \code{nblda} or \code{nblda_trained} object.
#'
#' @seealso \code{\link{trainNBLDA}}
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
#' nbldaTrained(fit)
#'
#'@export
setMethod(f = "nbldaTrained", signature = signature(object = "nblda"), function(object){
  object@result
})

#' @rdname nbldaTrained
#' @aliases nbldaTrained,nblda_trained-method
setMethod(f = "nbldaTrained", signature = signature(object = "nblda_trained"), function(object){
  object
})


###### inputs ######
#' @title Accessors for the 'input' slot.
#'
#' @description This slot stores the input data for trained model.
#'
#' @docType methods
#' @name inputs
#' @rdname inputs
#' @aliases inputs,nblda-method
#'
#' @param object an \code{nblda} object.
#'
#' @seealso \code{\link{trainNBLDA}}
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
#' inputs(fit)
#'
#'@export
setMethod(f = "inputs", signature = signature(object = "nblda"), function(object){
  object@input
})

###### selectedFeatures ######
#' @title Accessors for the 'selectedFeatures' slot.
#'
#' @description This slot, if not NULL, stores the selected features/variables for sparse model.
#'
#' @docType methods
#' @name selectedFeatures
#' @rdname selectedFeatures
#' @aliases selectedFeatures,nblda-method
#'
#' @param object an \code{nblda} or \code{nblda_trained} object.
#'
#' @return a list of selected features info including the followings:
#' \item{idx}{column indices of selected features/variables}
#' \item{names}{column names of selected features/variables if input data have pre-defined column names.}
#'
#' @note If \code{return.selected.features} = FALSE within \code{\link{nbldaControl}} or all features/variables
#' are selected and used in discrimination function, \code{idx} and \code{names} are returned NULL.
#'
#' @seealso \code{\link{trainNBLDA}}, \code{\linkS4class{nblda}}, \code{\linkS4class{nblda_trained}}
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 50, K = 2, param = 1, sdsignal = 0.5, DE = 0.6,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#' x <- t(counts$x + 1)
#' y <- counts$y
#' xte <- t(counts$xte + 1)
#' ctrl <- nbldaControl(folds = 2, repeats = 2, return.selected.features = TRUE,
#'                      transform = TRUE, phi.epsilon = 0.10)
#'
#' fit <- trainNBLDA(x = x, y = y, type = "mle", tuneLength = 10,
#'                   metric = "accuracy", train.control = ctrl)
#'
#' selectedFeatures(fit)
#'
#'@export
setMethod(f = "selectedFeatures", signature = signature(object = "nblda"), function(object){
  object@result@finalModel$selectedFeatures
})

#' @rdname selectedFeatures
#' @aliases selectedFeatures,nblda_trained-method
setMethod(f = "selectedFeatures", signature = signature(object = "nblda_trained"), function(object){
  object@finalModel$selectedFeatures
})


##### Plots ######
#' @rdname plot
#' @aliases plot,nblda-method
#'
#' @export
setMethod(f = "plot", signature = signature(x = "nblda"), definition = plot.nblda)

#' @rdname plot
#' @aliases plot,nblda_trained-method
#'
#' @export
setMethod(f = "plot", signature = signature(x = "nblda_trained"), definition = plot.nblda_trained)

##### Predictions ######
#' @rdname predict
#' @aliases predict,nblda-method
#'
#' @export
setMethod(f = "predict", signature = signature(object = "nblda"), definition = predict.nblda)
