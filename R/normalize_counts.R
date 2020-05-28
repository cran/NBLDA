#' @title Calculate Normalized Counts and Related Training Parameters.
#'
#' @description Fit training set to NBLDA model and estimate normalized counts. The related model parameters which are used
#' while normalizing training sets are also returned in order to normalize test sets using training set parameters.
#'
#' @param x a n-by-p data frame or matrix of count data. Samples should be in the rows.
#' @param type normalization methods. See \code{\link{control}} for details.
#'
#' @return a list with normalized counts and training set parameters used for normalizing raw counts.
#'
#' @author Dincer Goksuluk
#'
#' @note These functions are copied from \code{PoiClaClu} package and modified here to make "tmm" and "none" methods available.
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#' x <- counts$x
#' xte <- counts$xte
#'
#' x.out <- NullModel(x, "mle")
#' x.out$n ## Normalized counts using "mle" method
#'
#' xte.out <- NullModelTest(x.out, xte)
#' xte.out$n  # Normalized counts for test set using train set parameters.
#'
#' @name NullModel
#' @rdname NullModel
#'
#' @export
NullModel <- function(x, type = c("mle", "deseq", "quantile", "none", "tmm")){
  # x MUST be a n times p matrix - i.e. observations on the rows and features on the columns
  type <- match.arg(type)
  rowsum <- rowSums(x)
  colsum <- colSums(x)

  if (type == "mle"){
    sizes <- rowSums(x)/sum(x)  # size factors for each sample
    mle <- outer(sizes, colsum, "*")
    return(list(n = mle, sizes = sizes, lib.size = rowsum, gene.length = colsum, type = type))

  } else if (type == "quantile"){ ### BU KISIM DÜZENLENECEK.
    # This is quantile normalization idea of Bullard et al 2010 -- quantile-normalize using 75th quantile of observed counts for each sample, excluding zero-counts
    sample.qts <- apply(x, 1, quantile, .75)
    sample.qts.pmax <- pmax(sample.qts, 1) # Don't wait to standardize by 0... min allowed is 1
    sizes <- sample.qts.pmax/sum(sample.qts.pmax)   # size factors for each sample
    fit <- outer(sizes, colsum, "*")
    return(list(n = fit, sizes = sizes, lib.size = rowsum, gene.length = colsum, sample.qts.tr = sample.qts,
                rawsizetr = sample.qts.pmax, type = type))

  } else if (type == "deseq"){ #Trying to implement idea from Anders and Huber Genome Biology 2010 paper.
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))   # Geometric means for each features.
    sizes <- apply(counts, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
    rawsizestr <- sizes  ## This part will be used to normalize test samples using train sample information.
    sizes <- sizes/sum(sizes)  # size factors for each sample
    fit <- outer(sizes, colsum, "*")   # s.i * g.j
    return(list(n = fit, sizes = sizes, geomeans = geomeans, rawsizestr = rawsizestr,
                lib.size = rowsum, gene.length = colsum, type = type))

  } else if (type == "tmm"){
    counts <- t(x)
    null.out <- calcTMMnormFactors(object = counts)

    lib.size <- colSums(counts)  ## Library sizes
    gene.length <- rowSums(counts)   # Gene lengths

    nf <- null.out$norm.factors   ## Normalization factors
    sizes <- nf / sum(nf)
    fit <- outer(sizes, gene.length, "*")

    return(list(n = fit, refSample = null.out$refSample, geomean = null.out$geomean,
                sizes = sizes, normalization.factors = nf, lib.size = lib.size, type = type))

  } else if (type == "none"){
    fit <- x
    lsize <- rowSums(x)

    return(list(n = fit, lib.size = lsize, type = type))
  }
}


#' @param null.out an object returned from \code{\link{NullModel}}.
#' @param xte a n-by-p count matrix or data frame of test set. These counts are normalized using training set parameters.
#'
#' @rdname NullModel
#'
#' @importFrom stats median quantile
#' @export
NullModelTest <- function(null.out, xte = NULL){
  # xte MUST be a n times p matrix - i.e. observations on the rows and features on the columns

  type <- null.out$type

  if (is.null(xte)){
    stop("'xte' cannot be NULL. Which data set should be normalized?")
  }

  if (type == "mle"){
    sizeste <- rowSums(xte)/sum(null.out$lib.size)
    nste <- outer(sizeste, null.out$gene.length, "*")
    rawsize.xte <- NULL

  } else if (type == "quantile"){   ### BU KISIM DÜZENLENECEK.
    sizeste <- pmax(1, apply(xte, 1, quantile, .75)) # don't want to scale by anything less than 1...
    rawsize.xte <- sizeste
    sizeste <- sizeste/sum(null.out$sample.qts.tr)
    nste <- outer(sizeste, null.out$gene.length, "*")

  } else if (type == "deseq"){
    countste <- t(xte)
    geomeans <- null.out$geomeans
    geom.xte  = countste / geomeans
    rawsize.xte = apply(geom.xte, 2, function(x) median(x[geomeans > 0], na.rm = TRUE))

    sizeste <- rawsize.xte / sum(null.out$rawsizestr)  #### PoiClaClu paketi rawsize toplamları olarak train setin değerini kullanıyor. Bu değerlerin toplamı
    #### 1'e eşit çıkmıyor. Eşitliğin 1 olabilmesi için test setinin rawsize değerlerine bölmek gerekiyor.
    #### Burada bir problem olabilir mi???
    # sizete <- rawsize.xte / sum(rawsize.xte)

    nste <- outer(sizeste, null.out$gene.length, "*")

  } else if (type == "tmm"){
    referenceSample <- null.out$refSample

    ## Check if the feature names of reference sample are in the same order with rawCounts.
    if (identical(colnames(xte), names(referenceSample))){
      xte <- rbind(xte, referenceSample)
    } else {
      stop(warning("Reference sample either does not have same features or the features are not in the same order as test set. Calculation stops.."))
    }

    tmp.out <- calcTMMnormFactors(t(xte), refColumn = nrow(xte), geomean = null.out$geomean, test.sample = TRUE)
    nf <- tmp.out$norm.factors   ## Normalization factors
    # lsize <- countsDGE.nf$samples$lib.size  ## Library sizes

    xte <- xte[-nrow(xte), ]   ## Remove reference sample from test data.

    # sizeste <- nf / sum(nf)
    sizeste <- nf / sum(null.out$normalization.factors)   ## sum of normalization factors from training set. "deseq" ile benzer mantıkla.
    rawsize.xte <- nf

    nste <- outer(sizeste, colSums(xte), "*")
    lib.size <- rowSums(xte)  ## Library sizes
    gene.length <- colSums(xte)

  } else if (type == "none"){


  }

  return(list(nste = nste, sizeste = sizeste, type = type, rawsizete = rawsize.xte))
}



