GoodnessOfFit <- function(x, type){
  ns <- NullModel(x, type = type)$n
  return(sum(((x - ns)^2)/ns, na.rm=TRUE))
}

permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}

# importFrom caret createFolds
# foldIndex <- function(data = NULL, n = NULL, nFolds = 5, repeats = 2){
#   if (!is.null(data)){
#     n = nrow(data)
#   }
#
#   indIn <- indOut <- list()
#
#   for (j in 1:repeats){
#     tmpIn = createFolds(y = 1:n, k = nFolds, list = TRUE, returnTrain = TRUE)
#     tmpOut = lapply(tmpIn, function(x)c(1:n)[-x])
#
#     indIn = c(indIn, tmpIn)
#     indOut = c(indOut, tmpOut)
#   }
#
#   nms = paste(rep(paste("Fold", 1:nFolds, sep = ""), repeats),
#               rep(paste(".Rep", 1:repeats, sep = ""), c(rep(nFolds, repeats))), sep = "")
#   names(indIn) <- names(indOut) <- nms
#   return(list(indexIn = indIn, indexOut = indOut))
# }

# Create folds and determine the indices of samples in each fold.
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}


#' @importFrom stats sd
se <- function(vec){
  ## Standart error of a given vector of parameter estimation.
  return(sd(vec)/sqrt(length(vec)))
}

Soft <- function(x, a){
  return(sign(x)*pmax(abs(x) - a, 0))
}

# Find d.kj estimates
GetD <- function(ns, x, y, rho = NULL, beta){
  if (is.null(rho)){
    rho <- 0
    warning("'rho' is not given. It is set to rho = 0.")
  }

  if (missing(beta) || is.null(beta)){
    beta <- 1
  }

  uniq <- sort(unique(y))
  ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
  for (k in 1:length(uniq)){
    a <- colSums(x[y == uniq[k], ]) + beta
    b <- colSums(ns[y == uniq[k], ]) + beta
    ds[k, ] <- 1 + Soft(a/b - 1, rho/sqrt(b))
  }
  return(ds)
}


###### TMM Normalization functions ######
.calcFactorQuantile_modified <- function (data, lib.size, p = 0.75){
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
  f
}

.calcFactorWeighted_modified <- function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
                                          sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10){
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs))
    nO <- sum(obs)
  else nO <- libsize.obs
  if (is.null(libsize.ref))
    nR <- sum(ref)
  else nR <- libsize.ref
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  if (max(abs(logR)) < 1e-06)
    return(1)
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >=
                                                       loS & rank(absE) <= hiS)
  if (doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],
                                                   na.rm = TRUE)
  else f <- mean(logR[keep], na.rm = TRUE)
  if (is.na(f))
    f <- 0
  2^f
}

calcTMMnormFactors <- function(object, lib.size = NULL, refColumn = NULL, logratioTrim = 0.3, geomean = NULL,
                               sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, test.sample = FALSE, ...){
  # Args:
  #   test.sample: logical. TRUE ise test örneklemi için hesaplama yapıyor.
  #   geomean: train setinde hesaplanan f değerlerinin geometrik ortalaması.
  #   refColumn: referans sample'in kolon numarası. Train setinden alınan sample. object içerisine eklenmiş olarak gönderilmesi gerekiyor.

  x <- as.matrix(object)
  if (any(is.na(x))) stop("NA counts not permitted")
  if (is.null(lib.size)) lib.size <- colSums(x)
  if (any(is.na(lib.size))) stop("NA lib.sizes not permitted")

  allzero <- .rowSums(x > 0, nrow(x), ncol(x)) == 0

  if (any(allzero)) x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1) method = "none"

  f75 <- .calcFactorQuantile_modified(data = x, lib.size = lib.size, p = 0.75)
  if (is.null(refColumn)) refColumn <- which.min(abs(f75 - mean(f75)))
  if (length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
  f <- rep(NA, ncol(x))
  for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted_modified(obs = x[ ,i], ref = x[ ,refColumn], libsize.obs = lib.size[i],
                                                            libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                            sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)

  if (is.null(geomean)){
    geomean <- exp(mean(log(f)))
  }
  f <- f/geomean

  if (test.sample){
    f <- f[-refColumn]
  }
  return(list(norm.factors = f, geomean = geomean, refSample = x[ ,refColumn], type = "tmm"))
}
