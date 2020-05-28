# importFrom(sSeq,getAdjustDisp)
# importFrom(sSeq,getT)
# importFrom(sSeq,rowVars)

#' @importFrom stats quantile var cov
getT <- function(countsTable, sizeFactors = NULL, q.vec = NULL,
                 numPart = 1, propForSigma = c(0, 1), verbose = TRUE, shrinkTarget = NULL,
                 shrinkQuantile = NULL, shrinkVar = FALSE, eSlope = 0.05,
                 disp = NULL, dispXX = NULL, normalize = FALSE, lwd1 = 4.5,
                 cexlab1 = 1.2) {
  # This function is copied from the Bioconductor package "sSeq" in order to prevent
  # R-build error produced by "r-release-osx-x86_64" and "r-oldrel-osx-x86_64" machines.
  # The produced error is:
  #   Result: ERROR
  #      Package required but not available: ‘sSeq’

  if (!is.null(countsTable)) {
    counts = as.matrix(countsTable)
  }

  if (is.null(countsTable) & is.null(disp)) {
    stop("at least provide the initial dispersion estimates.")
  }

  if (is.null(sizeFactors) & !is.null(countsTable)) {
    sizeFactors = getNormFactor(countsTable)
  }

  if (is.null(eSlope)) {
    eSlope = 0.002
  } else {
    if (length(eSlope) > 1 & verbose)
      print("Note: only the first value in eSlope is used for tests.")
  }

  allAdjDisp = list()
  if (is.null(disp) & !is.null(countsTable)) {
    normc = as.matrix(t(t(counts) / sizeFactors))
    normc.m = rowMeans(normc)
    normc.v = rowVars(as.matrix(t(t(counts) / sqrt(sizeFactors))))
    s_st = mean(1 / sizeFactors)
    disp = (normc.v - normc.m) / (normc.m)^2
    normc.m[normc.m <= 0] = 1
    log.normc.m = log(normc.m)

  } else if (!is.null(dispXX)) {
    normc.m = dispXX
    normc.m[normc.m <= 0] = 1
    log.normc.m = log(normc.m)

  } else {
    normc.m = NULL
    log.normc.m = NULL
  }

  if (shrinkVar & is.null(disp)) {
    disp = normc.v
    if (verbose)
      print("Shrinkage estimates on variance are used.")
  } else {
    if (verbose)
      print("Shrinkage estimates on dispersion are used for the tests.")
  }

  disp[is.na(disp)] = 0
  disp[disp <= 0] = 0

  if (numPart == 1) {
    disp.m = mean(disp)
    asd.mle = round(mean((disp - disp.m)^2, na.rm = T), 4)
    rg.xx = quantile(disp[is.finite(disp)], prob = c(0.05, 0.995))
    xx = seq(rg.xx[1], rg.xx[2], length.out = 200)
    asd = rep(0, length(xx))

    for (i in 1:length(xx)) {
      allAdjDisp[[i]] = equalSpace(disp, log.normc.m, 1,
                                   propForSigma = propForSigma, shrinkTarget = xx[i],
                                   vb = FALSE)
      allAdjDisp[[i]] = pmax(1e-08, allAdjDisp[[i]])
      names(allAdjDisp[[i]]) = 1:length(disp)
      asd[i] = mean((allAdjDisp[[i]] - disp)^2, na.rm = T)
    }

    diff.q = diff.asd = rep(0, length(asd))
    maxASD = max(asd, na.rm = T)
    maxASD.pnt = which(asd == maxASD)

    for (i in 1:length(asd)) {
      diff.asd[i] = maxASD - asd[i]
      diff.q[i] = xx[maxASD.pnt] - xx[i]
    }

    numAdjPoints = 6
    len.asd = length(asd) - numAdjPoints + 1
    slope1 = rep(1, len.asd)

    if (normalize) {
      xx1 = xx / sd(xx)
      yy1 = asd / sd(asd)
      eSlope = eSlope * 5
    } else {
      xx1 = xx
      yy1 = asd
    }

    for (i in 1:len.asd) {
      slope1.xx = xx1[i:(i + numAdjPoints - 1)]
      slope1.yy = yy1[i:(i + numAdjPoints - 1)]
      slope1[i] = cov(slope1.xx, slope1.yy) / var(slope1.xx)
    }

    maxSlope1 = max(abs(slope1))
    maxSlope1.pnt = which(abs(slope1) == maxSlope1)
    sub.slope1 = abs(slope1)[maxSlope1.pnt:len.asd]
    sub.diff.asd = diff.asd[maxSlope1.pnt:length(diff.asd)]
    pred.diff = matrix(NA, nrow = length(sub.diff.asd), ncol = numAdjPoints)

    for (i in 1:length(sub.diff.asd)) {
      for (j in 1:numAdjPoints) {
        if (i - j >= 0) {
          pred.diff[i, j] = sub.diff.asd[i] / sub.slope1[i - j + 1]
        }
      }
    }

    max.pred = max(pred.diff, na.rm = T)
    max.rowInd = which(apply(pred.diff, 1, max, na.rm = T) == max.pred)
    temp.max.pnt = max.rowInd + maxSlope1.pnt - 1 - ceiling(numAdjPoints/2)
    max.pnt = rep(0, length(eSlope))

    for (k in 1:length(eSlope)) {
      max.pnt[k] = temp.max.pnt[1]
      while (-slope1[max.pnt[k] - ceiling(numAdjPoints/2)] <
             eSlope[k] & slope1[max.pnt[k] - ceiling(numAdjPoints/2)] < 0) {
        max.pnt[k] = max.pnt[k] - 1
      }
    }

    if (!is.null(shrinkQuantile)) {
      max.pnt = 1
      q.vec = shrinkQuantile
      adjDisp1 = equalSpace(disp, log.normc.m, numPart,
                            propForSigma = propForSigma, shrinkQuantile = q.vec[1],
                            vb = FALSE)
      adjDisp1 = pmax(1e-08, adjDisp1)
      names(adjDisp1) = 1:length(disp)
      asd.target = mean((adjDisp1 - disp)^2, na.rm = T)
      target = round(quantile(disp, prob = q.vec[1]), 3)
    }

    if (!is.null(shrinkTarget)) {
      max.pnt = 1
      disp.tm = c(disp[!is.na(disp)], shrinkTarget[1])
      q.vec = round(rank(disp.tm)[disp.tm == shrinkTarget[1]] / length(disp[!is.na(disp)]), 3)
      adjDisp1 = equalSpace(disp, log.normc.m, numPart,
                            propForSigma = propForSigma, shrinkQuantile = q.vec[1],
                            vb = FALSE)
      adjDisp1 = pmax(1e-08, adjDisp1)
      names(adjDisp1) = 1:length(disp)
      asd.target = mean((adjDisp1 - disp)^2, na.rm = T)
      target = shrinkTarget[1]
    }

    if (is.null(shrinkQuantile) & is.null(shrinkTarget)) {
      target = asd.target = q.vec = rep(0, length(eSlope))
      for (k in 1:length(eSlope)) {
        target[k] = xx[max.pnt[k]][1]
        asd.target[k] = asd[max.pnt[k]]
        disp.tm = c(disp[!is.na(disp)], target[k])
        q.vec[k] = round(rank(disp.tm)[disp.tm == target[k]] / length(disp[!is.na(disp)]), 3)
      }
    }

    if (verbose) {
      if (!is.null(shrinkQuantile) | !is.null(shrinkTarget)) {
        print(paste("The selected shrink target is",
                    target[1]))
        print(paste("The selected shrink quantile is",
                    q.vec[1]))
      }
      else {
        print(paste("The shrink target is", target[1]))
        print(paste("The shrink quantile is", q.vec[1]))
      }
    }

    return(list(q = q.vec[1], target = target[1]))

  } else if (numPart > 1) {
    if (is.null(log.normc.m)) {
      stop("Error in getT: log.normc.m can not be NULL.")
    }
    out = getTgroup(y = disp, x = log.normc.m, numPart = numPart,
                    verbose = verbose, eSlope = eSlope,
                    lwd1 = lwd1, cexlab1 = cexlab1)
    return(out)

  } else {
    stop("Error: numPart must be a non-negative integer.")
  }
}


#' @importFrom stats quantile var
getAdjustDisp <- function(obs, propForSigma = c(0.5, 1), shrinkTarget = NULL,
                          shrinkQuantile = NULL, verbose = TRUE){
  # This function is copied from the Bioconductor package "sSeq" in order to prevent
  # R-build error produced by "r-release-osx-x86_64" and "r-oldrel-osx-x86_64" machines.
  # The produced error is:
  #   Result: ERROR
  #      Package required but not available: ‘sSeq’

  obs[is.na(obs)] = 0
  if (is.null(shrinkTarget)) {
    upBound = quantile(obs, prob = shrinkQuantile, na.rm = T)
    if (verbose) {
      print(paste("shrink toward ", shrinkTarget, " (",
                  shrinkQuantile, "th quantile).", sep = ""))
    }
  } else {
    upBound = shrinkTarget
    if (verbose) {
      print(paste("shrink toward ", shrinkTarget, ".", sep = ""))
    }
  }
  if (is.null(propForSigma)) {
    subobs = obs[obs >= upBound & obs <= quantile(obs, prob = 0.999)]
    S.mat = var(subobs, na.rm = T)

  } else if (length(propForSigma) == 2) {
    subobs = obs
    rg = quantile(subobs[is.finite(subobs)], na.rm = T, prob = propForSigma)
    subobs = subobs[subobs >= rg[1] & subobs <= rg[2]]
    S.mat = var(subobs[is.finite(subobs)], na.rm = T)

  } else if (length(propForSigma) == 1 & is.numeric(propForSigma)) {
    S.mat = propForSigma

  } else if (is.na(propForSigma)) {
    subobs = obs[is.finite(obs)]
    S.mat = var(subobs[is.finite(subobs)], na.rm = T)

  } else {
    stop(paste("if don't know the empirical value on the variance of",
               "dispersion, please set it as NULL."))
  }
  cmp = data.frame(mean = mean(obs, na.rm = T), sigmasq.part = S.mat)
  mean.mat = rep(upBound, length(obs))
  dif.mat = obs - mean.mat
  dif2.mat = sum(dif.mat^2)
  deta = 1 - ((length(obs) - 2) * S.mat / (dif2.mat))
  jsDiff = pmax(0, deta) * dif.mat
  jsest = jsDiff + mean.mat
  return(list(adj = jsest, cmp = cmp))
}


#' @importFrom stats median
getNormFactor <- function (countsTable1){
  countsTable1.log = log(countsTable1)
  row.mean1.log = rowMeans(countsTable1.log)
  geo.dev1.log = countsTable1.log - row.mean1.log
  apply(geo.dev1.log, 2, function(x) {
    exp(median(x[is.finite(x)]))
  })
}

#' @importFrom stats var
rowVars <- function (x){
  apply(x, 1, var, na.rm = T)
}

equalSpace <- function (y, x = NULL, numcls = 1, propForSigma = c(0, 1), shrinkTarget = NULL,
                        shrinkQuantile = 0.975, vb = TRUE){
  if (numcls == 1 | is.null(x))
    return(getAdjustDisp(y, propForSigma = propForSigma,
                         shrinkTarget, shrinkQuantile, verbose = vb)$adj)

  if (!is.null(shrinkTarget) & length(shrinkTarget) != numcls) {
    print(paste("Warning: the number of shrink targes is unequal to the",
                "number of pre-decied groups. Only the first target is used."))
    shrinkTarget = shrinkTarget[1]
    numcls = 1
  }

  if (sum(is.na(x)) > 0)
    print("The NA values in the dependent variable were ignored.")

  if (length(y) != length(x))
    stop(paste("Error: check the input of equalSpace. y and x have",
               "unequal lengths in equalSpace function."))

  rgx = range(x[x > -Inf])
  cut = seq(from = rgx[1], to = rgx[2], length = numcls + 1)
  cls = rep(1, length(y))
  cls[x <= cut[2]] = 1
  cls[x > cut[numcls]] = numcls

  for (i in 2:(numcls - 1)) {
    cls[x > cut[i] & x <= cut[i + 1]] = i
  }

  sizes = tapply(rep(1, length(cls)), cls, sum)
  js = y
  mean.y = mean(y)

  for (i in 1:length(sizes)) {
    if (sizes[i] > 2) {
      x.sub = x[cls == i]
      if (!is.null(shrinkTarget)) {
        mixr = getAdjustDisp(y[cls == i], propForSigma = propForSigma,
                             shrinkTarget[i], shrinkQuantile, verbose = vb)
      } else {
        mixr = getAdjustDisp(y[cls == i], propForSigma = propForSigma,
                             shrinkTarget = NULL, shrinkQuantile = shrinkQuantile,
                             verbose = vb)
      }
      js[cls == i] = mixr$adj
    } else {
      js[cls == i] = mean.y
    }
  }
  return(js)
}

#' @importFrom stats quantile cov var
getTgroup <- function (y, x, numPart = 10, plotASD = FALSE, verbose = FALSE,
                       eSlope = 0.05, lwd1 = 4.5, cexlab1 = 1.2){

  rgx = range(x[is.finite(x)])
  cut = seq(from = rgx[1], to = rgx[2], length = numPart + 1)
  cls = rep(1, length(y))
  cls[x <= cut[2]] = 1
  cls[x > cut[numPart]] = numPart

  for (i in 2:(numPart - 1)) {
    cls[x > cut[i] & x <= cut[i + 1]] = i
  }

  sizes = tapply(rep(1, length(cls)), cls, sum)
  qall.vec = targetall = rep(1, numPart)

  for (gp in 1:numPart) {
    allAdjy = list()
    y1 = y[cls == gp]
    x1 = x[cls == gp]
    y.m = mean(y1)
    asd.mle = round(mean((y1 - y.m)^2, na.rm = T), 4)
    rg.xx = quantile(y[is.finite(y1)], prob = c(0.05, 0.995))
    xx = seq(rg.xx[1], rg.xx[2], length.out = 200)
    asd = rep(0, length(xx))

    for (i in 1:length(xx)) {
      allAdjy[[i]] = equalSpace(y = y1, x = x1, numcls = 1,
                                shrinkTarget = xx[i], vb = FALSE)
      allAdjy[[i]] = pmax(1e-08, allAdjy[[i]])
      names(allAdjy[[i]]) = 1:length(y1)
      asd[i] = mean((allAdjy[[i]] - y1)^2, na.rm = T)
    }

    diff.q = diff.asd = rep(0, length(asd))
    maxASD = max(asd, na.rm = T)
    maxASD.pnt = which(asd == maxASD)
    maxASD.pnt = max(maxASD.pnt)

    for (i in 1:length(asd)) {
      diff.asd[i] = maxASD - asd[i]
      diff.q[i] = xx[maxASD.pnt] - xx[i]
    }

    numAdjPoints = 6
    len.asd = length(asd) - numAdjPoints + 1
    slope1 = rep(1, len.asd)
    xx1 = xx
    y11 = asd

    for (i in 1:len.asd) {
      slope1.xx = xx1[i:(i + numAdjPoints - 1)]
      slope1.y1 = y11[i:(i + numAdjPoints - 1)]
      slope1[i] = cov(slope1.xx, slope1.y1) / var(slope1.xx)
    }

    maxSlope1 = max(abs(slope1))
    maxSlope1.pnt = which(abs(slope1) == maxSlope1)
    sub.slope1 = abs(slope1)[maxSlope1.pnt:len.asd]
    sub.diff.asd = diff.asd[maxSlope1.pnt:length(diff.asd)]
    pred.diff = matrix(NA, nrow = length(sub.diff.asd), ncol = numAdjPoints)

    for (i in 1:length(sub.diff.asd)) {
      for (j in 1:numAdjPoints) {
        if (i - j >= 0) {
          pred.diff[i, j] = sub.diff.asd[i] / sub.slope1[i - j + 1]
        }
      }
    }

    max.pred = max(pred.diff, na.rm = T)
    max.rowInd = which(apply(pred.diff, 1, max, na.rm = T) == max.pred)
    temp.max.pnt = max.rowInd + maxSlope1.pnt - 1 - ceiling(numAdjPoints / 2)
    max.pnt = rep(0, length(eSlope))

    for (k in 1:length(eSlope)) {
      max.pnt[k] = temp.max.pnt[1]
      tm1 = -slope1[max.pnt[k] - ceiling(numAdjPoints/2)]
      tm2 = -tm1
      while (!is.na(tm1) & tm1[1] < eSlope[k] & tm2[1] < 0) {
        max.pnt[k] = max.pnt[k] - 1
        tm1 = -slope1[max.pnt[k] - ceiling(numAdjPoints/2)]
        tm2 = -tm1
        if (length(tm1) == 0)
          break
      }
    }

    target = asd.target = q.vec = rep(0, length(eSlope))
    for (k in 1:length(eSlope)) {
      target[k] = xx[max.pnt[k]][1]
      asd.target[k] = asd[max.pnt[k]][1]
      y.tm = c(y1[!is.na(y1)], target[k])
      q.vec[k] = round(rank(y.tm)[y.tm == target[k]]/length(y1[!is.na(y1)]), 3)
    }

    if (verbose) {
      print(paste("In group", gp, "the average of the values on",
                  "X-axis for", sizes[gp], "genes is", mean(x1, na.rm = T)))
      print(paste("shrinkTarget ", target[1], " and shrinkQuantile ",
                  q.vec[1], ".", sep = ""))
    }
    qall.vec[gp] = q.vec[1]
    targetall[gp] = target[1]
  }

  return(list(q = qall.vec, target = targetall))
}
