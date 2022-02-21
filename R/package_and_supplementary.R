###### Documentation for package and supplementary files (if any) ######

#' @title Classifying count data using Poisson/Negative Binomial linear discriminant analysis
#'
#' @description
#' This package applies linear discriminant analysis using Poisson (PLDA) and Negative Binomial (NBLDA) distributions for the classification of count data, such as gene expression data from RNA-sequencing. PLDA algorithms have been proposed by Witten (2011) through an R package \code{PoiClaClu}, which is available at CRAN. Dong et al. (2016) proposed an extension of PLDA to negative Binomial distribution. However, the algorithm is not provided through an R package. Hence, we develop an R package \code{NBLDA} to make the proposed algorithm available through CRAN. Detailed information about mathematical backgrounds is available in the references given below.
#'
#' @seealso \url{https://CRAN.R-project.org/package=PoiClaClu}
#'
#' \tabular{ll}{
#'   Package: \tab NBLDA\cr
#'   Type: \tab Package\cr
#'   License: \tab GPL (>= 2)\cr
#' }
#'
#' @author Dincer Goksuluk, Gokmen Zararsiz, Selcuk Korkmaz, A. Ergun Karaagaoglu
#'
#' -----------------
#'
#' \strong{Maintainers:}
#'
#' Dincer Goksuluk (Correspondence), \email{dincer.goksuluk@hacettepe.edu.tr}
#'
#' Gokmen Zararsiz, \email{gokmenzararsiz@erciyes.edu.tr}
#'
#' Selcuk Korkmaz, \email{selcukorkmaz@hotmail.com}
#'
#' @docType package
#' @name NBLDA-package
#' @rdname NBLDA-package
#'
#' @references Witten, DM (2011). Classification and clustering of sequencing data using a Poisson model.
#' Ann. Appl. Stat. 5(4), 2493--2518. doi:10.1214/11-AOAS493.
#'
#' Dong, K., Zhao, H., Tong, T., & Wan, X. (2016). NBLDA: negative binomial linear discriminant analysis for RNA-Seq data.
#' BMC Bioinformatics, 17(1), 369. http://doi.org/10.1186/s12859-016-1208-1
#'
#' @keywords package
NULL


#' @title Cervical cancer data
#'
#' @description Cervical cancer data measures the gene expression levels of 714 miRNAs of human samples. There are 29 tumor and 29 non-tumor cervical samples, and these two groups correspond to two separate classes.
#'
#' @format A data frame with 58 observations and 715 variables (including the class labels).
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2880020/}
#'
#' @references Witten, D., et al. (2010) Ultra-high throughput sequencing-based small RNA discovery and discrete statistical biomarker analysis in a collection of cervical tumours and matched controls. BMC Biology, 8:58
#'
#' @docType data
#' @name cervical
#' @rdname cervical
#'
#' @keywords cervical data
#'
#' @examples
#' \dontrun{
#' data(cervical)
#' }
NULL



