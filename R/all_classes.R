# All classes..
###### 1. SubClass: nblda_trained  #######
#' @title \code{nblda_trained} object
#'
#' @description This object is the subclass for the NBLDA package. It stores the cross-validated results and the final model.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{crossValidated}:}{a list. Returns the results from cross-validation.}
#'   \item{\code{finalModel}:}{a list with the elements from the final model that is fitted using optimum model parameters from the cross-validated model.}
#'   \item{\code{control}:}{a list with controlling parameters for fitting NBLDA classifier.}
#' }
#'
#' @author Dincer Goksuluk
#'
#' @docType class
#' @name nblda_trained-class
#' @rdname nblda_trained-class
#'
#' @exportClass nblda_trained
setClass("nblda_trained", slots = c(crossValidated = "list", finalModel = "list", control = "list"))


###### 2. SubClass: nblda_input  #######
setClassUnion("count.data", c("matrix", "data.frame"))
setClassUnion("class.labels", c("numeric", "integer", "factor"))

#' @title \code{nblda_input} object
#'
#' @description This object is the subclass for the NBLDA package. It stores input objects, i.e., count data and class labels.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{x}:}{a data.frame or matrix containing the count data input for the NBLDA classifier.}
#'   \item{\code{y}:}{a vector of length equal to the number of rows of x. This is the class label of each subject. Should be either a numeric vector or factor.}
#' }
#'
#' @author Dincer Goksuluk
#'
#' @docType class
#' @name nblda_input-class
#' @rdname nblda_input-class
#'
#' @exportClass nblda_input
setClass("nblda_input", slots = c(x = "count.data", y = "class.labels"))

###### 3. Class: nblda #####
#' @title \code{nblda} object
#'
#' @description This object is the main class for the NBLDA package. It stores inputs, results, and call info for the trained model.
#'
#' @details Objects can be created by calls of the form \code{new("nblda", ...)}. This type of object is returned from \code{trainNBLDA} function of the \code{NBLDA} package. It is then used in \code{predict} function for predicting class labels of new samples.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{input}:}{an \code{nblda_input} object including the count matrix (or data.frame) and class labels.}
#'   \item{\code{result}:}{an \code{nblda_trained} object with elements from the cross-validated and final models.}
#'   \item{\code{call}:}{a call expression.}
#' }
#'
#' @author Dincer Goksuluk
#'
#' @docType class
#' @name nblda-class
#' @rdname nblda-class
#'
#' @seealso \code{\linkS4class{nblda_trained}}, \code{\linkS4class{nblda_input}}
#'
#' @exportClass nblda
setClass("nblda",
         slots = c(input = "nblda_input",   # raw and transformed data
                   result = "nblda_trained",  # crossValidatedResults, finalModel, control
                   call = "call"))

# setValidity("discrete.train", function(object){
#   TRUE    ##### Slots MUST BE validated
# })


# setOldClass(c("confusionMatrix", "train"))
# setClassUnion("MLSeq.train", c("train", "voom.train", "discrete.train"))
# setClassUnion("confMat", c("confusionMatrix"))

