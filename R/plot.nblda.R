
#' @title Plot Method for the \code{nblda} and \code{nblda_trained} Classes
#'
#' @description  This function is used to generate model performance plots using \code{\link[ggplot2:ggplot]{ggplot2}} functions.
#'
#' @param x a \code{nblda} object returned from \code{\link{trainNBLDA}} or \code{nblda_trained} object returned from \code{\link{nbldaTrained}}.
#' @param y same as \code{x} and not required to be defined. If \code{x} is missing or NULL, \code{nblda} or \code{nblda_trained} object is imported from \code{y}.
#' @param theme pre-defined plot themes. It can be defined outside \code{plot} function using ggplot's library. See examples.
#' @param ... further arguments to be passed to plotting function \code{\link[ggplot2]{ggplot}}.
#' @param metric which metric should be used in y-axis?
#' @param return should complete plot or a ggplot object from \code{ggplot} be returned? One may select "aes" in order to add plot layers
#' to returned ggplot aesthetics. See examples.
#'
#' @return A list of class \code{ggplot}.
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#'
#' @author Dincer Goksuluk
#'
#' @examples
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5,
#'                             DE = 0.8, allZero.rm = FALSE, tag.samples = TRUE)
#' x <- t(counts$x + 1)
#' y <- counts$y
#' xte <- t(counts$xte + 1)
#' ctrl <- nbldaControl(folds = 2, repeats = 2)
#'
#' fit <- trainNBLDA(x = x, y = y, type = "mle", tuneLength = 10,
#'                   metric = "accuracy", train.control = ctrl)
#'
#' plot(fit)
#'
#' # Use pre-defined theme
#' plot(fit, theme = "nblda")
#'
#' # Externally defining plot theme
#' plot(fit, theme = "default") + theme_dark(base_size = 14)
#'
#' # Return empty ggplot object and add layers.
#' plot(fit, theme = "nblda", return = "aes") +
#'   geom_point() + geom_line(linetype = 2)
#'
#' @name plot
#' @rdname plot
#'
#' @importFrom graphics plot
#' @import ggplot2
#'
#' @method plot nblda
plot.nblda <- function(x, y, ..., theme = c("nblda", "default"), metric = c("accuracy", "error", "sparsity"),
                       return = c("plot", "aes")){
  theme <- match.arg(theme)
  metric <- match.arg(metric)
  return <- match.arg(return)

  if (all(missing(x), missing(y))){
    stop("At least one of 'x' or 'y' should be given. Both can not be missing or NULL.")
  }

  if (missing(x) || is.null(x)){
    if (is.null(y)){
      stop("'y' can not be NULL when 'x' is missing or NULL.")
    }
    warning("'x' is not given. Plot object is imported from 'y'.")
    x <- y
  }

  if (!(class(x) %in% c("nblda", "nblda_trained"))){
    stop("'object' should be an object of classes 'nblda' or 'nblda_trained'.")
  }

  .theme_default <- theme_grey(base_size = 12)  # Default theme

  ## Set plot theme.
  if (theme == "nblda"){
    .theme <- theme_bw(base_size = 12) + theme(axis.title = element_text(face = "bold", colour = "#444444"),
                                               axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
                                               axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
                                               axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
                                               axis.text.y = element_text(margin = margin(0, 5, 0, 0)))
    line_color <- point_color <- "#3a8fff"
  } else {
    .theme <- .theme_default
    line_color <- point_color <- "black"
  }

  .data <- data.frame(nbldaTrained(x)@crossValidated$tuning.results)
  if (metric == "error"){
    aes_y <- "errors"
    y_lab <- "Average number of errors (cross-validated)"
  } else if (metric == "accuracy"){
    aes_y <- "accuracy"
    y_lab <- "Average classification accuracy (cross-validated)"
  } else {
    aes_y <- "nonzero"
    y_lab <- "Average number of selected features (cross-validated)"
  }

  p <- ggplot(.data, aes_string(x = "rho", y = aes_y), ...) +
    xlab("Threshold parameter (rho)") +
    ylab(y_lab) + .theme

  if (return == "plot"){
    if (nrow(.data) >= 2){
      p <- p + geom_line(color = line_color) +
        geom_point(color = point_color, pch = 21, fill = "white", size = 1.8)
    } else {
      p <- p + geom_point(color = point_color, pch = 21, fill = "white")
    }
  }
  p
}



#' @rdname plot
#' @method plot nblda_trained
plot.nblda_trained <- plot.nblda
