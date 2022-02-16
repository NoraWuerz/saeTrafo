plot_check <- function(x, label, color, cooks, range) {

  if (!any(inherits(x, which = TRUE, c("NER")))) {
    stop("First object needs to be of class saeTrafo, model.")
  }
  if (is.null(label) || (!(label == "orig" || label == "no_title" ||
                           label == "blank") && !inherits(label, "list"))) {
    stop(paste("label can be either one of the following characters 'orig',",
               "'no_title' or 'blank' or a list as specified in",
               "help(plot.saeTrafo)."))
  }
  if (length(color) != 2 || !is.vector(color)) {
    stop(paste("color needs to be a vector of length 2 defining the two colors",
               "for the diagnostic plot. See also help(plot.saeTrafo)."))
  }
  if (!(inherits(cooks, "logical") && length(cooks) == 1)) {
    stop(paste("cooks needs to be a logical value. Set na.rm to TRUE or FALSE.",
               "See also help(plot.saeTrafo)."))
  }
  if (!is.null(range) && !inherits(range, "numeric")) {
    stop(paste("range must be a sequence determining the range of the x-axis",
               "for plots of the optimal parameter.. Set na.rm to TRUE or",
               "FALSE. See also help(plot.saeTrafo)."))
  }
}
