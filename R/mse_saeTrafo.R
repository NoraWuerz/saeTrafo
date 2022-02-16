# function for extracting MSE and CVs from saeTrafo object

mse_saeTrafo <- function(object, indicator = "all", CV = FALSE) {

  if (is.null(object$MSE) && CV == TRUE) {
    stop(paste("No MSE estimates in saeTrafo object: arguments MSE and CV have",
               "to be FALSE or a new saeTrafo object with variance/MSE needs",
               "to be generated."))
  }

  if (inherits(object, "fh")) {
    object$MSE <- object$MSE[, c("Domain", "Direct", "FH")]
    object$ind <- object$ind[, c("Domain", "Direct", "FH")]
  }

  all_cv <- sqrt(object$MSE[, -1]) / object$ind[, -1]

  if (any(indicator == "fh") || any(indicator == "FH")) {
    ind <- object$MSE[, c("Domain", "FH")]
    ind_cv <- cbind(Domain = object$MSE[, 1], all_cv)
    ind_name <- "Fay-Herriot estimates"
  } else if (any(class(object) == "NER")) {
    ind <- object$MSE[, c("Domain", "Mean")]
    ind_cv <- data.frame(Domain = object$MSE[, 1], Mean = all_cv)
    ind_name <- "Nested error regression model estimates"
  } else if (any(indicator == "Direct") || any(indicator == "direct")) {
    ind <- object$MSE[, c("Domain", "Direct")]
    ind_cv <- cbind(Domain = object$MSE[, 1], all_cv)
    ind_name <- "Direct estimates used in Fay-Herriot approach"
  } else {
    selection <- colnames(object$MSE[-1]) %in% indicator
    ind <- object$MSE[, c(TRUE, selection)]
    ind_cv <- data.frame(Domain = object$MSE[, 1], all_cv[, selection])
    colnames(ind_cv) <- colnames(ind)
    ind_name <- paste(unique(indicator), collapse = ", ")
  }

  if (CV == FALSE) {
    mse_saeTrafo <- list(ind = ind, ind_name = ind_name)
  } else {
    mse_saeTrafo <- list(ind = ind, ind_cv = ind_cv, ind_name = ind_name)
  }

  class(mse_saeTrafo) <- "mse.saeTrafo"

  return(mse_saeTrafo)
}
