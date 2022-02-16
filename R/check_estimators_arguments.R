estimators_check <- function(object, indicator, MSE, CV) {

  if (is.null(object$MSE) && (MSE == TRUE || CV == TRUE)) {
    stop(paste("No MSE estimates in saeTrafo object: arguments MSE and CV have",
               "to be FALSE or a new emdi object with variance/MSE needs to be",
               "generated."))
  }
  if (!(inherits(MSE, "logical") && length(MSE) == 1)) {
    stop("MSE must be a logical value. Set MSE to TRUE or FALSE.")
  }
  if (!(inherits(CV, "logical") && length(CV) == 1)) {
    stop("CV must be a logical value. Set CV to TRUE or FALSE.")
  }
}
