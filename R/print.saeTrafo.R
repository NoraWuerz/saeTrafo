# Prints a saeTrafo object

#' @export

print.NER <- function(x, ...) {

  throw_class_error(x, "NER")

  cat("Nested error regression model (under transformations)\n")

  cat("\n")
  cat("Out-of-sample domains: ", x$framework$N_dom_unobs, "\n")
  cat("In-sample domains: ", x$framework$N_dom_smp, "\n")

  if (x$transformation == "log") {
    transform_method <- data.frame(
      Transformation  = x$transformation,
      Shift_parameter = round(x$transform_param$shift_par, 3),
      row.names       = ""
    )
  } else if (x$transformation == "log.shift") {
    transform_method <- data.frame(
      Transformation  = x$transformation,
      Method          = x$method,
      Optimal_lambda  = x$transform_param$optimal_lambda,
      row.names       = ""
    )
  } else if (x$transformation == "no") {
    transform_method <- NULL
  }

  cat("\n")
  if (is.null(transform_method)) {
    cat("Transformation: No transformation \n")
  } else {
    cat("Transformation:\n")
    print(transform_method)
  }
  cat("\n")
  cat("Model fit:\n")
  cat("For model fit lme methods are applicable to saeTrafoObject$model \n")
  cat("where transformed_data equals smp_data transformed by function \n")
  cat("data_transformation using above given transformation and lambda \n")
  cat("and where fixed/list(fixed) equals ")
  print(x$fixed)
  cat("\n")
}
