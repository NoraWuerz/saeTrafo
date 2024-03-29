# Extract model coefficients of saeTrafo objects -------------------------------

#' @aliases coefficients
#' @export
#' @method coef NER
#' @importFrom stats coef coefficients

coef.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  coef(object$model)
}

# Confidence intervals of an saeTrafo object -----------------------------------

#' @export
#' @method confint NER
#' @importFrom nlme intervals

confint.NER <- function(object, parm = NULL, level = 0.95, ...) {
  throw_class_error(object, "NER")
  if (!is.null(parm)) {
    confidence_intervals <- intervals(object$model, level = level)$fixed
    subset(confidence_intervals, rownames(confidence_intervals) %in% parm)
  } else {
    intervals(object$model, level = level)$fixed
  }
}

# Extracts family object of saeTrafo object ------------------------------------

#' @export
#' @method family NER
#' @importFrom stats family gaussian

family.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  gaussian(link = "identity")
}

# Extract fitted values of saeTrafo objects ------------------------------------

#' @aliases fitted.values
#' @export
#' @method fitted NER
#' @importFrom stats fitted

fitted.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  fitted(object$model, ...)
}

# Extract the model formula of an saeTrafo object ------------------------------

#' @export
#' @method formula NER

formula.NER <- function(x, ...) {
  throw_class_error(x, "NER")
  x$fixed
}

# Extract log-Likelihood of saeTrafo objects -----------------------------------
#' @export
#' @method logLik NER

logLik.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  message("Estimation approach used is reml: ", round(object$model$logLik, 5))
  invisible(object$model$logLik)
}

# Extract the number of `observations´ from a fit of an saeTrafo object -------
#' @export
#' @method nobs NER
#' @importFrom stats nobs

nobs.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  object$framework$N_smp
}

#-------------------------------------------------------------------------------
#' Predictions from saeTrafo objects
#'
#' Method \code{predict.NER} extracts the direct estimates, the empirical
#' best linear unbiased or empirical best predictors for all domains from an
#' saeTrafo object.
#'
#' @param object an object of type "saeTrafo".
#' @param ... additional arguments that are not used in this method.
#' @return Data frame with domain predictors.
#'
#' @examples
#' # Examples for Predictions from saeTrafo objects
#'
#' # Load Data
#' data("eusilcA_smp")
#' data("pop_area_size")
#' data("pop_mean")
#' data("pop_cov")
#'
#' # Nested error regression model
#' NER_model <- NER_Trafo(fixed = eqIncome ~ gender + eqsize + cash +
#'                        self_empl + unempl_ben + age_ben + surv_ben +
#'                        sick_ben + dis_ben + rent + fam_allow + house_allow +
#'                        cap_inv + tax_adj,
#'                        smp_domains = "district",
#'                        pop_area_size = pop_area_size,
#'                        pop_mean = pop_mean, pop_cov = pop_cov,
#'                        smp_data = eusilcA_smp)
#'
#' predict(NER_model)
#'
#' @seealso \code{\link{saeTrafoObject}}, \code{\link{NER_Trafo}}
#' @export
#' @method predict NER

predict.NER <- function(object, ...) {
  object$ind
}

# Extract residuals of saeTrafo objects ----------------------------------------

#' @aliases resid
#' @export
#' @method residuals NER
#' @importFrom stats residuals resid

residuals.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  residuals(object$model, ...)
}

# Extract residual standard deviation of saeTrafo objects ----------------------

#' @export
#' @method  sigma NER
#' @importFrom stats sigma

sigma.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  object$model$sigma
}

# Constructs a terms object from an saeTrafo object ----------------------------

#' @export
#' @method terms NER
#' @importFrom stats aov terms

terms.NER <- function(x, ...) {
  throw_class_error(x, "NER")
  terms(aov(x$fixed, x$framework$smp_data))
}

# Extract variance-covariance matrix of the main parameters --------------------
# of saeTrafo objects

#' @export
#' @method vcov NER
#' @importFrom stats vcov

vcov.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  vcov(object$model, ...)
}
