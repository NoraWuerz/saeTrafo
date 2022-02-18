#' Extract fixed effects from an saeTrafo object
#'
#' Methods \code{fixef.NER} extract the fixed effects from an saeTrafo object of
#' class "NER".
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return A vector containing the fixed effects is returned.
#' @details The alias \code{fixed.effects} can also be used instead of
#' \code{fixef}. The generic function \code{fixef} is imported from package
#' \pkg{nlme} and re-exported to make the S3-methods available, even though the
#' \pkg{nlme} package itself is not loaded or attached. For default
#' documentation,
#' see \code{\link[nlme]{fixed.effects}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{fixed.effects}}
#' @examples
#'
#' # Example to extract fixed effects
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' fixef(NER_model)
#'
#' @aliases fixed.effects
#' @name fixef
#' @importFrom nlme fixef fixed.effects
#' @export fixed.effects
#' @export fixef
fixef

#' @export fixef.NER
#' @export
#' @rdname fixef
fixef.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  object$model$coefficients$fixed
}

#' @export fixed.effects.NER
#' @export
#' @rdname fixef
fixed.effects.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  object$model$coefficients$fixed
}


#' Extract saeTrafo object data
#'
#' Methods \code{getData.NER} extract the data frame used to fit the model.
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return Data frame used to fit the model. For "NER" the (untransformed)
#' sample data is returned.
#' @details The generic function \code{getData} is imported from package
#' \pkg{nlme} and re-exported to make the S3-methods available, even though the
#' \pkg{nlme} package itself is not loaded or attached. For default
#' documentation, see \code{\link[nlme]{getData}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{getData}}
#' @examples
#'
#' # Example to extract object data
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' getData(NER_model)
#' @name getData
#' @importFrom nlme getData
#' @export getData
getData

#' @export getData.NER
#' @export
#' @rdname getData
#'
getData.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  if (object$transformation != "no") {
    message("The untransformed sample data set of the NER object is returned.")
  }
  object$framework$smp_data
}


#' Extract grouping factors from an saeTrafo object
#'
#' Method \code{getGroups.NER} extract grouping factors from a saeTrafo object.
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return A vector containing the grouping factors.
#' @details The generic function \code{getGroups} is imported from package nlme and
#' re-exported to make the S3-methods available, even though the nlme package
#' itself is not loaded or attached. For default documentation,
#' see \code{\link[nlme]{getGroups}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{getGroups}}
#'
#' @examples
#' # Example to extract grouping factors
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' getGroups(NER_model)
#' @name getGroups
#' @importFrom nlme getGroups
#' @export getGroups
getGroups

#' @export getGroups.NER
#' @export
#' @rdname getGroups
getGroups.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  object$framework$smp_domains_vec
}


#' Extract grouping formula from a saeTrafo object
#'
#' Methods \code{getGroupsFormula.NER} extract the grouping formula from an
#' saeTrafo object.
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return A one-sided formula.
#' @details The generic function \code{getGroupsFormula} is imported from package
#' \pkg{nlme} and re-exported to make the S3-methods available, even though the
#' \pkg{nlme} package itself is not loaded or attached. For default
#' documentation, see \code{\link[nlme]{getGroupsFormula}}.
#' @seealso \code{\link{NER_Trafo}} \code{\link[nlme]{getGroupsFormula}}
#'
#' @examples
#' # Example to extract grouping formula
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' getGroupsFormula(NER_model)
#' @name getGroupsFormula
#' @importFrom nlme getGroupsFormula
#' @export getGroupsFormula
getGroupsFormula

#' @export getGroupsFormula.NER
#' @export
#' @rdname getGroupsFormula
getGroupsFormula.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  eval(parse(text = paste("~", object$framework$smp_domains)))
}


#' Extract response variable from an saeTrafo object
#'
#' Methods \code{getResponse.NER} extract the response variable from a saeTrafo
#' object.
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return Vector containing the response variable.
#' @details The generic function \code{getResponse} is imported from package
#' nlme and re-exported to make the S3-methods available, even though the nlme
#' package itself is not loaded or attached. For default documentation,
#' see \code{\link[nlme]{getResponse}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{getResponse}}
#'
#' @examples
#' # Example to extract the response variable
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' getResponse(NER_model)
#' @name getResponse
#' @importFrom nlme getResponse
#' @export getResponse
getResponse

#' @export getResponse.NER
#' @export
#' @rdname getResponse
getResponse.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  makeXY(object$fixed, object$framework$smp_data)$y
}


#' Extract variance-covariance matrix from an saeTrafo object
#'
#' Method \code{getVarCov.NER} extract the variance-covariance matrix from a
#' fitted model of class "NER".
#
#' @param obj an object of type "NER".
#' @param individuals vector of levels of the in-sample domains can be specified
#' for the types "\code{conditional}" or "\code{marginal}".
#' @param type a character that determines the type of variance-covariance matrix.
#' Types that can be chosen
#' (i) random-effects variance-covariance matrix ("\code{random.effects}"),
#' (ii) conditional variance-covariance matrix ("\code{conditional}"),
#' (iii) marginal variance-covariance matrix ("\code{marginal}"). Defaults to
#' "\code{random.effects}".
#' @param ... additional arguments that are not used in this method.
#' @return A variance-covariance matrix or a list of variance-covariance matrices,
#' if more than one individual is selected. For method \code{getVarCov.NER}, the
#' dimensions of the matrices are 1 x 1 for type "\code{random.effects}" and
#' number of in-sample domains x number of in-sample domains for types
#' "\code{conditional}" and "\code{marginal}".
#' @details The generic function \code{getVarCov} is imported from package
#' \pkg{nlme} and re-exported to make the S3-methods available, even though the
#' \pkg{nlme} package itself is not loaded or attached. For default
#' documentation, see \code{\link[nlme]{getVarCov}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{getVarCov}}
#'
#' @examples
#' # Example to extract variance-covariance matrix
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' getVarCov(NER_model)
#' @name getVarCov
#' @importFrom nlme getVarCov
#' @export getVarCov
getVarCov

#' @export getVarCov.NER
#' @export
#' @rdname getVarCov

getVarCov.NER <- function(obj, individuals = 1, type = "random.effects", ...) {
  throw_class_error(obj, "NER")

  if (is.null(type) || !(type == "random.effects"
                         || type == "conditional"
                         || type == "marginal")) {
    stop(paste("The three options for type are ''random.effects'',",
               "''conditional'' or ''marginal''."))
  }

  getVarCov(obj$model, individuals = individuals, type = type)

}


#' Confidence intervals on coefficients of an saeTrafo object
#'
#' Method \code{intervals.NER} provide the approximate confidence intervals on
#' the coefficients (fixed effects) of an saeTrafo object.
#'
#' @param object an object of type "NER".
#' @param level an optional numeric value with the confidence level for the
#' intervals. Defaults to 0.95.
#' @param parm vector of names to specify which parameters are to be given
#' confidence intervals. If \code{NULL}, all parameters are taken into account.
#' Defaults to \code{NULL}.
#' @param ... additional arguments that are not used in this method.
#' @return A matrix with rows corresponding to the parameters and columns
#' containing the lower confidence limits (lower), the
#' estimated values (est.), and upper confidence limits (upper).
#' @details The generic function \code{intervals} is imported from package nlme
#' and re-exported to make the S3-methods available, even though the nlme
#' package itself is not loaded or attached. For default documentation,
#' see \code{\link[nlme]{intervals}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{intervals}}
#'
#' @examples
#' # Example to extract confidence intervals on coefficients
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' intervals(NER_model)
#' @name intervals
#' @importFrom nlme intervals
#' @export intervals
intervals

#' @export intervals.NER
#' @export
#' @rdname intervals
intervals.NER <- function(object, level = 0.95, parm = NULL, ...) {
  throw_class_error(object, "NER")
  if (!is.null(parm)) {
    confidence_intervals <- intervals(object$model, level = level)$fixed
    subset(confidence_intervals, rownames(confidence_intervals) %in% parm)
  } else {
    intervals(object$model, level = level)$fixed
  }
}


#' Extract random effects of saeTrafo object
#'
#' Method \code{ranef.NER} extract the fixed effects from an saeTrafo object.
#'
#' @param object an object of type "NER".
#' @param ... additional arguments that are not used in this method.
#' @return A vector containing the estimated random effects at domain level is
#' returned.
#' @details The alias \code{random.effects} can also be used instead of
#' \code{ranef}. The generic function \code{ranef} is imported from package
#' \pkg{nlme} and re-exported to make the S3-methods available, even though the
#' \pkg{nlme} package itself is not loaded or attached. For default
#' documentation, see \code{\link[nlme]{random.effects}}.
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{random.effects}}
#'
#' @examples
#' # Example to extract random effects
#'
#' #Load Data
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
#'                        smp_data = eusilcA_smp, MSE = TRUE)
#'
#' ranef(NER_model)
#' @name ranef
#' @aliases random.effects
#' @importFrom nlme ranef random.effects
#' @export random.effects
#' @export ranef
ranef

#' @export ranef.NER
#' @export
#' @rdname ranef
ranef.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  ranef(object$model)
}

#' @export random.effects.NER
#' @export
#' @rdname ranef
random.effects.NER <- function(object, ...) {
  throw_class_error(object, "NER")
  ranef(object$model)
}
