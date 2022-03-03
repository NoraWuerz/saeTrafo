#' Fitted saeTrafoObject
#'
#' An object of class saeTrafo that represents point predictions of
#' domain-specific means. Optionally, it also contains corresponding MSE
#' estimates. Objects of these classes have methods for various generic
#' functions. See Details for more information.
#'
#' @return
#' The following components are always included in an saeTrafo object but not
#' always filled:
#' \item{\code{call}}{the function call that produced the object.}
#' \item{\code{fixed}}{for details, see \code{fixed} in
#'       \code{\link{NER_Trafo}}.}
#' \item{\code{framework}}{a list with components that describe the data setup,
#'       e.g., number of domains in the sample.}
#' \item{\code{ind}}{data frame containing estimates for the mean per domain.}
#' \item{\code{method}}{character returning the method for the estimation
#'       approach used to fit the linear mixed model and for the the optimal
#'       lambda, in our case "reml".}
#' \item{\code{model}}{list containing a  selection of model components.}
#' \item{\code{MSE}}{data frame containing MSE estimates corresponding to the
#'       mean predictions in \code{ind} per domain if MSE is selected in
#'       function call. If \code{MSE}, \code{MSE} is \code{NULL}.}
#' \item{\code{transformation}}{character or list containing information about
#'       applied transformation.}
#' \item{\code{transform_param}}{a list with two elements, \code{optimal_lambda}
#'       and \code{shift_par}, where the first contains the optimal parameter
#'       for a transformation with transformation parameter or \code{NULL} for
#'       no and log transformation and the second the potential shift parameter
#'       for the log transformation and \code{NULL} for no transformation.}
#' \item{\code{successful_bootstraps}}{a numeric returning the number of
#'       successful bootstraps. If \code{MSE = FALSE} in the function call or
#'       \code{transformation = "no"}, \code{successful_bootstraps} is
#'       \code{NULL}.}
#' @details
#' Objects of class "saeTrafo" and subclass "NER_Trafo" have following methods:
#' \code{\link[saeTrafo]{compare_pred}}, \code{\link[saeTrafo]{estimators}},
#' \code{\link[saeTrafo]{plot.saeTrafo}},
#' \code{\link[saeTrafo]{predict.NER}},
#' \code{\link[saeTrafo]{qqnorm.saeTrafo}},
#' \code{\link[saeTrafo]{compare_plot}},
#' \code{\link[saeTrafo]{getData}},
#' \code{\link[saeTrafo]{getGroups}}, \code{\link[saeTrafo]{getGroupsFormula}},
#' \code{\link[saeTrafo]{getResponse}},
#' \code{plot} (for documentation, see \code{\link[saeTrafo]{plot.saeTrafo}}),
#' \code{print}, \code{qqnorm} (for documentation, see
#' \code{\link[saeTrafo]{qqnorm.saeTrafo}}) and
#' \code{summary} (for documentation, see
#' \code{\link[saeTrafo]{summaries.saeTrafo}}),
#' \code{coef} (for default documentation, see \code{\link[stats]{coef}}),
#' \code{confint} (for default documentation, see \code{\link[stats]{confint}}),
#' \code{family} (for default documentation, see \code{\link[stats]{family}}),
#' \code{fitted} (for default documentation, see
#' \code{\link[stats]{fitted.values}}),
#' \code{\link[saeTrafo]{fixef}},
#' \code{formula} (for default documentation, see \code{\link[stats]{formula}}),
#' \code{\link[saeTrafo]{getVarCov}},
#' \code{\link[saeTrafo]{intervals}},
#' \code{logLik} (for default documentation, see \code{\link[stats]{logLik}}),
#' \code{nobs} (for default documentation, see \code{\link[stats]{nobs}}),
#' \code{\link[saeTrafo]{ranef}},
#' \code{residuals} (for default documentation, see
#' \code{\link[stats]{residuals}}),
#' \code{terms} (for default documentation, see \code{\link[stats]{terms}}),
#' \code{vcov} (for default documentation, see \code{\link[stats]{vcov}})
#' \code{sigma} (for default documentation, see \code{\link[stats]{sigma}})
#' @seealso \code{\link{NER_Trafo}}, \code{\link[nlme]{lme}},
#' \code{ \link[nlme]{lmeObject}}
#'
#' @name saeTrafoObject
NULL
