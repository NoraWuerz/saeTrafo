#' Summarizes an emdiObject
#'
#' Additional information about the data and model in small area estimation
#' methods and components of an emdi object are extracted. The generic function
#' summary has methods for classes "direct", "ebp" and "fh" and the returned object
#' is suitable for printing  with the \code{print}.
#' @param object an object of type "direct", "ebp" or "fh", representing point
#' and MSE estimates. Objects differ depending on the estimation method.
#' @param ... additional arguments that are not used in this method.
#' @return an object of type "summary.direct", "summary.ebp" or "summary.fh" with
#' information about the sample and population data, the usage of transformation, normality
#' tests and information of the model fit.
#' @references
#' Lahiri, P. and Suntornchost, J. (2015), Variable selection for linear mixed
#' models with applications in small area estimation, The Indian Journal of
#' Statistics 77-B(2), 312-320. \cr \cr
#' Marhuenda, Y., Morales, D. and Pardo, M.C. (2014). Information criteria for
#' Fay-Herriot model selection. Computational Statistics and Data Analysis 70,
#' 268-280. \cr \cr
#' Nakagawa S, Schielzeth H (2013). A general and simple method for obtaining R2
#' from generalized linear mixed-effects models. Methods in Ecology and Evolution,
#' 4(2), 133-142.
#' @seealso \code{\link{emdiObject}}, \code{\link{direct}}, \code{\link{ebp}},
#' \code{\link{fh}}, \code{\link[MuMIn]{r.squaredGLMM}}, \code{\link[moments]{skewness}},
#' \code{\link[moments]{kurtosis}}, \code{\link[stats]{shapiro.test}}
#' @name saeTrafo_summaries
#' @order 1
NULL
