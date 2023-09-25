#' Nested error regression Model under transformations
#'
#' Function \code{NER_Trafo} estimates small area means based on the
#' (transformed) nested error regression (NER) model
#' (\cite{Battese et al., 1988}).
#' In contrast to the empirical best predictor of \cite{Molina and Rao (2010)},
#' which is implemented in the package \pkg{emdi} (\code{\link[emdi]{ebp}}), no
#' unit-level population data are required.\cr \cr
#' \code{NER_Trafo} supports the log as well as the data-driven log-shift
#' transformation. Especially for skewed variables, (data-driven)
#' transformations are useful to meet the model assumptions for the error terms.
#' If a transformation is chosen and aggregates (means and covariance) are
#' simultaneously provided for the population, point estimates are produced by
#' the method of \cite{Wuerz et al. (2022)}, which uses kernel density
#' estimation to resolve the issue of not having access to population
#' micro-data.
#' In the case that population data are available at unit-level and the log or
#' log-shift transformation is selected, the bias-correction of
#' \cite{Berg and Chandra (2014)} and \cite{Molina and Martín (2018)} is
#' applied. For this data situation, more methods and options are provided in
#' the package \pkg{emdi}.
#' If only population means are available and the log or log-shift
#' transformation is selected, a bias-correction due to
#' the transformation is added but for the lack of access to population data no
#' correction is available. Therefore, a part of the bias is
#' disregarded.\cr \cr
#' Additionally, analytically mean squared errors (MSE) are calculated in the
#' case of no transformation following \cite{Prasad and Rao (1990)}.
#' For the log and log-shift transformation, a parametric bootstrap procedure
#' proposed by \cite{Wuerz et al. (2022)} following
#' \cite{Gonzalez-Manteiga et al. (2008)} is applied. Please note that this can
#' only be determined if covariance data are also provided.
#' If population data is available on unit-level a bootstrap procedure as
#' described in \cite{Molina and Martín (2018)} is applied.
#'
#'
#' @param fixed a two-sided linear formula object describing the
#' fixed-effects part of the nested error linear regression model with the
#' dependent variable on the left of a ~ operator and the explanatory
#' variables on the right, separated by + operators. The argument corresponds
#' to the argument \code{fixed} in function \code{\link[nlme]{lme}}.
#' @param pop_area_size a named numeric vector containing the number of
#' individuals within each domain. This numeric vector is named with the
#' domain names.
#' @param pop_mean a named list. Each element of the list contains the
#' population means for the p covariates for a specicfic domain. The list is
#' named with the respective domain name. The numeric vector within the list is
#' named with the covariate names. The covariates right of the ~ operator in
#' \code{fixed} need to comprise.
#' @param pop_cov a named list. Each element of the list contains the
#' domain-specific covariance matrice for p covariates for a specicfic domain.
#' The list is named with the respective domain name. The matrix within the list
#' has row and column names with the respective covariate names. The covariates
#' right of the ~ operator in \code{fixed} need to comprise. If \code{pop_cov}
#' is not available only a bias-correction due to the transformation is
#' added but for the lack of access to population data a correction is not
#' possible. Additionally, no MSE could be determined.
#' @param pop_data a data frame that needs to comprise the variables
#' named on the right of the ~ operator in \code{fixed}, i.e. the explanatory
#' variables, and \code{pop_domains}.
#' Please note, if population data is available other methods using unit-level
#' population data, like \code{\link[emdi]{ebp}}, could be applied.
#' @param pop_domains a character string containing the name of a variable that
#' indicates domains in the population data. The variable can be numeric or
#' a factor but needs to be of the same class as the variable named in
#' \code{smp_domains}. Only needed if population data are given.
#' @param smp_data a data frame that needs to comprise all variables named in
#' \code{fixed} and \code{smp_domains}.
#' @param smp_domains a character string containing the name of a variable
#' that indicates domains in the sample data. The variable can be numeric or a
#' factor but needs to be of the same class as the variable named in
#' \code{pop_domains}.
#' @param threshold a numeric value indicating the threshold for using pooled
#' domain data (for domains with sample sizes below the threshold) or non pooled
#' domain data (for domains with sample sizes above the threshold) for the
#' density estimation within the approach of \cite{Wuerz et al. (2022)}.
#' Defaults to 30.
#' @param transformation a character string. Three different transformation
#' types for the dependent variable can be chosen (i) no transformation ("no");
#' (ii) log transformation ("log"); (iii) Log-Shift transformation
#' ("log.shift"). Defaults to \code{"log.shift"}.
#' @param interval a string equal to 'default' or a numeric vector containing a
#' lower and upper limit determining an interval for the estimation of the
#' optimal parameter for the log-shift transformation. The interval is passed to
#' function  \code{\link[stats]{optimize}} for the optimization. Defaults to
#' an interval based on the range of y. If the convergence fails, it is often
#' advisable to choose a smaller more suitable interval. For right skewed
#' distributions, the negative values may be excluded.
#' @param MSE optional logical. If \code{TRUE}, MSE estimates are provided.
#' Defaults to \code{FALSE}.
#' @param B a number determining the number of bootstrap replications in the
#' parametric bootstrap approach. The number must be greater than 1. Defaults to
#' 50. For practical applications, values larger than 200 are recommended.
#' @param seed an integer to set the seed for the random number generator. For
#' the usage of random number generation, see Details. If seed is set to
#' \code{NULL}, seed is chosen randomly. Defaults to \code{123}.
#' @param parallel_mode modus of parallelization, defaults to an automatic
#' selection of a suitable mode, depending on the operating system, if the
#' number of \code{cpus} is chosen higher than 1. For details, see
#' \code{\link[parallelMap]{parallelStart}}.
#' @param cpus number determining the kernels that are used for the
#' parallelization. Defaults to 1. For details, see
#' \code{\link[parallelMap]{parallelStart}}.
#' @return An object of class "NER", "saeTrafo" that provides estimators for
#' regional means optionally corresponding MSE estimates. Several generic
#' functions have methods for the returned object. For a full list and
#' descriptions of the components of objects of class "saeTrafo", see
#' \code{\link{saeTrafoObject}}.
#' @details For the parametric bootstrap and the density estimation
#' approach random number generation is used. Thus, a seed is set by the
#' argument \code{seed}. \cr \cr
#' @references
#' Battese, G.E., Harter, R.M. and Fuller, W.A. (1988). An Error-Components
#' Model for Predictions of County Crop Areas Using Survey and Satellite Data.
#' Journal of the American Statistical Association, Vol.83, No. 401,
#' 28-36. \cr \cr
#' Berg, E. and Chandra, H. (2014). Small area prediction for a unit-level
#' lognormal model. Computational Statistics & Data Analysis, Vol.78,
#' 159–175.\cr \cr
#' González-Manteiga, W., Lombardía, M. J., Molina, I., Morales, D. and
#' Santamaría, L. (2008). Analytic and bootstrap approximations of prediction
#' errors under a multivariate Fay–Herriot model. Computational Statistics &
#' Data Analysis, Vol. 52, No. 12, 5242-5252. \cr \cr
#' Molina, I. and Martín, N. (2018). Empirical best prediction under a nested
#' error model with log transformation. The Annals of Statistics, Vol.46, No. 5,
#' 1961–1993. \cr \cr
#' Molina, I. and Rao, J.N.K. (2010). Small area estimation of poverty
#' indicators. The Canadian Journal of Statistics, Vol. 38, No.3,
#' 369-385. \cr \cr
#' Prasad, N.N., Rao, J.N.K. (1990). The estimation of the mean squared error of
#' small-area estimators. Journal of the American statistical association,
#' Vol. 85, No. 409, 163-171. \cr \cr
#' Wuerz, N., Schmid, T., and Tzavidis, N. (2022) Estimating regional income
#' indicators under transformations and access to limited population auxiliary
#' information. Journal of the Royal Statistical Society: Series A
#' (Statistics in Society), Vol. 185, No. 4, 1679-1706.
#' @examples
#'
#' # Examples for (transformed) nested error regression model
#'
#' # Load Data
#' data("eusilcA_pop")
#' data("eusilcA_smp")
#' data("pop_area_size")
#' data("pop_mean")
#' data("pop_cov")
#'
#' # formula object for all examples
#' formula <- eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben +
#'                       age_ben + surv_ben + sick_ben + dis_ben + rent +
#'                       fam_allow + house_allow + cap_inv + tax_adj
#'
#' # For all four examples, no MSEs/variances are determined in order to avoid
#' # long run times. These can be obtained with MSE = TRUE.
#'
#' # Example 1: No transformation - classical NER
#' NER_model_1 <- NER_Trafo(fixed = formula, transformation = "no",
#'                          smp_domains = "district", smp_data = eusilcA_smp,
#'                          pop_area_size = pop_area_size, pop_mean = pop_mean)
#'
#' # Example 2: Log-shift transformation and population aggregates
#' # (means and covariances) with changed threshold
#' NER_model_2 <- NER_Trafo(fixed = formula,
#'                          smp_domains = "district", smp_data = eusilcA_smp,
#'                          pop_area_size = pop_area_size, pop_mean = pop_mean,
#'                          pop_cov = pop_cov, threshold = 50)
#'
#' # Example 3: Log-shift transformation and population data
#' # A bias-corrections which need unit-level population data are applied
#' \donttest{NER_model_3 <- NER_Trafo(fixed = formula,
#'                          smp_domains = "district", smp_data = eusilcA_smp,
#'                          pop_data = eusilcA_pop, pop_domains = "district")
#'}
#'
#' # Example 4: Log-shift transformation and population aggregates
#' # (only means (!) - Therefore, no MSE estimation is available, bias is
#' # disregarded)
#' NER_model_4 <- NER_Trafo(fixed = formula,
#'                          smp_domains = "district", smp_data = eusilcA_smp,
#'                          pop_area_size = pop_area_size, pop_mean = pop_mean)
#'
#' @seealso \code{\link{saeTrafoObject}}, \code{\link[nlme]{lme}},
#' \code{\link{estimators.saeTrafo}},  \code{\link{plot.saeTrafo}},
#' \code{\link{summaries.saeTrafo}}
#' @export
#' @importFrom emdi data_transformation

NER_Trafo <- function(fixed,
                      pop_area_size = NULL,
                      pop_mean = NULL,
                      pop_cov = NULL,
                      pop_data = NULL,
                      pop_domains = NULL,
                      smp_data,
                      smp_domains,
                      threshold = 30,
                      B = 50,
                      transformation = "log.shift",
                      interval = "default",
                      MSE = FALSE,
                      parallel_mode = ifelse(grepl("windows", .Platform$OS.type),
                                             "socket", "multicore"),
                      cpus = 1,
                      seed = 123) {

  NER_check1(fixed, pop_data, pop_domains, smp_data, smp_domains, pop_mean,
             pop_area_size)
  NER_check2(threshold, transformation, interval, MSE, B, cpus, seed)

  # Save function call ---------------------------------------------------------
  call <- match.call()
  if (inherits(call$fixed, "name")) {
    call$fixed <- fixed
  }

  # Data manipulation and notational framework ---------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # The function framework_NER can be found in script framework_NER.R
  framework <- framework_NER(pop_area_size = pop_area_size,
                             pop_mean      = pop_mean,
                             pop_cov       = pop_cov,
                             pop_data      = pop_data,
                             pop_domains   = pop_domains,
                             smp_data      = smp_data,
                             smp_domains   = smp_domains,
                             fixed         = fixed
  )


  # Point Estimation -----------------------------------------------------------
  # The function point_estim can be found in script point_estimation.R
  point_estim <- point_estim(framework      = framework,
                             fixed          = fixed,
                             transformation = transformation,
                             threshold      = threshold,
                             interval       = interval,
                             keep_data      = TRUE
  )

  # MSE Estimation -------------------------------------------------------------

  if (MSE == TRUE) {

    if (is.null(pop_cov) &&  is.null(pop_data) && MSE == TRUE &&
        transformation != "no") {
      stop(paste("No MSE estimator available. For MSE estimation a covariance",
                 "matrices or population data are needed otherwise the Prasad",
                 "Rao MSE (transformation == \"no\") is avaliable."))
    }

    # The function parametric_bootstrap can be found in script mse_estimation.R

    mse_estimates <- mse(framework      = framework,
                         point_estim    = point_estim,
                         fixed          = fixed,
                         transformation = transformation,
                         interval       = interval,
                         threshold      = threshold,
                         B              = B,
                         cpus           = cpus,
                         parallel_mode  = parallel_mode
    )

    NER_out <- list(
      ind                   = point_estim$ind,
      MSE                   = mse_estimates$MSE,
      transform_param       = point_estim[c("optimal_lambda", "shift_par")],
      model                 = point_estim$model,
      framework             = framework[c("N_dom_unobs", "N_dom_smp", "N_smp",
                                          "N_pop", "smp_domains", "smp_data",
                                          "smp_domains_vec", "pop_area_size",
                                          "pop_mean.mat",  "pop_cov.mat")],
      transformation        = transformation,
      method                = "reml",
      fixed                 = fixed,
      call                  = call,
      successful_bootstraps = mse_estimates$successful_bootstraps
    )
  } else {
    NER_out <- list(
      ind                   = point_estim$ind,
      MSE                   = NULL,
      transform_param       = point_estim[c("optimal_lambda", "shift_par")],
      model                 = point_estim$model,
      framework             = framework[c("N_dom_unobs", "N_dom_smp", "N_smp",
                                          "N_pop", "smp_domains", "smp_data",
                                          "smp_domains_vec", "pop_area_size",
                                          "pop_mean.mat",  "pop_cov.mat")],
      transformation        = transformation,
      method                = "reml",
      fixed                 = fixed,
      call                  = call,
      successful_bootstraps = NULL
    )
  }

  class(NER_out) <- c("saeTrafo", "NER")
  return(NER_out)
}
