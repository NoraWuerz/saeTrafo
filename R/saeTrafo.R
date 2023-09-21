#' The R Package saeTrafo for Estimating unit-level Small Area Models under
#' Transformations
#'
#' The package \pkg{saeTrafo} supports estimating regional means based on the
#' nested error regression model (\cite{Battese et al., 1988}).
#' Therefore, point estimation and mean squared error estimation
#' (\cite{Prasad and Rao, 1990}) for the classical model is offered.
#' In addition to the classical model, the logarithmic and the data-driven
#' log-shift transformation are provided.
#' The core function \code{\link{NER_Trafo}} allows several options to enter
#' population data: Either individual population data or only aggregates can be
#' entered. If full population data is given, the method of
#' \cite{Molina and Martín (2018)} is applied.
#' Compared to other small area packages, these transformations are
#' accessible in the absence of population micro-data. Only population
#' aggregates (mean values, population sizes and preferably also covariances)
#' need to be supplied. The methodology for point and mean squared error
#' estimates is described in \cite{Wuerz et al. (2022)} and is made available
#' in a user-friendly way within \pkg{saeTrafo}.
#'
#' @details
#' The estimation function is called \code{\link{NER_Trafo}}. For this
#' function, several methods are available such as
#' \code{\link{estimators.saeTrafo}} and
#' \code{\link{summaries.saeTrafo}}. For a full list, please see
#' \code{\link{saeTrafoObject}}.
#' Furthermore, functions \code{\link{map_plot}} and \code{\link{write.excel}}
#' help to visualize and export results. An overview of all currently provided
#' functions can be requested by \code{library(help=saeTrafo)}.
#'
#' @references
#' Battese, G.E., Harter, R.M. and Fuller, W.A. (1988). An Error-Components
#' Model for Predictions of County Crop Areas Using Survey and Satellite Data.
#' Journal of the American Statistical Association, Vol.83, No. 401,
#' 28-36. \cr \cr
#' Molina, I. and Martín, N. (2018). Empirical best prediction under a nested
#' error model with log transformation. The Annals of Statistics, Vol.46, No. 5,
#' 1961–1993. \cr \cr
#' Prasad, N.N., Rao, J.N.K. (1990). The estimation of the mean squared error of
#' small-area estimators. Journal of the American statistical association,
#' Vol.85, No. 409, 163-171. \cr \cr
#' Wuerz, N., Schmid, T., Tzavidis, N. (2022). Estimating regional income
#' indicators under transformations and access to limited population auxiliary
#' information. Unpublished.
#' @docType package
#' @name saeTrafo
#' @aliases saeTrafo-package
NULL
