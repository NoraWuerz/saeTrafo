#' Aggregates from simulated eusilc population data: domain-specific covariances
#'
#' This data contains aggregates from \code{\link{eusilcA_pop}}
#' which is based on \code{\link[simFrame]{eusilcP}} from package \pkg{simFrame}.
#'
#' @format A named list. Each element of the list contains the
#' domain-specific covariance matrice for p covariates for a specicfic domain.
#' The list is named with the respective domain name. The matrix within the list
#' has row and column names with the respective covariate names. The covarites
#' right of the ~ operator in \code{fixed} need to comprise.
#'
#' @docType data
"pop_cov"
