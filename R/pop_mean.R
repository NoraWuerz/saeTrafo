#' Aggregates from simulated eusilc population data: domain-specific means
#'
#' This data contains aggregates from \code{\link{eusilcA_pop}}
#' which is based on \code{\link[simFrame]{eusilcP}} from package
#' \pkg{simFrame}.
#'
#' @format A named list. Each element of the list contains the
#' population means for the p covariates for a specicfic domain. The list is
#' named with the respective domain name. The numeric vector within the list is
#' named with the covariate names. The covariates right of the ~ operator in
#' \code{fixed} need to comprise.
#'
#' @docType data
"pop_mean"
