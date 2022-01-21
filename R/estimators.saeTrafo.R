#' Presents point, MSE and CV estimates
#'
#' Function \code{estimators} is a generic function used to present point and
#' mean squared error (MSE) estimates and calculated coefficients of variation
#' (CV).
#' @param object an object for which point and/or MSE estimates and/or
#' calculated CV's are desired.
#' @param MSE optional logical. If \code{TRUE}, MSE estimates for selected indicators
#' per domain are added to the data frame of point estimates. Defaults to
#' \code{FALSE}.
#' @param CV optional logical. If \code{TRUE}, coefficients of variation for selected
#' indicators per domain are added to the data frame of point estimates.
#' Defaults to \code{FALSE}.
#' @param ... arguments to be passed to or from other methods.
#' @return
#' The return of \code{estimators} depends on the class of its argument. The
#' documentation of particular methods gives detailed information about the
#' return of that method.
#' @export

estimators <- function(object, MSE, CV, ...) UseMethod("estimators")


#' Presents point, MSE and/or CV estimates of an saeTrafoObject
#'
#' Method \code{estimators.saeTrafo} presents point and MSE estimates for regional
#' disaggregated indicators. Coefficients of variation are calculated
#' using these estimators. This method enables to select for which indicators
#' the estimates shall be returned. The returned object is suitable for
#' printing with the \code{print.estimators.saeTrafo} method.
#' @param object an object of type "saeTrafo", representing point and,
#' if chosen, MSE estimates.
#' @param MSE optional logical. If \code{TRUE}, MSE estimates for selected indicators
#' per domain are added to the data frame of point estimates. Defaults to
#' \code{FALSE}.
#' @param CV optional logical. If \code{TRUE}, coefficients of variation for selected
#' indicators per domain are added to the data frame of point estimates.
#' Defaults to \code{FALSE}.
#' @param ... other parameters that can be passed to function \code{estimators}.
#' @return
#' The return of \code{estimators.saeTrafo} is an object of type "estimators.saeTrafo" with
#' point and/or MSE estimates and/or calculated CV's per domain obtained from
#' \code{saeTrafoObject$ind} and, if chosen, \code{saeTrafoObject$MSE}. These objects
#' contain two elements, one data frame \code{ind} and a character naming the
#' indicator or indicator group \code{ind_name}.
#' @details Objects of class "estimators.saeTrafo" have methods for following generic
#' functions: \code{head} and \code{tail} (for default documentation, see
#' \code{\link[utils]{head}}),  \code{as.matrix} (for default documentation, see
#' \code{\link[base]{matrix}}), \code{as.data.frame} (for default documentation,
#' see \code{\link[base]{as.data.frame}}), \code{subset} (for default
#' documentation, see \code{\link[base]{subset}}).
#' @seealso \code{\link{saeTrafoObject}},  \code{\link{NER_Trafo}}
#' @rdname estimators
#' @export

estimators.saeTrafo <- function(object,  MSE = FALSE, CV = FALSE, ...) {

  indicator <- c("Mean")

  #estimators_check(object = object, indicator = indicator,
  #                 MSE = MSE, CV = CV)

  # Only point estimates
  all_ind <- point_emdi(object = object, indicator = indicator)
  selected <- colnames(all_ind$ind)[-1]

  if ( MSE == TRUE || CV == TRUE ) {
    all_precisions <- mse_emdi(object = object, indicator = indicator, CV = TRUE)
    colnames(all_precisions$ind) <- paste0(colnames(all_precisions$ind), "_MSE")
    colnames(all_precisions$ind_cv) <- paste0(colnames(all_precisions$ind_cv), "_CV")
    combined <- data.frame(all_ind$ind, all_precisions$ind, all_precisions$ind_cv)
    endings <- c("","_MSE", "_CV")[c(TRUE,MSE,CV)]

    combined <- combined[,c("Domain",paste0(rep(selected,each = length(endings)),
                                            endings))]
  } else {
    combined <- all_ind$ind
  }

  estimators_emdi <- list(ind = combined, ind_name = all_ind$ind_name)

  class(estimators_emdi) <- "estimators.saeTrafo"

  return(estimators_emdi)
}

# Prints estimators.emdi objects
#' @export

print.estimators.saeTrafo <- function(x,...) {
  cat(paste0("Indicator/s: ", x$ind_name, "\n"))
  print(x$ind)
}


# Tail/head functions ----------------------------------------------------------


#' @importFrom utils head
#' @export
# CV estimators

head.estimators.saeTrafo <- function(x, n = 6L, addrownums = NULL, ...) {
  head(x$ind, n = n, addrownums = addrownums, ...)
}

#' @importFrom utils tail
#' @export

tail.estimators.saeTrafo <- function(x, n = 6L, keepnums = TRUE, addrownums = NULL, ...) {
  tail(x$ind, n = n, keepnums = keepnums, ...)
}


# Transforms estimators.saeTrafo objects into a matrix object
#' @export

as.matrix.estimators.saeTrafo <- function(x,...) {
  as.matrix(x$ind[,-1])
}

# Transforms estimators.saeTrafo objects into a dataframe object
#' @export

as.data.frame.estimators.saeTrafo <- function(x,...) {
  as.data.frame(x$ind, ...)
}

# Subsets an estimators.saeTrafo object
#' @export

subset.estimators.saeTrafo <- function(x, ...) {
  x <- as.data.frame(x)
  subset(x = x,  ...)
}



