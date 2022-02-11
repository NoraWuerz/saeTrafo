# Presents point estimates of an saeTrafoObject

point_saeTrafo <- function(object, indicator = "all") {

  if (is.null(object$ind)) {
    stop("No estimates in object: method point not applicable")
  }

  if (inherits(object, "fh")) {
    object$ind["Out"] <- NULL
  }

  if (any(indicator == "all") || any(indicator == "All")) {
    ind <- object$ind
    ind_name <- "All indicators"
  } else if (any(indicator == "fh") || any(indicator == "FH")) {
    ind <- object$ind[, c("Domain", "FH")]
    ind_name <- "Fay-Herriot estimates"
  } else if (any(indicator == "fh_bench") || any(indicator == "FH_Bench")) {
    ind <- object$ind[, c("Domain", "FH_Bench")]
    ind_name <- "Benchmarked Fay-Herriot estimates"
  } else if (any(indicator == "Direct") || any(indicator == "direct")) {
    ind <- object$ind[, c("Domain", "Direct")]
    ind_name <- "Direct estimates used in Fay-Herriot approach"
  } else {
    selection <- colnames(object$ind[-1]) %in% indicator
    ind <- object$ind[, c(TRUE, selection)]
    ind_name <- paste(unique(indicator), collapse = ", ")
  }

  point_saeTrafo <- list(ind = ind, ind_name = ind_name)

  return(point_saeTrafo)
}
