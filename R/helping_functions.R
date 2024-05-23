# Auxiliary functions for the package saeTrafo ---------------------------------

# add a 0 row and 0 column to a matrix
crbind_0 <- function(x) {
  return(cbind(0, rbind(0, x)))
}

# add a 1 to a vector
c_1 <- function(x) {
  return(c(1, x))
}

# generate covariance names from covariate name vector
cov_names <- function(x) {
  return(paste(
    matrix(x, nrow = length(x), ncol = length(x)),
    t(matrix(x, nrow = length(x), ncol = length(x)))
  ))
}

# returns matrix or vector only with entries which match to the names in var
only_mod_vars <- function(x, var) {
  if (!(is.matrix(x) | is.data.frame(x))) {
    return(x[var])
  } else {
    (return(x[var, var]))
  }
}

# for vector or matices unobserved domains are included in the respective object
include_dom_unobs <- function(x, obs_dom) {
  if (is.matrix(x)) {
    tmp <- matrix(0, length(obs_dom), ncol(x))
    tmp[obs_dom, ] <- x
    return(tmp)
  } else {
    tmp <- rep(0, length(obs_dom))
    tmp[obs_dom] <- x
    return(tmp)
  }

}

# check if the object is of class saeTrafo
throw_class_error <- function(object, subclass) {
  if (!inherits(object, "saeTrafo")) {
    error_string <- paste0(subclass, " object has to be created by the ",
                           "saeTrafo package for saeTrafo methods to work.")
    stop(error_string)
  }
}

# Extract response and covariates from data
#' @importFrom stats model.frame model.matrix model.response
makeXY <- function(formula, data) {
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)

  list(y = y,
       x = x)
}

# Extract variable names from all domains
all_in_domain <- function(x, mod_vars) {
  if (is.matrix(x)) {
    all(mod_vars %in% row.names(x)) && all(mod_vars %in% colnames(x))
  } else {
    all(mod_vars %in% names(x))
  }
}

# likelihoods function for plots
likelihoods_f <- function(lam, fixed, smp_data, smp_domains, transformation) {
  result <- NULL
  try(result <- -as.numeric(generic_opt(lambda         = lam,
                                        fixed          = fixed,
                                        smp_data       = smp_data,
                                        smp_domains    = smp_domains,
                                        transformation = transformation
                 )),
      silent = TRUE)
  if (is.null(result)) {result <- NA}
  result
}
