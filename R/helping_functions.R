crbind_0 <- function(x) {
  return(cbind(0, rbind(0, x)))
}

c_1 <- function(x) {
  return(c(1, x))
}

cov_names <- function(x) {
  return(paste(
    matrix(x, nrow = length(x), ncol = length(x)),
    t(matrix(x, nrow = length(x), ncol = length(x)))
  ))
}

only.mod_vars <- function(x, var) {
  if (!is.matrix(x)) {
    return(x[var])
  } else {
    (return(x[var, var]))
  }
}
