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

include_dom_unobs <- function(x, obs_dom) {
  if(is.vector(x) | is.numeric(x)){
    tmp <- rep(0, length(obs_dom))
    tmp[obs_dom] <- x
    return(tmp)
  }
  if(is.matrix(x)){
    tmp <- matrix(0, length(obs_dom), ncol(x))
    tmp[obs_dom,] <- x
    return(tmp)
  }
}

throw_class_error <- function(object, subclass){
  if(!inherits(object, "saeTrafo")){
    error_string <- paste0(subclass, " object has to be created by the saeTrafo package for saeTrafo methods to work.")
    stop(error_string)
  }
}


#' @importFrom stats model.frame model.matrix model.response

makeXY <- function(formula, data){
  mf <- stats::model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)

  list(y = y,
       x = x)
}


