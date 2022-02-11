# Back transformation function -------------------------------------------------

back_transformation <- function(y, transformation, lambda, shift) {

  back_transformed <- if (transformation == "no") {
    no_transform_back(y = y)
  } else if (transformation == "log") {
    log_transform_back(y = y, shift = shift)
  } else if (transformation == "log.shift") {
    log_shift_opt_back(y = y, lambda = lambda)
  }

  return(y = back_transformed)
}

# Function std_data_transformation only returns a data frame with transformed
# dependent variable.

std_data_transformation <- function(fixed=fixed,
                                    smp_data,
                                    transformation,
                                    lambda) {

  y_vector <- as.matrix(smp_data[paste(fixed[2])])

  std_transformed <- if (transformation == "log.shift") {
    as.data.frame(log_shift_opt_std(y = y_vector, lambda = lambda))
  } else if (transformation == "log") {
    smp_data[paste(fixed[2])]
  } else if (transformation == "no") {
    smp_data[paste(fixed[2])]
  }

  smp_data[paste(fixed[2])] <- std_transformed
  return(transformed_data = smp_data)
}

# Transformation types ---------------------------------------------------------

# No transformation

# Transformation: no transformation
no_transform <- function(y, shift = NULL) {
  return(list(y = y, shift = NULL))
}

# Back transformation: no transformation
no_transform_back <- function(y) {
  return(y = y)
}

# Log transformation

# Transformation: log
log_transform <- function(y, shift = 0) {
  min <- min(y)
  if (min <= 0) {
    shift <- abs(min) + 1
    y <- y + shift
  }
  y <- log(y)
  return(list(y = y, shift = shift))
}

# Back transformation: log
log_transform_back <- function(y, shift = 0) {
  y <- exp(y) - shift
  return(y = y)
}


# The log-shift transformation

#  Transformation: log_shift_opt
log_shift_opt <- function(y, lambda = lambda, shift = NULL) {

  with_shift <-  function(y, lambda) {
    min <- min(y + lambda)
    if (min <= 0) {
      lambda <- lambda + abs(min) + 1
    } else {
      lambda <- lambda
    }
    return(lambda)
  }

  lambda <- with_shift(y = y, lambda = lambda) # Shift parameter

  log_trafo <- function(y, lambda = lambda) {
    y <- log(y + lambda)
    return(y)
  }
  yt <- log_trafo(y = y, lambda = lambda)
  return(list(y = yt, shift = NULL))
}

# Standardized transformation: Log_shift_opt
geometric.mean <- function(x) {
  exp(mean(log(x)))
}

log_shift_opt_std <- function(y, lambda) {

  with_shift <-  function(y, lambda) {
    min <- min(y + lambda)
    if (min <= 0) {
      lambda <- lambda + abs(min(y)) + 1
    } else {
      lambda <- lambda
    }
    return(lambda)
  }

  lambda <- with_shift(y = y, lambda = lambda) # Shift parameter

  log_trafo_std <- function(y, lambda = lambda) {
    gm <- geometric.mean(y + lambda)
    y <- gm * log(y + lambda)
    return(y)
  }

  y <- log_trafo_std(y = y, lambda = lambda)
  return(y)
}

# Back transformation: log_shift_opt
log_shift_opt_back <- function(y, lambda) {

  log_shift_opt_back <- function(y, lambda = lambda) {
    y <-  exp(y) - lambda
    return(y = y)
  }

  y <- log_shift_opt_back(y = y, lambda = lambda)
  return(y = y)
}
