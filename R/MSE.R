# MSE estimation function
mse <- function(framework,
                point_estim,
                fixed,
                transformation,
                interval = c(-1, 2),
                threshold,
                B,
                cpus,
                parallel_mode) {

  if (transformation == "no") {

    mse_out <- mse_prasad_rao(
      framework   = framework,
      point_estim = point_estim,
      fixed       = fixed
    )

    successful_bootstraps <- NULL
    MSE <- data.frame(
      Domain = names(framework$pop_area_size),
      Mean   = mse_out
    )

    return(list(
      MSE                   = MSE,
      successful_bootstraps = successful_bootstraps
    ))
  }

  if (transformation == "log" | transformation == "log.shift") {
    start_time <- Sys.time()

    if (cpus > 1) {
      cpus <- min(cpus, parallel::detectCores())
      parallelMap::parallelStart(
        mode = parallel_mode,
        cpus = cpus, show.info = FALSE
      )

      if (parallel_mode == "socket") {
        parallel::clusterSetRNGStream()
      }
      parallelMap::parallelLibrary(packages = c("base", "nlme", "emdi", "stats", "sfsmisc", "saeTrafo"))

      if (!is.null(framework$pop_cov)) {
        mse_out <- simplify2array(parallelMap::parallelLapply(
          xs             = seq_len(B),
          fun            = mse_bc_agg_wrapper,
          B              = B,
          framework      = framework,
          point_estim    = point_estim,
          fixed          = fixed,
          transformation = transformation,
          interval       = interval,
          threshold      = threshold,
          start_time     = start_time
        ))

        parallelMap::parallelStop()
      } else { # hier ev. Formel ueberlegen (!!!!)
        mse_out <- NULL
      }
    } else {
      if (!is.null(framework$pop_cov)) {
        mse_out <- simplify2array(lapply(
          X              = seq_len(B),
          FUN            = mse_bc_agg_wrapper,
          B              = B,
          framework      = framework,
          point_estim    = point_estim,
          fixed          = fixed,
          transformation = transformation,
          interval       = interval,
          threshold      = threshold,
          start_time     = start_time
        ))
      } else { # hier ev. Formel ueberlegen (!!!!)
        mse_out <- NULL
      }
    }

    successful_bootstraps <- sum(!is.na(mse_out)) / framework$N_dom_pop
    MSE <- data.frame(
      Domain = names(framework$pop_area_size),
      Mean   = apply(mse_out, MARGIN = 1, FUN = mean, na.rm = T)
    )

    return(list(
      MSE                   = MSE,
      successful_bootstraps = successful_bootstraps
    ))
  }

  message("\r", "Bootstrap completed", "\n")
  if (.Platform$OS.type == "windows") {
    flush.console()
  }
}

# parametric bootstrap mse estimation ------------------------------------------
# wrapper function for the parametric bootstrap showing running time information
mse_bc_agg_wrapper <- function(i,
                               B,
                               framework,
                               point_estim,
                               fixed,
                               transformation,
                               interval,
                               threshold,
                               start_time) {

  tmp <- mse_bc_agg(
    framework      = framework,
    point_estim    = point_estim,
    fixed          = fixed,
    transformation = transformation,
    interval       = interval,
    threshold      = threshold
  )

  if (i %% 10 == 0) {
    if (i != B) {
      delta <- difftime(Sys.time(), start_time, units = "secs")
      remaining <- (delta / i) * (B - i)
      remaining <- unclass(remaining)
      remaining <- sprintf(
        "%02d:%02d:%02d:%02d",
        remaining %/% 86400, # days
        remaining %% 86400 %/% 3600, # hours
        remaining %% 3600 %/% 60, # minutes
        remaining %% 60 %/% 1 # seconds
      )

      message(
        "\r", i, " of ", B, " Bootstrap iterations completed \t Approximately ",
        remaining, " remaining \n"
      )
      if (.Platform$OS.type == "windows") flush.console()
    }
  }
  return(tmp)
}

# function for the parametric bootstrap mse estimtion
# The parametric bootstrap mse estimator can be found in Wuerz et. al.
mse_bc_agg <- function(framework,
                       point_estim,
                       fixed,
                       transformation,
                       interval = c(-1, 2),
                       threshold) {

  # generate random effects for each bootstrap replication
  u_d_b <- rnorm(
    n    = framework$N_dom_pop,
    mean = 0,
    sd   = sqrt(point_estim$model_par$sigmau2est)
  )

  u_di_b <- unlist(base::mapply(
    rep,
    x     = u_d_b,
    times = include_dom_unobs(
      x       = framework$n_smp,
      obs_dom = framework$dist_obs_dom
    )
  ))

  # true means
  Y_mean_b_orig <- back_transformation(
    y              = framework$pop_mean.mat %*% point_estim$model_par$betas + u_d_b,
    transformation = transformation,
    shift          = point_estim$shift_par,
    lambda         = point_estim$optimal_lambda
  )

  Y_mean_b_orig_korr <- rep(NA, length = framework$N_dom_pop)

  for (i in 1:framework$N_dom_pop) {
    if (transformation == "log") {
      Y_mean_b_orig_korr[i] <-
        Y_mean_b_orig[i] *
          as.numeric(exp(0.5 * (framework$pop_cov.mat[i, ] %*%
            as.numeric(point_estim$model_par$betas
              %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est)))
    }
    if (transformation == "log.shift") {
      Y_mean_b_orig_korr[i] <-
        (Y_mean_b_orig[i] + point_estim$optimal_lambda) *
        as.numeric(exp(0.5 * (framework$pop_cov.mat[i, ]
        %*% as.numeric(point_estim$model_par$betas
            %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est))) -
        point_estim$optimal_lambda
    }
  }

  # construct bootstrap sample
  e_di_b <- stats::rnorm(
    n    = framework$N_smp,
    mean = 0,
    sd   = sqrt(point_estim$model_par$sigmae2est)
  )

  Y_smp_b <- stats::model.matrix(fixed, framework$smp_data) %*%
    point_estim$model_par$betas + u_di_b + e_di_b

  smp_data_b <- framework$smp_data
  smp_data_b[paste(fixed[2])] <- back_transformation(
    y              = Y_smp_b,
    transformation = transformation,
    shift          = point_estim$shift_par,
    lambda         = point_estim$optimal_lambda
  )

  framework_b <- framework
  framework_b$smp_data <- smp_data_b

  # calculate mean with proposed method
  Y_estim_b <- rep(NA, length = framework$N_dom_pop)
  try(
    Y_estim_b <- point_estim(
      framework      = framework_b,
      fixed          = fixed,
      transformation = transformation,
      threshold      = threshold,
      interval       = interval
    )$ind$Mean
  )

  return((Y_estim_b - Y_mean_b_orig_korr)^2)
}

# Prasad Rao analytical MSE estimator-------------------------------------------
# Calculates an analytical MSE estimator for the nested error regression model.
# The analytical MSE can be found in Prasad, Rao (1990).
mse_prasad_rao <- function(framework,
                           point_estim,
                           fixed,
                           transformation,
                           interval,
                           threshold) {

  out_g1 <- g1(
    sigmau2   = point_estim$model_par$sigmau2est,
    sigmae2   = point_estim$model_par$sigmae2est,
    n_smp     = include_dom_unobs(
      x       = framework$n_smp,
      obs_dom = framework$dist_obs_dom
    )
  )

  out_g2 <- g2(
    sigmau2 = point_estim$model_par$sigmau2est,
    sigmae2 = point_estim$model_par$sigmae2est,
    n_smp   = include_dom_unobs(
      x       = framework$n_smp,
      obs_dom = framework$dist_obs_dom
    ),
    Xmean = framework$pop_mean.mat,
    X     = model.matrix(fixed, framework$smp_data),
    area  = framework$smp_domains_vec
  )

  out_g3 <- g3(
    sigmau2 = point_estim$model_par$sigmau2est,
    sigmae2 = point_estim$model_par$sigmae2est,
    n_smp   = include_dom_unobs(
      x       = framework$n_smp,
      obs_dom = framework$dist_obs_dom
    ),
    X       = model.matrix(fixed, framework$smp_data),
    area     = framework$smp_domains_vec,
    pop_area = row.names(framework$pop_mean.mat)
  )

  mse1 <- out_g1 + out_g2 + 2 * out_g3
}

# Components g1, g2, g3 for the function mse_prasad_rao
g1 <- function(sigmau2,
               sigmae2,
               n_smp) {
  tmp <- ((sigmau2) / (sigmau2 + sigmae2 / n_smp)) * (sigmae2 / n_smp)
  tmp[is.nan(tmp)] <- 0
  return(tmp)
}

g2 <- function(sigmau2,
               sigmae2,
               n_smp,
               Xmean,
               X,
               area) {

  gamma_d <- sigmau2 / (sigmau2 + (sigmae2 / n_smp))

  u_area <- row.names(Xmean)
  p <- ncol(X)
  m <- nrow(Xmean)
  m_in <- length(unique(area))

  xmean <- include_dom_unobs(matrix(unlist(nlme::gapply(
        object = data.frame(X),
        FUN    = colMeans,
        groups = area
      )),
      nrow  = m_in,
      ncol  = p,
      byrow = TRUE
    ),
    obs_dom = row.names(Xmean) %in% area
  )

  X_vec <- Xmean - gamma_d * xmean

  sum_Mitte <- matrix(0, p, p)

  for (i in 1:m) {
    x_areawise <- X[area == u_area[i], ]
    V_inv <- sigmae2^(-1) * (diag(n_smp[i]) - (gamma_d[i] / n_smp[i]) *
      matrix(1, n_smp[i], n_smp[i]))

    sum_Mitte <- sum_Mitte + t(x_areawise) %*% V_inv %*% x_areawise
  }

  g2_res <- rep(NA, m)

  for (i in 1:m) {
    g2_res[i] <- t(X_vec[i, ]) %*% solve(sum_Mitte) %*% X_vec[i, ]
  }
  return(g2_res)
}

g3 <- function(sigmau2,
               sigmae2,
               n_smp,
               X,
               area,
               pop_area) {

  t <- length(n_smp)
  k <- (ncol(X) - 1)
  term <- (sum(n_smp) - t - k + 1)^(-1)

  sum_ni2_xi_xi <- matrix(0, k + 1, k + 1)

  for (i in 1:t) {
    if (n_smp[i] != 0) {
      x_areawise <- X[area == pop_area[i], ]
      smp_mean_xi <- apply(x_areawise, MARGIN = 2, FUN = mean)
      sum_ni2_xi_xi <- sum_ni2_xi_xi + n_smp[i]^2 * smp_mean_xi %*% t(smp_mean_xi)
    }
  }

  n_star <- sum(n_smp) - sum(diag(solve(t(X) %*% X) %*% sum_ni2_xi_xi))

  M <- diag(nrow(X)) - X %*% solve(t(X) %*% X) %*% t(X)
  n_star_star <- sum(diag((M %*% diag(nrow(X)) %*% t(diag(nrow(X))))^2))

  var_sig.e <- 2 * term * sigmae2^2
  var_sig.u <- 2 * n_star^(-2) * (term * (t - 1) * (sum(n_smp) - k) * sigmae2^2 +
    2 * n_star * sigmae2 * sigmau2 +
    n_star_star * sigmau2^2)
  cov_sig.e_sig.u <- -(t - 1) * n_star^(-1) * var_sig.e

  res_g3 <- n_smp^(-2) * (sigmau2 + sigmae2 / n_smp)^(-3) *
    (sigmae2^2 * var_sig.u +
       sigmau2^2 * var_sig.e -
       2 * sigmae2 * sigmau2 * cov_sig.e_sig.u)
  res_g3[is.nan(res_g3)] <- 0
  return(res_g3)
}