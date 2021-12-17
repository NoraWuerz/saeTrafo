
mse_boot <- function(framework,
                     point_estim,
                     fixed,
                     transformation,
                     interval,
                     threshold,
                     B) {


  message('\r', "Bootstrap started")

  Y_mean_b_est <- matrix(NA, nrow = framework$N_dom_pop, ncol = B)
  Y_mean_b_orig <- matrix(NA, nrow = framework$N_dom_pop, ncol = B)

  # lapply statt for-Schleife (!!!!)
  for(b in 1: B){

    message(paste("Bootstrap", b, "started"))

    # generate random effects for each bootstrap replication
    u_d_b <- rnorm(n    = framework$N_dom_pop,
                   mean = 0,
                   sd   = sqrt(point_estim$model_par$sigmau2est)
    )

    u_di_b <- unlist(mapply(
      rep,
      x     = u_d_b,
      times = include_dom_unobs(x       = framework$n_smp,
                                obs_dom = framework$dist_obs_dom)
    ))

    # true means
    Y_mean_b_orig[,b] <- back_transformation(y              = framework$pop_mean.mat %*% point_estim$model_par$betas + u_d_b,
                                             transformation = transformation,
                                             shift          = point_estim$shift_par,
                                             lambda         = point_estim$optimal_lambda
    )

    # construct bootstrap sample
    e_di_b  <- rnorm(n    = framework$N_smp,
                     mean = 0,
                     sd   = sqrt(point_estim$model_par$sigmae2est)
    )

    Y_smp_b <- model.matrix(fixed, framework$smp_data) %*% point_estim$model_par$betas + u_di_b + e_di_b

    smp_data_b <- framework$smp_data
    smp_data_b[paste(fixed[2])] <- back_transformation(y              = Y_smp_b,
                                                       transformation = transformation,
                                                       shift          = point_estim$shift_par,
                                                       lambda         = point_estim$optimal_lambda
    )

    framework_b <- framework
    framework_b$smp_data <- smp_data_b

    # step 5: calculate mean with proposed method
    point_estim_b <- point_estim(framework      = framework_b,
                                 fixed          = fixed,
                                 transformation = transformation,
                                 threshold      = threshold,
                                 interval       = interval
    )


  }
  # auserhalb der Funktion wird nur Y_mean_b_orig und Y_mean_b_est benoetigt

  # step 6: calculate MSE
  Y_mean_b_orig_korr <- matrix(NA, nrow = framework$N_dom_pop, ncol = B)

  for(i in 1:framework$N_dom_pop){

    if(transformation == "log"){
      Y_mean_b_orig_korr[i,] <-
        Y_mean_b_orig[i,] *
        as.numeric(exp(0.5 * (framework$pop_cov.mat[i,] %*%
                                as.numeric(point_estim$model_par$betas
                                           %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est)))}
    if(transfromation == "log.shift"){
      Y_mean_b_orig_korr[i,] <-
        (Y_mean_b_orig[i,] + point_estim$optimal_lambda) *
        as.numeric(exp(0.5 * (framework$pop_cov.mat[i,]
                              %*% as.numeric(point_estim$model_par$betas
                                             %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est))) -
        point_estim$optimal_lambda}

  }
    # vorgeschlagene MSE wie in unserem Paper
    MSE_korr <- apply((point_estim_b$ind$Mean - Y_mean_b_orig_korr)^2, MARGIN = 1, FUN = mean, na.rm = T)
    return(MSE_korr)
  }




