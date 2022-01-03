# wrapper-function

mse <- function(framework,
                point_estim,
                fixed,
                transformation,
                interval = c(-1,2),
                threshold,
                B) {

  if(transformation == "no") {

    mse_out <- mse_prasad_rao(framework      = framework,
                              point_estim    = point_estim,
                              fixed          = fixed,
                              transformation = transformation,
                              interval       = interval,
                              threshold      = threshold,
                              B              = B)

  }

  if (transformation == "log" | transformation == "log.shift") {

    if (!is.null(framework$pop_cov)) {

      mse_out <- mse_bc_agg(framework      = framework,
                            point_estim    = point_estim,
                            fixed          = fixed,
                            transformation = transformation,
                            interval       = interval,
                            threshold      = threshold,
                            B              = B)

    }else{

      mse_out <- NULL

    }
  }

}


mse_bc_agg <- function(framework,
                     point_estim,
                     fixed,
                     transformation,
                     interval = c(-1,2),
                     threshold,
                     B) {


  message('\r', "Bootstrap started")

  Y_estim_b <- matrix(NA, nrow = framework$N_dom_pop, ncol = B)
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

    # calculate mean with proposed method
    try(
    Y_estim_b[,b] <- point_estim(framework      = framework_b,
                                 fixed          = fixed,
                                 transformation = transformation,
                                 threshold      = threshold,
                                 interval       = interval
    )$ind$Mean)


  }

  # step 6: calculate MSE
  Y_mean_b_orig_korr <- matrix(NA, nrow = framework$N_dom_pop, ncol = B)

  for(i in 1:framework$N_dom_pop){

    if(transformation == "log"){
      Y_mean_b_orig_korr[i,] <-
        Y_mean_b_orig[i,] *
        as.numeric(exp(0.5 * (framework$pop_cov.mat[i,] %*%
                                as.numeric(point_estim$model_par$betas
                                           %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est)))}
    if(transformation == "log.shift"){
      Y_mean_b_orig_korr[i,] <-
        (Y_mean_b_orig[i,] + point_estim$optimal_lambda) *
        as.numeric(exp(0.5 * (framework$pop_cov.mat[i,]
                              %*% as.numeric(point_estim$model_par$betas
                                             %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est))) -
        point_estim$optimal_lambda}

  }

  successful_bootstraps <- sum(!is.na(Y_estim_b)) / framework$N_dom_pop

  MSE <- apply((Y_estim_b - Y_mean_b_orig_korr)^2, MARGIN = 1, FUN = mean, na.rm = T)

  return(list(MSE                   = MSE,
              successful_bootstraps = successful_bootstraps
  ))
}





