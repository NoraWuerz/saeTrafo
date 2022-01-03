
mse_par <- function(framework,
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

    start_time = Sys.time()
    if (cpus > 1) {

      cpus <- min(cpus, parallel::detectCores())
      parallelMap::parallelStart(mode = parallel_mode,
                                 cpus = cpus, show.info = FALSE)

      if (parallel_mode == "socket") {
        parallel::clusterSetRNGStream()
      }
      parallelMap::parallelLibrary("nlme", "emdi", "stats", "sfsmisc")

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
         start_time     = start_time))

       parallelMap::parallelStop()

       # successful_bootstraps <- sum(!is.na(Y_estim_b)) / framework$N_dom_pop
       #
       # MSE <- apply((Y_estim_b - Y_mean_b_orig_korr)^2, MARGIN = 1, FUN = mean, na.rm = T)
       #
       #   return(list(MSE                   = MSE,
       #               successful_bootstraps = successful_bootstraps
       #   ))

     }else{ # hier ev. Formel ueberlegen (!!!!)

       mse_out <- NULL

     }

    }else{
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
          start_time     = start_time))

      }else{ # hier ev. Formel ueberlegen (!!!!)

        mse_out <- NULL

      }
    }
  }

  message('\r', "Bootstrap completed", "\n")
  if (.Platform$OS.type == "windows") {
    flush.console()
  }

  mses <- apply(mses, c(1,2), mean)
  mses <- data.frame(Domain = unique(framework$pop_domains_vec), mses)

  return(mses)
}



mse_bc_agg_wrapper <- function(i,
                               B,
                               framework,
                               point_estim,
                               fixed,
                               transformation,
                               interval,
                               threshold,
                               start_time) {

  tmp <- mse_bc_agg(framework      = framework,
                    point_estim    = point_estim,
                    fixed          = fixed,
                    transformation = transformation,
                    interval       = interval,
                    threshold      = threshold
  )

  if (i %% 10 == 0) {
    if (i != B) {
      delta <- difftime(Sys.time(), start_time, units = "secs")
      remaining <- (delta/i)*(B - i)
      remaining <- unclass(remaining)
      remaining <- sprintf("%02d:%02d:%02d:%02d",
                           remaining %/% 86400,  # days
                           remaining %% 86400 %/% 3600,  # hours
                           remaining %% 3600 %/% 60,  # minutes
                           remaining %% 60 %/% 1) # seconds)

      message('\r', i, " of ", B, " Bootstrap iterations completed \t Approximately ",
              remaining, " remaining \n")
      if (.Platform$OS.type == "windows") flush.console()
    }
  }
  return(tmp)
}

mse_bc_agg <- function(framework,
                       point_estim,
                       fixed,
                       transformation,
                       interval = c(-1,2),
                       threshold) {



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
  Y_mean_b_orig <- back_transformation(y              = framework$pop_mean.mat %*% point_estim$model_par$betas + u_d_b,
                                       transformation = transformation,
                                       shift          = point_estim$shift_par,
                                       lambda         = point_estim$optimal_lambda
  )

  Y_mean_b_orig_korr <- rep(NA, length = framework$N_dom_pop)

  for(i in 1:framework$N_dom_pop){

      if(transformation == "log"){
        Y_mean_b_orig_korr[i] <-
          Y_mean_b_orig[i] *
          as.numeric(exp(0.5 * (framework$pop_cov.mat[i,] %*%
                                  as.numeric(point_estim$model_par$betas
                                             %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est)))}
      if(transformation == "log.shift"){
        Y_mean_b_orig_korr[i] <-
          (Y_mean_b_orig[i] + point_estim$optimal_lambda) *
          as.numeric(exp(0.5 * (framework$pop_cov.mat[i,]
                                %*% as.numeric(point_estim$model_par$betas
                                               %*% t(point_estim$model_par$betas)) + point_estim$model_par$sigmae2est))) -
          point_estim$optimal_lambda}

  }

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
    Y_estim_b <- point_estim(framework      = framework_b,
                             fixed          = fixed,
                             transformation = transformation,
                             threshold      = threshold,
                             interval       = interval
    )$ind$Mean)

  return(list(Y_estim_b          = Y_estim_b,
              Y_mean_b_orig_korr = Y_mean_b_orig_korr))
}



