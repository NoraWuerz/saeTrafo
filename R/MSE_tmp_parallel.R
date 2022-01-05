
mse_par <- function(framework,
                point_estim,
                fixed,
                transformation,
                interval = c(-1,2),
                threshold,
                B,
                cpus,
                parallel_mode) {

  if(transformation == "no") {

    mse_out <- mse_prasad_rao(framework      = framework,
                              point_estim    = point_estim,
                              fixed          = fixed,
                              transformation = transformation,
                              interval       = interval,
                              threshold      = threshold)

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

     }else{ # hier ev. Formel ueberlegen (!!!!)

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
          start_time     = start_time))

      }else{ # hier ev. Formel ueberlegen (!!!!)

        mse_out <- NULL

      }
    }

    successful_bootstraps <- sum(!is.na(mse_out)) / framework$N_dom_pop
    MSE <- data.frame(Domain = names(framework$pop_area_size),
                      Mean   = apply(mse_out, MARGIN = 1, FUN = mean, na.rm = T))

    return(list(MSE                   = MSE,
                successful_bootstraps = successful_bootstraps
    ))

  }

  message('\r', "Bootstrap completed", "\n")
  if (.Platform$OS.type == "windows") {
    flush.console()
  }
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
  Y_estim_b <- rep(NA, length = framework$N_dom_pop)
  try(
    Y_estim_b <- point_estim(framework      = framework_b,
                             fixed          = fixed,
                             transformation = transformation,
                             threshold      = threshold,
                             interval       = interval
    )$ind$Mean)

  return((Y_estim_b - Y_mean_b_orig_korr)^2)
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

mse_prasad_rao <- function(framework,
                           point_estim,
                           fixed,
                           transformation,
                           interval,
                           threshold) {


  # fuer out-of sample ist diese Komponeten 0, da dies die Unsicherheit der
  # Schaetzung von u abbildet
  g1 <- function(sigmau2,
                  sigmae2,
                  n){
    return(((sigmau2)/(sigmau2 + sigmae2/n)) * (sigmae2/n))
  }

  a1 <- g1(sigmau2 = point_estim$model_par$sigmau2est,
           sigmae2 = point_estim$model_par$sigmae2est,
           n       = include_dom_unobs(x       = framework$n_smp,
                                       obs_dom = framework$dist_obs_dom))



  ##################################
  nd <- include_dom_unobs(x = framework$n_smp, obs_dom = framework$dist_obs_dom)
  sig.u <- sqrt(point_estim$model_par$sigmau2est)
  sig.e <- sqrt(point_estim$model_par$sigmae2est)
  Xmean <- framework$pop_mean.mat # Unterschied mehrere x
  x <- model.matrix(fixed, framework$smp_data)
  saind <- framework$smp_domains_vec
  ###
  Xmean <- cbind(1,Xmean)
  X <- cbind(1,x)
  xmean <- cbind(1,tapply(x,saind, mean))
  area <- saind
  sigmae <- sig.e
  sigmau <- sig.u
  n <- nd

  g2 <- function(Xmean,X,xmean,area,sigmae,sigmau,n){
    UAREA<-unique(area)
    Mitte <- matrix(0,ncol=ncol(Xmean),nrow=ncol(Xmean))
    for(i in 1:nrow(Xmean)){ # Schleife geht alle Areas durch
      Xtmp<-X[area==UAREA[i],] # suche nun alle in-sample in area i
      Sig<-Vinv(nrow(Xtmp),sigmae,sigmau)
      # Mitte + hier wierd dann summiert
      Mitte<-Mitte+t(Xtmp)%*%(matrix(rep(Sig$nDiag * colSums(Xtmp),nrow(Xtmp)),byrow=TRUE,nrow=nrow(Xtmp),ncol=ncol(Xtmp)) + (Sig$Diag - Sig$nDiag) * Xtmp)
    }
    MitteInv<-solve(Mitte)
    lambda <- sigmau^2 / (sigmau^2 + sigmae^2/n)
    ret  <-numeric()
    for(i in 1:nrow(Xmean)){
      X1 <- Xmean[i,] - lambda[i] * xmean[i,]
      ret[i]<-t(X1)%*%MitteInv%*%(X1)
    }
    return(ret)
  }

  Vinv<-function(p,sigmae,sigmau){
    nDiagonale  <- - sigmau / (sigmae^2 + p * sigmau*sigmae)
    Diagonale <-  (sigmae + (p-1)*sigmau ) / (sigmae^2 + p * sigmau * sigmae)
    return(list(Diag=Diagonale,nDiag=nDiagonale))
  }

  g3 <- function(sigmae,sigmau,n){
    w   <-sigmae^2 + n*sigmau^2
    a   <-sum(n^2/w^2)*sum((n-1)/sigmae^4 + 1/w^2) - sum(n/w^2)^2
    Iuu <- 2/a * sum((n-1)/sigmae^4+1/w^2)
    Iee <- 2/a * sum(n^2/w^2)
    Iue <- -2/a * sum(n/w^2)
    1/n^2 * (sigmau^2+sigmae^2/n)^(-3 ) * (sigmae^4*Iuu + sigmae^4*Iee -2*sigmau^2*sigmae^2*Iue)
  }



  a2<-g2(cbind(1,Xmean),cbind(1,x),cbind(1,tapply(x,saind, mean)),saind,(sig.e),(sig.u),nd)
  a3<-g3((sig.e),(sig.u),nd)
  mse1   <- a1+a2+2*a3
}






