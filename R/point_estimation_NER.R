# Internal documentation -------------------------------------------------------

# Point estimation function

# This function implements the transformation of data, estimation of the nested
# error linear regression model and calculates different estimators.


point_estim <- function (framework,
                         fixed,
                         transformation,
                         threshold = threshold,
                         interval = interval,
                         L,
                         keep_data  = FALSE
) {

  # Transformation of data -------------------------------------------------------

  # Estimating the optimal parameter by optimization
  # browser()
  # Optimal parameter function returns the minimum of the optimization
  # functions from generic_opt; the minimum is the optimal lambda.
  # The function can be found in the script optimal_parameter.R

  optimal_lambda <- optimal_parameter(generic_opt    = generic_opt,
                                      fixed          = fixed,
                                      smp_data       = framework$smp_data,
                                      smp_domains    = framework$smp_domains,
                                      transformation = transformation,
                                      interval       = interval
  )

  # Data_transformation function returns transformed data and shift parameter.
  # The function can be found in the script transformation_functions.R
  # browser()
  transformation_par <- data_transformation(fixed          = fixed,
                                            smp_data       = framework$smp_data,
                                            transformation = transformation,
                                            lambda         = optimal_lambda
  )
  shift_par <- transformation_par$shift

  # Model estimation, model parameter and parameter of generating model --------

  # Estimation of the nested error linear regression model
  # See Molina and Rao (2010) p. 374
  # lme function is included in the nlme package which is imported.

  mixed_model <- nlme::lme(fixed  = fixed,
                           data   = transformation_par$transformed_data ,
                           random = as.formula(paste0("~ 1 | as.factor(",
                                                      framework$smp_domains, ")")),
                           method = "REML",
                           keep.data = keep_data)


  # Function model_par extracts the needed parameters theta from the nested
  # error linear regression model. It returns the beta coefficients (betas),
  # sigmae2est, sigmau2est and the random effect (rand_eff).

  est_par <- model_par(mixed_model = mixed_model,
                       framework   = framework,
                       fixed       = fixed,
                       transformation_par = transformation_par
  )

  # WEITER !!!
  # Schätzen von Ed schoener schreiben
  Res <- syn_est(framework = framework,
                 est_par = est_par,
                 threshold = threshold
  )

}


# All following functions are only internal ------------------------------------

# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).

model_par <- function(framework,
                      mixed_model,
                      fixed,
                      transformation_par) {
  # browser()
  if(is.null(framework$weights)) {
    # fixed parametersn
    betas <- nlme::fixed.effects(mixed_model)
    # Estimated error variance
    sigmae2est <- mixed_model$sigma^2
    # VarCorr(fit2) is the estimated random error variance
    sigmau2est <- as.numeric(nlme::VarCorr(mixed_model)[1,1])
    # Random effect: vector with zeros for all domains, filled with
    # browser()
    rand_eff <- rep(0, length(unique(framework$pop_domains_vec)))
    # random effect for in-sample domains (dist_obs_dom)
    rand_eff[framework$dist_obs_dom] <- (random.effects(mixed_model)[[1]])

    return(list(betas      = betas,
                sigmae2est = sigmae2est,
                sigmau2est = sigmau2est,
                rand_eff   = rand_eff
    )
    )
  } else {
    # fixed parameters
    betas <- nlme::fixed.effects(mixed_model)
    # Estimated error variance
    sigmae2est<-mixed_model$sigma^2
    # VarCorr(fit2) is the estimated random error variance
    sigmau2est <- as.numeric(nlme::VarCorr(mixed_model)[1,1])

    # Calculations needed for pseudo EB

    weight_sum    <- rep(0, framework$N_dom_smp)
    mean_dep      <- rep(0, framework$N_dom_smp)
    mean_indep    <- matrix(0 ,nrow = framework$N_dom_smp, ncol = length(betas))
    delta2        <- rep(0,framework$N_dom_smp)
    gamma_weight  <- rep(0,framework$N_dom_smp)
    num           <- matrix(0, nrow = length(betas), ncol = 1)
    den           <- matrix(0, nrow = length(betas), ncol = length(betas))

    for (d in 1:framework$N_dom_smp){
      domain  <- as.character(unique(framework$smp_domains_vec)[d])

      # Domain means of of the dependent variable
      dep_smp       <- transformation_par$transformed_data[[as.character(mixed_model$terms[[2]])]][
        framework$smp_domains_vec == domain]
      weight_smp    <- transformation_par$transformed_data[[as.character(framework$weights)]][
        framework$smp_domains_vec == domain]
      weight_sum[d] <- sum(weight_smp)
      indep_smp     <- model.matrix(fixed, framework$smp_data)[framework$smp_domains_vec == domain,]

      # weighted mean of the dependent variable
      mean_dep[d] <- sum(weight_smp * dep_smp) / weight_sum[d]

      # weighted means of the auxiliary information
      for (k in 1:length(betas)){
        mean_indep[d,k] <- sum(weight_smp * indep_smp[,k]) / weight_sum[d]
      }

      delta2[d]       <- sum(weight_smp^2) / (weight_sum[d]^2)
      gamma_weight[d] <- sigmau2est / (sigmau2est + sigmae2est * delta2[d])
      weight_smp_diag <- diag(weight_smp)
      dep_var_ast     <- dep_smp - gamma_weight[d] * mean_dep[d]
      indep_weight    <- t(indep_smp) %*% weight_smp_diag
      indep_var_ast   <- indep_smp - matrix(rep(gamma_weight[d] * mean_indep[d,], framework$n_smp[d]),
                                            nrow = framework$n_smp[d], byrow = TRUE)

      num <- num + (indep_weight %*% dep_var_ast)
      den <- den + (indep_weight %*% indep_var_ast)

    }


    betas    <- solve(den) %*% num
    # Random effect: vector with zeros for all domains, filled with
    rand_eff <- rep(0, length(unique(framework$pop_domains_vec)))
    # random effect for in-sample domains (dist_obs_dom)
    rand_eff[framework$dist_obs_dom] <- gamma_weight * (mean_dep - mean_indep %*% betas)


    return(list(betas      = betas,
                sigmae2est = sigmae2est,
                sigmau2est = sigmau2est,
                rand_eff   = rand_eff,
                gammaw     = gamma_weight,
                delta2     = delta2
    )
    )
  }

} # End model_par

