# Point estimation function for saeTrafo objects ------------------------------

# This function implements the transformation of data, estimation of the nested
# error linear regression model and calculates different estimators.
#' @importFrom emdi data_transformation
#' @importFrom nlme lme random.effects
#' @importFrom stats model.matrix reformulate
#'

point_estim <- function(framework,
                        fixed,
                        transformation,
                        threshold = threshold,
                        interval = interval,
                        keep_data = FALSE) {

  # Transformation of data -----------------------------------------------------

  # Estimating the optimal parameter by optimization
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
  # The data_transformation is in the package emdi which is imported.

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
  mixed_model <- lme(fixed     = fixed,
                     data      = transformation_par$transformed_data,
                     random    = as.formula(paste0(
                       "~ 1 | as.factor(", framework$smp_domains, ")")),
                     method    = "REML",
                     keep.data = keep_data
  )


  # Function model_par extracts the needed parameters theta from the nested
  # error linear regression model. It returns the beta coefficients (betas),
  # sigmae2est, sigmau2est and the random effect (rand_eff).
  est_par <- model_par(mixed_model        = mixed_model,
                       framework          = framework,
                       fixed              = fixed,
                       transformation_par = transformation_par
  )


  # Compute area-specific means ------------------------------------------------

  ind <- data.frame(Domain = names(framework$pop_area_size), Mean = NA)

  # if transformation == "no" compute the BHF (Battese et al., 1988)
  if (transformation == "no") {

    value_in_sample <-
      tapply(X     = framework$smp_data[as.character(fixed[2])][, 1],
             INDEX = framework$smp_domains_vec,
             FUN   = mean
    )

    gamma_est_in <- est_par$sigmau2est /
      (est_par$sigmau2est + est_par$sigmae2est / framework$n_smp)

    ind$Mean[framework$dist_obs_dom] <-
      gamma_est_in *
      (value_in_sample + (framework$pop_mean.mat[framework$dist_obs_dom, ]
                          %*% est_par$betas)[, 1] -
        tapply(X     = (model.matrix(fixed, framework$smp_data) %*%
                          est_par$betas)[, 1],
               FUN   = mean,
               INDEX = framework$smp_domains_vec
        )) +
      (1 - gamma_est_in) *
      (framework$pop_mean.mat[framework$dist_obs_dom, ] %*% est_par$betas)[, 1]

    ind$Mean[!framework$dist_obs_dom] <-
      framework$pop_mean.mat[!framework$dist_obs_dom, ] %*% est_par$betas

  }

  # if transformation == "log" or "log.shift" and ...
  # ... pop_cov is available compute bc-agg from Wuerz et. al.
  # ... pop_cov is not available compute bc-naive-agg see Wuerz et. al.
  # ... pop_data is available compute bc

  if (transformation == "log" | transformation == "log.shift") {

    n_smp_long <- include_dom_unobs(x       = framework$n_smp,
                                    obs_dom = framework$dist_obs_dom
    )

    rand_eff_long <- include_dom_unobs(
      x       = est_par$rand_eff[framework$dist_obs_dom],
      obs_dom = framework$dist_obs_dom
    )

    gamma_est_d <- est_par$sigmau2est /
      (est_par$sigmau2est + est_par$sigmae2est / n_smp_long)
    bc_d <- (est_par$sigmau2est * (1 - gamma_est_d) + est_par$sigmae2est) / 2


    if (is.null(framework$pop_data)) {

      if (!is.null(framework$pop_cov)) {

        synthetic <- syn_est(framework = framework,
                             est_par   = est_par,
                             fixed     = fixed,
                             threshold = threshold
        )

        if (transformation == "log") {
          ind$Mean <- c(1 / framework$pop_area_size *
            (synthetic * exp(rand_eff_long + bc_d)) - shift_par)
        }
        if (transformation == "log.shift") {
          ind$Mean <- c(1 / framework$pop_area_size *
            (synthetic * exp(rand_eff_long + bc_d)) - optimal_lambda)
        }
      } else {
        est_dr <- (framework$pop_mean.mat %*% est_par$betas)[, 1] +
          rand_eff_long + bc_d
        est_ds <- include_dom_unobs(
          x      = tapply(X = framework$smp_data[as.character(fixed[2])][, 1],
                          INDEX = framework$smp_domains_vec,
                          FUN   = mean
                    ),
          obs_dom = framework$dist_obs_dom
        )

        ind$Mean <- c(1 / framework$n_pop *
          (n_smp_long * est_ds + (framework$n_pop - n_smp_long) *
              back_transformation(y              = est_dr,
                                  transformation = transformation,
                                  shift          = shift_par,
                                  lambda         = optimal_lambda
              )
           ))
      }

    } else {

      mod_vars <- all.vars(fixed)
      mod_vars <- mod_vars[mod_vars != as.character(fixed[2])]

      est_pop <- (model.matrix(reformulate(mod_vars), framework$pop_data)
                 %*% est_par$betas)[, 1] +
        unlist(mapply(FUN   = rep,
                      x     = rand_eff_long + bc_d,
                      times = table(framework$pop_domains_vec)
        ))
      est_pop_back <- back_transformation(y              = est_pop,
                                          transformation = transformation,
                                          shift          = shift_par,
                                          lambda         = optimal_lambda
      )

      est_ds <- include_dom_unobs(
        x      = tapply(X     = framework$smp_data[as.character(fixed[2])][, 1],
                        INDEX = framework$smp_domains_vec,
                        FUN   = mean
        ),
        obs_dom = framework$dist_obs_dom
      )

      est_ds_est <- (model.matrix(reformulate(mod_vars), framework$smp_data)
                  %*% est_par$betas)[, 1] +
        unlist(mapply(FUN   = rep,
                      x     = est_par$rand_eff[framework$dist_obs_dom] +
                        bc_d[framework$dist_obs_dom],
                      times = table(framework$smp_domains_vec)
        ))
      est_ds_est_back <- back_transformation(y              = est_ds_est,
                                             transformation = transformation,
                                             shift          = shift_par,
                                             lambda         = optimal_lambda
      )

      ind$Mean <- c(1 / framework$n_pop *
        (n_smp_long * est_ds
         +
           framework$n_pop * tapply(X     = est_pop_back,
                                    INDEX = framework$pop_domains_vec,
                                    FUN   = mean)
         -
           n_smp_long *
            include_dom_unobs(
              x       = as.numeric(tapply(X     = est_ds_est_back,
                                          INDEX = framework$smp_domains_vec,
                                          FUN   = mean)),
              obs_dom = framework$dist_obs_dom
            )
        ))
    }
  }

  return(list(ind            = ind,
              optimal_lambda = optimal_lambda,
              shift_par      = shift_par,
              model_par      = est_par,
              model          = mixed_model
  ))
}

# Functions to extract and calculate model parameter----------------------------

# Function model_par extracts the needed parameters theta from the nested
# error linear regression model. It returns the beta coefficients (betas),
# sigmae2est, sigmau2est and the random effect (rand_eff).
#' @importFrom nlme fixed.effects VarCorr

model_par <- function(framework, mixed_model, fixed, transformation_par) {

    # fixed parametersn
    betas <- fixed.effects(mixed_model)
    # Estimated error variance
    sigmae2est <- mixed_model$sigma^2
    # VarCorr(fit2) is the estimated random error variance
    sigmau2est <- as.numeric(VarCorr(mixed_model)[1, 1])
    # Random effect: vector with zeros for all domains, filled with
    rand_eff <- rep(0, length(unique(framework$pop_domains_vec)))
    # random effect for in-sample domains (obs_dom)
    rand_eff[framework$dist_obs_dom] <- (random.effects(mixed_model)[[1]])

    return(list(betas      = betas,
                sigmae2est = sigmae2est,
                sigmau2est = sigmau2est,
                rand_eff   = rand_eff
    ))
}

# Function for synthetic part estimation ---------------------------------------

# Function for small area estimation specific kernel-density estimation
# method for the synthetic part see Wuerz et. al.
#' @importFrom stats density bw.SJ model.matrix sd
#' @importFrom sfsmisc integrate.xy

syn_est <- function(framework, est_par, fixed, threshold) {

  y_est <- model.matrix(fixed, framework$smp_data) %*% est_par$betas

  x_mean_d <- framework$pop_mean.mat %*% est_par$betas
  x_sd_d <- sqrt(framework$pop_cov.mat %*%
                   as.numeric(est_par$betas %*% t(est_par$betas)))

  area_smp <- include_dom_unobs(x       = framework$n_smp,
                                obs_dom = framework$dist_obs_dom
  )

  # get the standardised predicted values
  data_smp_z <- rep(NA, framework$N_smp)

  for (i in 1:framework$N_dom_pop) {
    pos <- framework$smp_domains_vec %in% names(framework$pop_area_size)[i]
    if (sum(pos) > 1) {
      data_smp_z[pos] <- (y_est[pos] - mean(y_est[pos])) / sd(y_est[pos])
      if (any(((y_est[pos] - mean(y_est[pos]))) == 0)) {
        warning(paste("Generated Bootstrap-sample for area", i,
                      "contain the following inputs for the dependent variable",
                      paste(y_est[pos], collapse = " "),
                      "and therefore y - mean(y) = 0 and thus",
                      "y - mean(y) / sd(y) == Nan"))
      }
    }
    if (sum(pos) == 1) {
      data_smp_z[pos] <- 0
    }
  }

  # adjuste all standardised predicted values and use all or only the data from
  # one area according to the respective sample size
  est_synthetic <- c()

  for (i in 1:framework$N_dom_pop) {

    if (area_smp[i] > threshold) {
      pos <- framework$smp_domains_vec %in% names(framework$pop_area_size)[i]
      input <- data_smp_z[pos] * x_sd_d[i] + x_mean_d[i]
    } else {
      input <- data_smp_z * x_sd_d[i] + x_mean_d[i]
    }

    est_synthetic_density <- density(x      = input,
                                     bw     = bw.SJ(x = input, method = "dpi"),
                                     kernel = "epanechnikov"
    )

    est_synthetic[i] <- integrate.xy(x  = est_synthetic_density$x,
                                     fx = est_synthetic_density$y *
                                       exp(est_synthetic_density$x)
    )
  }

  # compute and return the estimated total
  return(framework$pop_area_size * est_synthetic)
}
