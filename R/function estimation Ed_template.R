# Name

syn_est <- function(framework, est_par, fixed, threshold) {

  y_est <- model.matrix(fixed, framework$smp_data) %*% est_par$betas

  x_mean_d <- framework$pop_mean.mat %*% est_par$betas
  x_sd_d <- sqrt(framework$pop_cov.mat %*%
                   as.numeric(est_par$betas %*% t(est_par$betas)))

  area_smp <- include_dom_unobs(x       = framework$n_smp,
                                obs_dom = framework$obs_dom
  )

  # get the standardised predicted values
  data_smp_z <- rep(NA, framework$N_smp)

  for (i in 1:framework$N_dom_pop) {
    pos <- framework$smp_domains_vec %in% names(framework$pop_area_size)[i]
    if (sum(pos) > 1) {
      data_smp_z[pos] <- (y_est[pos] - mean(y_est[pos])) / sd(y_est[pos])
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
                                     bw     = bw.SJ(input, method = "dpi"),
                                     kernel = "epanechnikov"
    )

    est_synthetic[i] <-
      sfsmisc::integrate.xy(x  = est_synthetic_density$x,
                            fx = est_synthetic_density$y * exp(est_synthetic_density$x)
    )
  }

  # compute and return the estimated total
  return(framework$pop_area_size * est_synthetic)

}
