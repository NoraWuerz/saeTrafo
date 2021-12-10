
syn_est <- function(framework, est_par, threshold) {
  data_smp <- data.frame(
    y_est = model.matrix(
      fixed,
      framework$smp_data
    )
    %*% est_par$betas,
    domain = framework$smp_domains_vec
  )
  domains <- framework$smp_domains
  area_size_names <- names(framework$pop_area_size)
  area_size <- framework$pop_area_size
  x_mean_d <- framework$pop_mean.mat %*% est_par$betas
  x_sd_d <- sqrt(framework$pop_cov.mat %*%
    as.numeric(est_par$betas %*% t(est_par$betas)))
  ###################


  area_smp <- c()
  area_smp[framework$obs_dom] <- framework$n_smp
  area_smp[!framework$obs_dom] <- 0


  # 1. Transformation der Dichte -----

  data_smp_z <- NA

  for (i in 1:framework$N_dom_pop) {
    if (length(which(data_smp$domain == area_size_names[i])) > 1) {
      data_smp_z[which(data_smp$domain == area_size_names[i])] <-
        (data_smp$y_est[which(data_smp$domain == area_size_names[i])] - mean(data_smp$y_est[which(data_smp$domain == area_size_names[i])])) / sd(data_smp$y_est[which(data_smp$domain == area_size_names[i])])
    }
    if (length(which(data_smp$domain == area_size_names[i])) == 1) {
      data_smp_z[which(data_smp$domain == area_size_names[i])] <- 0
    }
  }

  data_smp_kor_pop_x_kor_sd <- list()

  for (i in 1:framework$N_dom_pop) {
    data_smp_kor_pop_x_kor_sd[[i]] <- data_smp_z * x_sd_d[i] + x_mean_d[i]
  }

  # 3. Bestimmen der geschätzen Dichten (aus sample) und E-Wert Berechnung ----

  expectation_mod_kor_pop_2 <- c()

  for (i in 1:framework$N_dom_pop) {
    if (area_smp[i] > threshold) {
      x_tmp <- data_smp$y_est[data_smp$domain == area_size_names[i]]
      x_tmp <- (x_tmp - mean(x_tmp)) / sd(x_tmp)
      x_tmp <- x_tmp * x_sd_d[i] + x_mean_d[i]
    } else {
      x_tmp <- data_smp_kor_pop_x_kor_sd[[i]]
    }

    density_xy_mod_kor_pop <- density(x_tmp, bw = bw.SJ(x_tmp, method = "dpi"), kernel = "epanechnikov")


    expectation_mod_kor_pop_2[i] <- sfsmisc::integrate.xy(density_xy_mod_kor_pop$x, density_xy_mod_kor_pop$y * exp(density_xy_mod_kor_pop$x))
    rm(x_tmp)
  }

  # 4. Berechnen des Totalwerts ----

  E_d_density_est_kor_pop_2 <- area_size * expectation_mod_kor_pop_2

  # Ergebnisse ----

  return(E_d_density_est_kor_pop_2)
}
