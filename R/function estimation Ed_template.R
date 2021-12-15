
syn_est <- function(framework, est_par, fixed, threshold) {

  y_est <- model.matrix(fixed, framework$smp_data)  %*% est_par$betas

  x_mean_d <- framework$pop_mean.mat %*% est_par$betas
  x_sd_d <- sqrt(framework$pop_cov.mat %*%
    as.numeric(est_par$betas %*% t(est_par$betas)))

  area_smp <- include_dom_unobs(framework$n_smp, framework$obs_dom)

  num_area <- framework$N_dom_pop

  area_smp <- c()
  for(i in 1:num_area){
    area_smp[i] <- sum(data_smp$domain == area_size_names[i])
  }


  # Schätzung-----

  # 1. Transformation der Dichte -----

  data_smp_z <- NA

  for (i in 1:framework$N_dom_pop) {
    if (length(which(framework$smp_domains_vec == names(framework$pop_area_size)[i])) > 1) {
      data_smp_z[which(framework$smp_domains_vec == names(framework$pop_area_size)[i])] <-
        (y_est[which(framework$smp_domains_vec == names(framework$pop_area_size)[i])] -
           mean(y_est[which(framework$smp_domains_vec == names(framework$pop_area_size)[i])])) /
        sd(y_est[which(framework$smp_domains_vec == names(framework$pop_area_size)[i])])
    }
    if (length(which(framework$smp_domains_vec == names(framework$pop_area_size)[i])) == 1) {
      data_smp_z[which(framework$smp_domains_vec == names(framework$pop_area_size)[i])] <- 0
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
      x_tmp <- y_est[framework$smp_domains_vec == names(framework$pop_area_size)[i]]
      x_tmp <- (x_tmp - mean(x_tmp)) / sd(x_tmp)
      x_tmp <- x_tmp * x_sd_d[i] + x_mean_d[i]
    } else {
      x_tmp <- data_smp_kor_pop_x_kor_sd[[i]]
    }


    density_xy_mod_kor_pop <- density(x_tmp, bw = bw.SJ(x_tmp, method = "dpi"), kernel = "epanechnikov")

    if(area_smp[i] <= gewichtete_den_grenze & area_smp[i] > gewichtete_den_grenze_u){
      #print("withd different weights")
      gew_k <- ((1/length(which(data_smp$idD == area_size_names[i])) - 1/nrow(data_smp))/
                  (gewichtete_den_grenze - gewichtete_den_grenze_u))
      gewichte_tmp <- rep(1/nrow(data_smp),nrow(data_smp))
      gewichte_tmp[which(data_smp$idD == area_size_names[i])] <- gew_k * (area_smp[i] - gewichtete_den_grenze_u)+
        1/nrow(data_smp)
      gewichte_tmp <- gewichte_tmp / sum(gewichte_tmp)
    }else if(area_smp[i] > gewichtete_den_grenze){
      gewichte_tmp <- rep(1/nrow(data_smp[data_smp$idD == area_size_names[i],]), nrow(data_smp[data_smp$idD == area_size_names[i],]))
    }else{
      gewichte_tmp <- rep(1/nrow(data_smp), nrow(data_smp))
    }

    print(x_tmp)

    density_xy_mod_kor_pop     <- density(x_tmp, bw = bw.SJ(x_tmp, method = "dpi"),
                                          kernel = "epanechnikov", weights = gewichte_tmp )


    expectation_mod_kor_pop_2[i] <- sfsmisc::integrate.xy(density_xy_mod_kor_pop$x, density_xy_mod_kor_pop$y * exp(density_xy_mod_kor_pop$x))
    rm(x_tmp)
  }

  # 4. Berechnen des Totalwerts ----

  E_d_density_est_kor_pop_2 <- framework$pop_area_size * expectation_mod_kor_pop_2

  # Ergebnisse ----

  return(E_d_density_est_kor_pop_2)
}
