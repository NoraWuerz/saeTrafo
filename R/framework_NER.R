# Internal documentation -------------------------------------------------------

# The function notation defines the notational framework for NER_Trafo


framework_NER <- function(fixed, pop_area_size, pop_mean, pop_cov, pop_domains,
                          smp_data, smp_domains) {

  # Reduction of number of variables
  mod_vars <- all.vars(fixed)
  mod_vars <- mod_vars[mod_vars != as.character(fixed[2])]
  smp_vars <- c(as.character(fixed[2]), mod_vars, smp_domains, weights)
  pop_vars <- c(mod_vars, pop_domains)
  smp_data <- smp_data[, smp_vars]

# weiter
  pop_data <- pop_data[, pop_vars]


  # Deletion of NA
  if (na.rm == TRUE) {
    pop_data <- na.omit(pop_data)
    smp_data <- na.omit(smp_data)
  } else if (any(is.na(pop_data)) || any(is.na(smp_data))){
    stop('EBP does not work with missing values. Set na.rm = TRUE in function
          ebp.')
  }


  # Order of domains
  pop_data <- pop_data[order(pop_data[[pop_domains]]),]
  pop_data[[pop_domains]] <- factor(pop_data[[pop_domains]],
                                    levels = unique(pop_data[[pop_domains]]))
  pop_domains_vec <- pop_data[[pop_domains]]

  smp_data <- smp_data[order(smp_data[[smp_domains]]),]
  smp_data[[smp_domains]] <- factor(smp_data[[smp_domains]],
                                    levels = unique(pop_data[[pop_domains]]))
  smp_domains_vec <- smp_data[[smp_domains]]
  smp_domains_vec <- droplevels(smp_domains_vec)



  # Number of households in population
  N_pop <- length(pop_domains_vec)
  # Number of households in sample
  N_smp <- length(smp_domains_vec)
  # Number of out-of-sample households
  N_unobs <- N_pop - N_smp
  # Number of domains in the population
  N_dom_pop <- length(unique(pop_domains_vec))
  # Number of domains in the sample
  N_dom_smp <- length(unique(smp_domains_vec))
  # Number of out-of-sample domains
  N_dom_unobs <- N_dom_pop - N_dom_smp
  # Number of households in population per domain
  n_pop <- as.vector(table(pop_domains_vec))
  # Number of households in sample per domain
  smp_domains_vec_tmp <- as.numeric(smp_domains_vec)
  n_smp <- as.vector(table(smp_domains_vec_tmp))

  # Indicator variables that indicate if domain is in- or out-of-sample
  obs_dom <- pop_domains_vec %in% unique(smp_domains_vec)
  dist_obs_dom <- unique(pop_domains_vec) %in% unique(smp_domains_vec)



  return(list(pop_data         = pop_data,
              pop_domains_vec  = pop_domains_vec,
              smp_data         = smp_data,
              smp_domains_vec  = smp_domains_vec,
              smp_domains      = smp_domains,
              N_pop            = N_pop,
              N_smp            = N_smp,
              N_unobs          = N_unobs,
              N_dom_pop        = N_dom_pop,
              N_dom_smp        = N_dom_smp,
              N_dom_unobs      = N_dom_unobs,
              n_pop            = n_pop,
              n_smp            = n_smp,
              obs_dom          = obs_dom,
              dist_obs_dom     = dist_obs_dom,
              indicator_list   = indicator_list,
              indicator_names  = indicator_names,
              threshold        = threshold,
              weights          = weights
              )
         )
}



