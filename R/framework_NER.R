# Internal documentation -------------------------------------------------------

# The function notation defines the notational framework for NER_Trafo


framework_NER <- function(fixed, pop_area_size, pop_mean, pop_cov, pop_data,
                          pop_domains, smp_data, smp_domains) {

  # Reduction of number of variables
  mod_vars <- all.vars(fixed)
  mod_vars <- mod_vars[mod_vars != as.character(fixed[2])]
  smp_vars <- c(as.character(fixed[2]), mod_vars, smp_domains)
  smp_data <- smp_data[, smp_vars]

  if (!is.null(pop_data)) {

    fw_check_pop(
      pop_data = pop_data, mod_vars = mod_vars, pop_domains = pop_domains,
      smp_data = smp_data, fixed = fixed, smp_domains = smp_domains
    )

    pop_vars <- c(mod_vars, pop_domains)
    pop_data <- pop_data[, pop_vars]

  } else {

    fw_check_agg(
       pop_area_size = pop_area_size, pop_mean = pop_mean, pop_cov = pop_cov,
       mod_vars = mod_vars, smp_data = smp_data, fixed = fixed,
       smp_domains = smp_domains
     )

    # hier weiter !!!
    # Means aufarbeiten
    # 1. means ordnen mod_vars und nur die aus mod_vars verwenden
    #lapply

    # matrix erstellen
    pop_mean.mat <- matrix(unlist(lapply(pop_mean, c_1)), ncol = length(mod_vars), byrow = TRUE)
    row.names(pop_mean.mat) <- names(pop_mean)
    colnames(pop_mean.mat) <- c("intercept", mod_vars)

    # Cov aufbereiten
    pop_cov.mat <- matrix(unlist(lapply(pop_cov, crbind_0)),
                          ncol = (length(mod_vars) + 1)^2, byrow = TRUE)
    row.names(pop_cov.mat) <- names(pop_cov)
    colnames(pop_cov.mat) <- cov_names(c("intercept", mod_vars))

  }



  # Order of domains
  pop_data <- pop_data[order(pop_data[[pop_domains]]), ]
  pop_data[[pop_domains]] <- factor(pop_data[[pop_domains]],
    levels = unique(pop_data[[pop_domains]])
  )
  pop_domains_vec <- pop_data[[pop_domains]]

  smp_data <- smp_data[order(smp_data[[smp_domains]]), ]
  smp_data[[smp_domains]] <- factor(smp_data[[smp_domains]],
    levels = unique(pop_data[[pop_domains]])
  )
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



  return(list(
    pop_data = pop_data,
    pop_domains_vec = pop_domains_vec,
    smp_data = smp_data,
    smp_domains_vec = smp_domains_vec,
    smp_domains = smp_domains,
    N_pop = N_pop,
    N_smp = N_smp,
    N_unobs = N_unobs,
    N_dom_pop = N_dom_pop,
    N_dom_smp = N_dom_smp,
    N_dom_unobs = N_dom_unobs,
    n_pop = n_pop,
    n_smp = n_smp,
    obs_dom = obs_dom,
    dist_obs_dom = dist_obs_dom
  ))
}
