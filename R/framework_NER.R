# Internal documentation -------------------------------------------------------
# The function notation defines the notational framework for NER_Trafo

#' @importFrom stats cov model.matrix

framework_NER <- function(fixed, pop_area_size, pop_mean, pop_cov, pop_data,
                          pop_domains, smp_data, smp_domains) {

  # Reduction of number of variables
  mod_vars <- all.vars(fixed)
  mod_vars <- mod_vars[mod_vars != as.character(fixed[2])]
  smp_vars <- c(as.character(fixed[2]), mod_vars, smp_domains)


  # no population data available
  if (!is.null(pop_data)) {

    fw_check_pop(
      pop_data = pop_data, mod_vars = mod_vars, pop_domains = pop_domains,
      smp_data = smp_data, fixed = fixed, smp_domains = smp_domains
    )

    smp_data <- smp_data[, smp_vars]
    pop_area_size <- as.numeric(table(pop_data[pop_domains]))
    names(pop_area_size) <- names(table(pop_data[pop_domains]))

    N_dom_pop <- length(pop_area_size)

    pop_mean.mat <- matrix(data     = NA,
                           nrow     = length(pop_area_size),
                           ncol     = length(mod_vars) + 1,
                           dimnames = list(names(pop_area_size),
                                           c("intercept", mod_vars))
    )
    pop_cov.mat <- matrix(data     = NA,
                          nrow     = length(pop_area_size),
                          ncol     = (length(mod_vars) + 1)^2,
                          dimnames = list(names(pop_area_size),
                                          cov_names(c("intercept", mod_vars)))
    )
    for (i in 1:N_dom_pop) {
      pos <- pop_data[pop_domains] == names(pop_area_size)[i]
      pop_mean.mat[i, ] <- apply(X      = model.matrix(fixed, pop_data[pos, ]),
                                 MARGIN = 2,
                                 FUN    = mean
      )
      pop_cov.mat[i, ] <- c(cov(model.matrix(fixed, pop_data[pos, ]),
                                model.matrix(fixed, pop_data[pos, ]))
      )
    }

    pop_vars <- c(mod_vars, pop_domains)
    pop_data <- pop_data[, pop_vars]

    pop_data <- pop_data[order(pop_data[[pop_domains]]), ]
    pop_data[[pop_domains]] <- factor(x      = pop_data[[pop_domains]],
                                      levels = unique(pop_data[[pop_domains]])
    )
    pop_domains_vec <- pop_data[[pop_domains]]

    smp_data <- smp_data[order(smp_data[[smp_domains]]), ]
    smp_data[[smp_domains]] <- factor(x      = smp_data[[smp_domains]],
                                      levels = unique(pop_data[[pop_domains]])
    )
    smp_domains_vec <- smp_data[[smp_domains]]
    smp_domains_vec <- droplevels(smp_domains_vec)

    # Number of households in population
    N_pop <- length(pop_domains_vec)
    # Number of households in population per domain
    n_pop <- as.vector(table(pop_domains_vec))

    # Indicator variables that indicate if domain is in- or out-of-sample
    obs_dom <- pop_domains_vec %in% unique(smp_domains_vec)
    dist_obs_dom <- unique(pop_domains_vec) %in% unique(smp_domains_vec)

  } else {
    # Population data available - therefore, aggregates (pop_cov and pop_mean)
    # will be produced

    fw_check_agg(
      pop_area_size = pop_area_size, pop_mean = pop_mean, pop_cov = pop_cov,
      mod_vars = mod_vars, smp_data = smp_data, fixed = fixed,
      smp_domains = smp_domains
    )

    smp_data <- smp_data[, smp_vars]
    pop_mean <- lapply(X   = pop_mean,
                       FUN = only_mod_vars,
                       var = mod_vars
    )
    pop_mean.mat <- matrix(data  = unlist(lapply(X   = pop_mean,
                                                 FUN = c_1)),
                           ncol  = length(mod_vars) + 1,
                           byrow = TRUE
    )
    row.names(pop_mean.mat) <- names(pop_mean)
    colnames(pop_mean.mat) <- c("intercept", mod_vars)

    if (!is.null(pop_cov)) {

      pop_cov <- lapply(X   = pop_cov,
                        FUN = only_mod_vars,
                        var = mod_vars
      )
      pop_cov.mat <- matrix(data  = unlist(lapply(X   = pop_cov,
                                                  FUN = crbind_0)),
                            ncol  = (length(mod_vars) + 1)^2,
                            byrow = TRUE
      )
      row.names(pop_cov.mat) <- names(pop_cov)
      colnames(pop_cov.mat) <- cov_names(c("intercept", mod_vars))
    } else {
      pop_cov.mat <- NULL
    }

    smp_data <- smp_data[order(smp_data[[smp_domains]]), ]
    smp_data[[smp_domains]] <-
      factor(smp_data[[smp_domains]], levels = names(pop_area_size))
    smp_domains_vec <- smp_data[[smp_domains]]
    smp_domains_vec <- droplevels(smp_domains_vec)

    # Number of households in population
    N_pop <- sum(pop_area_size)
    # Number of domains in the population
    N_dom_pop <- length(pop_area_size)
    # Number of households in population per domain
    n_pop <- as.vector(pop_area_size)

    # Indicator variables that indicate if domain is in- or out-of-sample
    dist_obs_dom <- unique(names(pop_area_size)) %in% unique(smp_domains_vec)

    pop_domains_vec <- NULL
    obs_dom <- NULL
  }

  # Number of households in sample
  N_smp <- length(smp_domains_vec)
  # Number of out-of-sample households
  N_unobs <- N_pop - N_smp
  # Number of domains in the sample
  N_dom_smp <- length(unique(smp_domains_vec))
  # Number of out-of-sample domains
  N_dom_unobs <- N_dom_pop - N_dom_smp
  # Number of households in sample per domain
  smp_domains_vec_tmp <- as.numeric(smp_domains_vec)
  n_smp <- as.vector(table(smp_domains_vec_tmp))

  return(list(pop_data        = pop_data,
              pop_domains_vec = pop_domains_vec,
              pop_area_size   = pop_area_size,
              pop_mean.mat    = pop_mean.mat,
              pop_cov.mat     = pop_cov.mat,
              smp_data        = smp_data,
              smp_domains_vec = smp_domains_vec,
              smp_domains     = smp_domains,
              N_pop           = N_pop,
              N_smp           = N_smp,
              N_unobs         = N_unobs,
              N_dom_pop       = N_dom_pop,
              N_dom_smp       = N_dom_smp,
              N_dom_unobs     = N_dom_unobs,
              n_pop           = n_pop,
              n_smp           = n_smp,
              obs_dom         = obs_dom,
              dist_obs_dom    = dist_obs_dom
  ))
}
