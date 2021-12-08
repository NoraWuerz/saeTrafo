# Functions called in framework_NER
fw_check_pop <- function(pop_data, mod_vars, pop_domains, smp_data,
                          fixed, smp_domains) {

  if (!all(mod_vars %in% colnames(pop_data))) {
    stop(paste0("Variable ", mod_vars[which(!(mod_vars %in% colnames(smp_data)))], " is not contained in pop_data.
                Please provide valid variable names for the explanatory variables."))
  }
  if (!(pop_domains %in% colnames(pop_data))) {
    stop(paste0("The domain variable ", pop_domains, " is not contained in pop_data.
                Please provide valid variable name for pop_domains."))
  }
  if (!all(mod_vars %in% colnames(smp_data))) {
    stop(paste0("Variable ", mod_vars[which(!(mod_vars %in% colnames(smp_data)))], " is not contained in smp_data.
                 Please provide valid variable names for the explanatory variables."))
  }
  if (!(smp_domains %in% colnames(smp_data))) {
    stop(paste0("The domain variable ", smp_domains, " is not contained in smp_data.
                 Please provide valid variable name for smp_domains."))
  }
  if (!((as.character(fixed[2])) %in% colnames(smp_data))) {
    stop(paste0("Variable ", as.character(fixed[2]), " is not contained in smp_data.
                Please provide valid variable name for the dependent variable."))
  }

  if (!is.numeric(smp_data[[paste(fixed[2])]])) {
    stop(paste0(as.character(fixed[2]), " must be the name of a variable that
               is a numeric vector."))
  }

  if (dim(pop_data)[1] < dim(smp_data)[1]) {
    stop("The population data set cannot have less observations than the
         sample data set.")
  }

  if (any(is.na(pop_data)) || any(is.na(smp_data))) {
    stop("NER_Trafo does not work with missing values.")
  }

  print("SAE methods for full population data and transformations are additionally provided in the R package emdi, which offers further functionalities")
}

fw_check_agg <- function(pop_area_size, pop_mean, pop_cov, mod_vars,
                         smp_data, fixed, smp_domains) {
  # check: Domain-names in all three objects
  # check: mod_vars includes in pop_cov and pop_mean
    # alle gleich benannt?

  if (!all(mod_vars %in% colnames(smp_data))) {
    stop(paste0("Variable ", mod_vars[which(!(mod_vars %in% colnames(smp_data)))], " is not contained in smp_data.
                 Please provide valid variable names for the explanatory variables."))
  }
  if (!(smp_domains %in% colnames(smp_data))) {
    stop(paste0("The domain variable ", smp_domains, " is not contained in smp_data.
                 Please provide valid variable name for smp_domains."))
  }
  if (!((as.character(fixed[2])) %in% colnames(smp_data))) {
    stop(paste0("Variable ", as.character(fixed[2]), " is not contained in smp_data.
                Please provide valid variable name for the dependent variable."))
  }

  if (!is.numeric(smp_data[[paste(fixed[2])]])) {
    stop(paste0(as.character(fixed[2]), " must be the name of a variable that
               is a numeric vector."))
  }

  if (any(is.na(smp_data))) {
    stop("NER_Trafo does not work with missing values.")
  }
}
