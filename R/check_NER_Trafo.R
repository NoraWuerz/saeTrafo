# Function called in NER_Trafo

NER_check1 <- function(fixed, pop_data, pop_domains, smp_data, smp_domains,
                       pop_mean, pop_area_size) {

  # check population data
  if (!is.null(pop_data)) {
    if (!is.character(pop_domains) || length(pop_domains) != 1) {
      stop(paste("Pop_domains must be a vector of lenght 1 and of class",
                 "character specifying the variable name of a numeric or",
                 "factor variable indicating domains in the population data.",
                 "See also help(NER_Trafo)."))
    }
    if (!is.data.frame(pop_data)) {
      stop("Population data must be a data frame. See also help(NER_Trafo)")
    }
  }

  # check sample data
  if (!is.data.frame(smp_data)) {
    stop(paste("Smp_data must be a data frame containing sample data. See also",
               "help(NER_Trafo)."))
  }
  if (!is.character(smp_domains) || length(smp_domains) != 1) {
    stop(paste("Smp_domains must be a vector of lenght 1 and of class",
               "character specifying the variable (name) of a numeric or",
               "factor variable indicating domains in the sample data. See",
               "also help(NER_Trafo)."))
  }

  # check domains
  if (!is.null(pop_data)) {
    if (!all(unique(as.character(smp_data[[smp_domains]])) %in%
             unique(as.character(pop_data[[pop_domains]])))) {
      stop(paste("The sample data contains domains that are not contained in",
                 "the population data."))
    }
  }

}

NER_check2 <- function(threshold, transformation, interval, MSE, B, cpus,
                       seed) {
  # check threshold
  if (!is.null(threshold) &&
      !(is.numeric(threshold) && length(threshold) == 1)) {
    stop(paste("threshold needs to be a single numeric value. If it is NULL it",
               "is set to 30. See also help(NER_Trafo)."))
  }
  if (is.null(transformation) ||
      !(transformation == "log" || transformation == "log.shift" ||
        transformation == "no")) {
    stop(paste("The three options for transformation are ''no'', ''log'' or",
               "''log.shift''."))
  }

  # check interval
  if (any(interval != "default") & (!is.vector(interval, mode = "numeric") ||
                                    length(interval) != 2 ||
                                    !(interval[1] < interval[2]))) {
    stop(paste("Interval needs to be a numeric vector of length 2 defining a",
               "lower and upper limit for the estimation of the optimal",
               "transformation parameter. The value of the lower limit needs",
               "to be smaller than the upper limit. You can also choose",
               "'default'. See also help(NER_Trafo)."))
  }

  # check MSE
  if (!is.logical(MSE) || length(MSE) != 1) {
    stop(paste("MSE must be a logical value. Set MSE to TRUE or FALSE. See",
               "also help(NER_Trafo)."))
  }
  if (MSE == TRUE && !(is.numeric(B) && length(B) == 1  && B > 1)) {
    stop(paste("If MSE is set to TRUE, a single numeric value for the number",
               "of bootstrap sample needs to be chosen that is greater than 1.",
               "See also help(NER_Trafo)."))
  }
  if (!is.numeric(cpus) || !(is.numeric(cpus) && length(cpus) == 1)) {
    stop(paste("Cpus must be a single number determining the number of kernels",
               "for the parallelization."))
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || !(is.numeric(seed) && length(seed) == 1))) {
    stop(paste("The seed must be a single value, interpreted as an integer,",
               "or NULL See also help(NER_Trafo)."))
  }
}

# Function called in framework_NER.R and check inputs in the case of available
# pop_data
fw_check_pop <- function(pop_data, mod_vars, pop_domains, smp_data,
                         fixed, smp_domains) {

  if (!all(mod_vars %in% colnames(pop_data))) {
    stop(paste0("Variable ",
                mod_vars[which(!(mod_vars %in% colnames(smp_data)))],
                " is not contained in pop_data. Please provide valid variable",
                " names for the explanatory variables."))
  }
  if (!(pop_domains %in% colnames(pop_data))) {
    stop(paste0("The domain variable ",
                pop_domains,
                " is not contained in pop_data. Please provide valid variable",
                " name for pop_domains."))
  }
  if (!all(mod_vars %in% colnames(smp_data))) {
    stop(paste0("Variable ",
                mod_vars[which(!(mod_vars %in% colnames(smp_data)))],
                " is not contained in smp_data. Please provide valid variable",
                " names for the explanatory variables."))
  }
  if (!(smp_domains %in% colnames(smp_data))) {
    stop(paste0("The domain variable ",
                 smp_domains,
                 " is not contained in smp_data. Please provide valid variable",
                 " name for smp_domains."))
  }
  if (is.null(fixed)  || !inherits(fixed, "formula")) {
    stop("Fixed must be a formula object. See also help(NER_Trafo).")
  }
  if (!((as.character(fixed[2])) %in% colnames(smp_data))) {
    stop(paste0("Variable ",
                as.character(fixed[2]),
                " is not contained in smp_data. Please provide valid variable",
                " name for the dependent variable."))
  }
  if (!is.numeric(smp_data[[paste(fixed[2])]])) {
    stop(paste0(as.character(fixed[2]),
                " must be the name of a variable that is a numeric vector."))
  }
  if (dim(pop_data)[1] < dim(smp_data)[1]) {
    stop(paste("The population data set cannot have less observations than the",
               "sample data set."))
  }
  if (any(is.na(pop_data)) || any(is.na(smp_data))) {
    stop("NER_Trafo does not work with missing values.")
  }
  print(paste("More SAE methods for full population data and transformations",
              "are offered in the R package emdi."))
}

# Function called in framework_NER.R and check inputs in the case of
# no pop_data
fw_check_agg <- function(pop_area_size, pop_mean, pop_cov, mod_vars,
                         smp_data, fixed, smp_domains) {

  if (sum(pop_area_size) == nrow(smp_data)) {
    stop(paste("The sample size is equal to the population size. Please check",
               "pop_area_size and smp_data."))
  }

  if (is.null(pop_mean) || is.null(pop_area_size)) {
    stop(paste("Population mean and population area sizes (and optional",
               "population covariances) must be available or population data",
               "must be provided. See also help(NER_Trafo)."))
  }
  if (!all(names(pop_area_size) %in% names(pop_mean))) {
    stop("Domains in pop_area_size not equal to domains in pop_mean.")
  }
  if (!all(unlist(lapply(pop_mean,
                         mod_vars = mod_vars,
                         FUN = all_in_domain)))) {
    stop(paste0("Population mean vector for domain ",
                names(pop_mean)[!unlist(lapply(pop_mean,
                                               mod_vars = mod_vars,
                                               FUN = all_in_domain))],
                " has variables that are not contained in smp_data. Please",
                " provide valid variable names for the explanatory variables."))
  }
  if (!is.null(pop_cov)) {
    if (!all(names(pop_area_size) %in% names(pop_cov))) {
      stop("Domains in pop_area_size not equal to domains in pop_cov.")
    }
    if (!all(unlist(lapply(pop_cov,
                           mod_vars = mod_vars,
                           FUN = all_in_domain)))) {
      stop(paste0("Population covariance matrix for domain ",
                  names(pop_mean)[!unlist(lapply(pop_cov,
                                                 mod_vars = mod_vars,
                                                 FUN = all_in_domain))],
                  " has variables that are not contained in smp_data. Please",
                  " provide valid variable names for the explanatory",
                  " variables."))
    }
  }
  if (!all(mod_vars %in% colnames(smp_data))) {
    stop(paste0("Variable ",
                mod_vars[which(!(mod_vars %in% colnames(smp_data)))],
                " is not contained in smp_data. Please provide valid variable",
                " names for the explanatory variables."))
  }
  if (!(smp_domains %in% colnames(smp_data))) {
    stop(paste0("The domain variable ",
                smp_domains,
                " is not contained in smp_data. Please provide valid variable",
                " name for smp_domains."))
  }
  if (!all(smp_data[smp_domains][, 1] %in% names(pop_area_size))) {
    stop("More than one sample domain is not inclueded in the population data.")
  }
  if (is.null(fixed)  || !inherits(fixed, "formula")) {
    stop("Fixed must be a formula object. See also help(NER_Trafo).")
  }
  if (!((as.character(fixed[2])) %in% colnames(smp_data))) {
    stop(paste0("Variable ",
                as.character(fixed[2]),
                " is not contained in smp_data. Please provide valid variable",
                " name for the dependent variable."))
  }
  if (!is.numeric(smp_data[[paste(fixed[2])]])) {
    stop(paste0(as.character(fixed[2]),
                " must be the name of a variable that is a numeric vector."))
  }
  if (any(is.na(smp_data))) {
    stop("NER_Trafo does not work with missing values.")
  }
}
