#' Nested error regression Model (Battese) with transformations
#'


NER_Trafo <- function(fixed,
                      pop_area_size = NULL,
                      pop_mean = NULL,
                      pop_cov = NULL,
                      pop_data = NULL,
                      pop_domains,
                      smp_data,
                      smp_domains,
                      threshold = 30,
                      B = 50,
                      transformation = "log.shift",
                      interval = 'default',
                      MSE = FALSE,
                      seed = 123) {
  # NER_check1()
  # NER_check2()

  # Save function call ---------------------------------------------------------

  call <- match.call()
  if (inherits(call$fixed, "name")) {
    call$fixed <- fixed
  }



  # Data manipulation and notational framework ---------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # The function framework_NER can be found in script framework_NER.R
  framework <- framework_NER(
    pop_area_size = pop_area_size,
    pop_mean = pop_mean,
    pop_cov = pop_cov,
    pop_data = pop_data,
    pop_domains = pop_domains,
    smp_data = smp_data,
    smp_domains = smp_domains,
    fixed = fixed
  )

  NER_out <- framework


  # Point Estimation -----------------------------------------------------------
  # The function point_estim can be found in script point_estimation.R
  point_estim <- point_estim(
     framework = framework,
     fixed = fixed,
     transformation = transformation,
     threshold = threshold,
     interval = interval,
     keep_data = TRUE
  )



  # MSE Estimation -------------------------------------------------------------

  # if (MSE == TRUE) {
  #
  #   # The function parametric_bootstrap can be found in script mse_estimation.R
  #   mse_estimates <- parametric_bootstrap(
  #     framework = framework,
  #     point_estim = point_estim,
  #     fixed = fixed,
  #     transformation = transformation,
  #     interval = interval,
  #     B = B,
  #   )
  #
  #
  #
  #   ebp_out <- list(
  #     ind = point_estim$ind,
  #     MSE = mse_estimates,
  #     transform_param = point_estim[c("optimal_lambda", "shift_par")],
  #     model = point_estim$model,
  #     framework = framework[c(
  #       "N_dom_unobs",
  #       "N_dom_smp",
  #       "N_smp",
  #       "N_pop",
  #       "smp_domains",
  #       "smp_data",
  #       "smp_domains_vec",
  #       "pop_domains_vec"
  #     )],
  #     transformation = transformation,
  #     method = "reml",
  #     fixed = fixed,
  #     call = call,
  #     successful_bootstraps = NULL
  #   )
  # } else {
  #   ebp_out <- list(
  #     ind = point_estim$ind,
  #     MSE = NULL,
  #     transform_param = point_estim[c("optimal_lambda", "shift_par")],
  #     model = point_estim$model,
  #     framework = framework[c(
  #       "N_dom_unobs",
  #       "N_dom_smp",
  #       "N_smp",
  #       "N_pop",
  #       "smp_domains",
  #       "smp_data",
  #       "smp_domains_vec",
  #       "pop_domains_vec",
  #       "response"
  #     )],
  #     transformation = transformation,
  #     method = "reml",
  #     fixed = fixed,
  #     call = call,
  #     successful_bootstraps = NULL
  #   )
  # }

  class(NER_out) <- c("NER", "SAE_Trafo")
  return(NER_out)
}
