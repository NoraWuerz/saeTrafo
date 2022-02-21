#' @rdname write.excel
#' @importFrom readODS write_ods
#' @export

write.ods <- function(object,
                      file      = "ods_output.ods",
                      MSE       = FALSE,
                      CV        = FALSE,
                      split     = FALSE) {

  writeods_check(object = object, file = file, split = split)

  wb <- gsub(".ods", "", file)

  if (inherits(object, "NER"))  {
    add_summary_ods_NER(object = object, wb = wb)
  }

  if (!split & (MSE | CV)) {
    add_estims_ods(object       = object,
                   wb           = wb,
                   MSE          = MSE,
                   CV           = CV
    )
  } else {
    add_pointests_ods(wb           = wb,
                      object       = object
    )
    if (MSE || CV) {
      add_precisions_ods(object = object,
                         MSE = MSE,
                         wb = wb,
                         CV = CV
      )
    }
  }

}

add_summary_ods_NER <- function(object, wb, headlines_cs) {

  su <- summary(object)

  df_nobs <- data.frame(Count = c(su$out_of_smp,
                                  su$in_smp, su$size_pop,
                                  su$size_smp)
  )
  rownames(df_nobs) <- c("out of sample domains",
                         "in sample domains",
                         "out of sample observations",
                         "in sample observations")

  df_nobs <- cbind.data.frame(rownames(df_nobs), df_nobs)
  write_ods(x = df_nobs, path = paste0(wb, "_sumObs", ".ods"))

  df_size_dom <- as.data.frame(su$size_dom)
  df_size_dom <- cbind.data.frame(rownames(df_size_dom), df_size_dom)
  write_ods(x = df_size_dom, path = paste0(wb, "_sumDom", ".ods"))

  if (!is.null(su$transform)) {
    write_ods(x = su$transform, path = paste0(wb, "_sumTrafo", ".ods"))
  }

  su$normality <-  cbind.data.frame(rownames(su$normality), su$normality)
  write_ods(x = su$normality, path = paste0(wb, "_sumNorm", ".ods"))

  su$coeff_determ <-  cbind.data.frame("Coefficients of determination",
                                       su$coeff_determ)
  write_ods(x = su$coeff_determ, path = paste0(wb, "_sumCoefDet", ".ods"))

  return(NULL)
}

add_pointests_ods <- function(object, wb, headlines_cs) {

  data <- point_saeTrafo(object = object)$ind
  data[, 1] <- iconv(x = data[, 1], from = "", to = "UTF-8")
  write_ods(x = data, path = paste0(wb, "_pointEstim", ".ods"))

  return(NULL)
}

add_precisions_ods <- function(object, MSE, wb, headlines_cs, CV) {

  precisions <- mse_saeTrafo(object = object, CV = TRUE)
  if (MSE) {
    precisions$ind[, 1] <-
      iconv(x <- precisions$ind[, 1], from = "", to = "UTF-8")
    write_ods(x = precisions$ind, path = paste0(wb, "_precMSE", ".ods"))
  }
  if (CV) {
    precisions$ind_cv[, 1] <-
      iconv(x <- precisions$ind_cv[, 1], from = "", to = "UTF-8")
    write_ods(x = precisions$ind_cv, path = paste0(wb, "_precCV", ".ods"))
  }
  return(NULL)
}

#' @importFrom readODS write_ods

add_estims_ods <- function(object, wb, headlines_cs, MSE, CV) {

  data <- estimators(object = object, MSE = MSE, CV = CV)$ind
  data[, 1] <- iconv(x = data[, 1], from = "", to = "UTF-8")

  write_ods(x = data, path = paste0(wb, "_estim", ".ods"))
  return(NULL)
}
