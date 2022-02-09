#' Compare predictions of model objects
#'
#' Function \code{compare_pred} is a generic function used to compare
#' predictions of two model objects.
#'
#' @param object1 an object of type "saeTrafo".
#' @param object2 an object of type "saeTrafo" or "emdi" (\code{\link{emdi}}).
#' @param MSE if \code{TRUE}, MSE estimates are returned. Defaults to
#' \code{FALSE} and than point estimates are returned.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @name compare_pred

compare_pred <- function(object1, object2, MSE = FALSE, ...)
  UseMethod("compare_pred")


#' Compare predictions of saeTrafo objects
#'
#' Method \code{compare_pred.saeTrafo} compares predictions of two saeTrafo
#' objects or a saeTrafo object and a emdi object.
#'
#' @param object1 an object of type "saeTrafo".
#' @param object2 an object of type "saeTrafo" or "emdi" (\code{\link{emdi}}).
#' @param MSE if \code{TRUE}, MSE estimates are returned. Defaults to
#' \code{FALSE} and than point estimates are returned.
#' @param ... further arguments passed to or from other methods.
#' @return Data frame containing the point estimates or the MSE estimates
#' (if \code{MSE} is set to \code{TRUE}\) of both objects. If column names are
#' duplicated, the suffixes "_1" and "_2" are added to their names. "_1" and
#' "_2" standing for object1 and object2, respectively.
#' @seealso \code{\link{emdi}}, \code{\link{NER_Trafo}},
#' \code{\link{saeTrafo.object}}
#' @examples
#' \donttest{
#' # Example comparing two saeTrafo objects
#' # (the second object can be of the class emdi)
#' formel <- eqIncome ~ gender + eqsize + cash + self_empl + unempl_ben +
#' age_ben + surv_ben + sick_ben + dis_ben + rent + fam_allow + house_allow +
#' cap_inv + tax_adj
#'
#' res_ls <- NER_Trafo(fixed = formel, smp_domains = "district",
#' pop_area_size = pop_area_size, pop_mean = pop_mean, pop_cov = pop_cov,
#' smp_data = eusilcA_smp, MSE = F, B = 3, cpus = 4, seed = 3)
#' }
#' @export
#' @rdname compare_pred
#' @method compare_pred saeTrafo


compare_pred.saeTrafo <- function(object1, object2, MSE = FALSE, ...) {

  if (!inherits(object1, "saeTrafo") || !(inherits(object2, "emdi") |
                                          inherits(object2, "saeTrafo"))) {
    stop("Object 1 must be of class saeTrafo object 2 of class saeTrafo or
         emdi.")
  }

  if ((length(object1$ind$Domain) == length(object2$ind$Domain)) &&
      (!all(as.character(object1$ind$Domain) %in%
            as.character(object2$ind$Domain)))) {
    stop("It is only possible to compare objects with the same domains.")
  }

  if ((length(object1$ind$Domain) < length(object2$ind$Domain)) &&
      !all(as.character(object1$ind$Domain) %in%
           as.character(object2$ind$Domain))) {
    stop("The first object contains domains that are not contained in the second
         object. It is only possible to compare emdi objects with the same
         domains.")
  }

  if ((length(object2$ind$Domain) < length(object1$ind$Domain)) &&
      !all(as.character(object2$ind$Domain) %in%
           as.character(object1$ind$Domain))) {
    stop("The second object contains domains that are not contained in the first
         object. It is only possible to compare objects with the same domains.")
  }

  if ((MSE == TRUE) && (is.null(object1$MSE) || is.null(object2$MSE))) {
    stop("If MSE is set to TRUE, both objects need to contain MSE estimates.")
  }


  if (MSE == FALSE) {
    object1data <- get("ind", object1)
    object2data <- get("ind", object2)
  } else if (MSE == TRUE) {
    object1data <- get("MSE", object1)
    object2data <- get("MSE", object2)
  }

  if (inherits(object1, "fh")) {
    object1data <- object1data[, -4] # remove column Out
  }

  if (inherits(object2, "fh")) {
    object2data <- object2data[, -4] # remove column Out
  }

  order_direct_ebp <- c("Domain", "Mean_1", "Mean_2", "Head_Count_1",
                        "Head_Count_2", "Poverty_Gap_1", "Poverty_Gap_2",
                        "Gini_1", "Gini_2", "Quintile_Share_1",
                        "Quintile_Share_2", "Quantile_10_1", "Quantile_10_2",
                        "Quantile_25_1", "Quantile_25_2", "Median_1",
                        "Median_2", "Quantile_75_1", "Quantile_75_2",
                        "Quantile_90_1", "Quantile_90_2")

  if ((inherits(object1, "ebp") && inherits(object2, "ebp")) ||
      (inherits(object1, "direct") && inherits(object2, "direct")) ||
      (inherits(object1, "direct") && inherits(object2, "ebp")) ||
      (inherits(object1, "ebp") && inherits(object2, "direct")) ||
      (inherits(object1, "NER") && inherits(object2, "NER")) ||
      (inherits(object1, "direct") && inherits(object2, "NER")) ||
      (inherits(object1, "NER") && inherits(object2, "direct")) ||
      (inherits(object1, "ebp") && inherits(object2, "NER")) ||
      (inherits(object1, "NER") && inherits(object2, "ebp"))
      ) {

    if (dim(object1data)[2] > 11) {
      colnames(object1data)[12:dim(object1data)[2]] <-
        paste0(names(object1data[12:dim(object1data)[2]]), "_1")
    }

    if (dim(object2data)[2] > 11) {
      colnames(object2data)[12:dim(object2data)[2]] <-
        paste0(names(object2data[12:dim(object2data)[2]]), "_2")
    }

    data <- merge(object1data, object2data, by.x = "Domain", by.y = "Domain",
                  suffixes = c("_1", "_2"))

    if (dim(data)[2] == 21) {
      data <- data[, order_direct_ebp]
    } else if (dim(data)[2] > 21) {
      custom_indicators <- colnames(data)[
        (which(!colnames(data) %in% order_direct_ebp))
        ]
      data <- data[, c(order_direct_ebp, custom_indicators)]
    }

  } else if (inherits(object1, "fh") && inherits(object2, "fh")) {
    data <- merge(object1data, object2data, by.x = "Domain", by.y = "Domain",
                  suffixes = c("_1", "_2"))
    data <- data[, c("Domain", "Direct_1", "Direct_2", "FH_1", "FH_2")]
  } else if ((inherits(object1, "direct") && inherits(object2, "fh")) ||
             (inherits(object1, "fh") && inherits(object2, "direct"))) {
    data <- merge(object1data, object2data, by.x = "Domain", by.y = "Domain",
                  all = TRUE, suffixes = c("_1", "_2"))
  } else if ((inherits(object1, "fh") && inherits(object2, "ebp")) ||
             (inherits(object1, "ebp") && inherits(object2, "fh")) ||
             (inherits(object1, "fh") && inherits(object2, "NER")) ||
             (inherits(object1, "NER") && inherits(object2, "fh"))) {
    data <- merge(object1data, object2data, by.x = "Domain", by.y = "Domain",
                  all = TRUE, suffixes = c("_1", "_2"))
  }
  data
}
