compare_plot_check <- function(model, indicator, label, color, shape,
                                line_type, gg_theme) {

  if (is.null(label) || (!(label == "orig" || label == "no_title" ||
                           label == "blank"))) {
    stop(paste("label can be one of the following characters 'orig',",
               "'no_title' or 'blank'."))
  }
  if (length(color) != 2 || !is.vector(color)) {
    stop(paste("color needs to be a vector of length 2 defining the two colors",
               "for the scatter and line plot. See also help(compare_plot)."))
  }
  if (length(shape) != 2 || !is.vector(shape)) {
    stop(paste("shape needs to be a vector of length 2 defining the two shapes",
               "for the estimates in the line plots. See also",
               "help(compare_plot)."))
  }
  if (length(line_type) != 2 || !is.vector(shape)) {
    stop(paste("line_type needs to be a vector of length 2 defining the types",
               "for the lines in the line plots. See also help(compare_plot)."))
  }
  if (!all(line_type %in% c("twodash", "solid", "longdash", "dotted", "dotdash",
                            "dashed", "blank"))) {
    stop(paste("An element in argument line_type is not a valid option. See",
               "help(compare_plot) for valid options."))
  }
}

compare_plot_check2 <- function(ind_direct, ind_model) {

  if (!any(ind_direct$Domain %in% ind_model$Domain)) {
    stop(paste("Domain identifiers between direct and model estimates never",
               "match. Please verify you are comparing estimates obtained on",
               "the same sample."))
  }
  if (!all(ind_direct$Domain %in% ind_model$Domain)) {
    warning(paste("Not all domains contained in the direct estimation have",
                  "been found in the model estimation. Following plots will",
                  "only contain results for estimates available in both",
                  "objects"))
  }
  if (!all(ind_model$Domain %in% ind_direct$Domain)) {
    message(paste("Not all domains contained in the model estimation have been",
                  "found in the direct estimation. Following plots will only",
                  "contain results for estimates available in both objects."))
  }
}


compare_pred_check <- function(object1, object2, MSE) {

  if (!inherits(object1, "saeTrafo") || !(inherits(object2, "emdi") |
                                          inherits(object2, "saeTrafo"))) {
    stop(paste("Object 1 must be of class saeTrafo object 2 of class saeTrafo",
               "or emdi."))
  }
  if ((length(object1$ind$Domain) == length(object2$ind$Domain)) &&
      (!all(as.character(object1$ind$Domain) %in%
            as.character(object2$ind$Domain)))) {
    stop("It is only possible to compare objects with the same domains.")
  }
  if ((length(object1$ind$Domain) < length(object2$ind$Domain)) &&
      !all(as.character(object1$ind$Domain) %in%
           as.character(object2$ind$Domain))) {
    stop(paste("The first object contains domains that are not contained in",
               "the second object. It is only possible to compare objects with",
               "the same domains."))
  }
  if ((length(object2$ind$Domain) < length(object1$ind$Domain)) &&
      !all(as.character(object2$ind$Domain) %in%
           as.character(object1$ind$Domain))) {
    stop(paste("The second object contains domains that are not contained in",
               "the first object. It is only possible to compare objects with",
               "the same domains."))
  }
  if ((MSE == TRUE) && (is.null(object1$MSE) || is.null(object2$MSE))) {
    stop("If MSE is set to TRUE, both objects need to contain MSE estimates.")
  }
}
