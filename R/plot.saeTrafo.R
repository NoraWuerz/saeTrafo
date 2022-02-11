#' Plots for an saeTrafo object
#'
#' Diagnostic plots of the nested error regression model
#' (see also \code{\link{NER_Trafo}}) are obtained. These include Q-Q plots and
#' density plots of residuals and random effects, a Cook's distance plot for
#' detecting outliers and the. The return depends on the
#' transformation such that a plot for the optimal parameter is
#' only returned in case a transformation with transformation parameter is chosen.
#' The range of the x-axis is optional but necessary to change if there are convergence problems.
#' All plots are obtained by \code{\link[ggplot2]{ggplot}}.
#' @param x an object of type "emdi", either "ebp" or "fh", representing point
#' and, if chosen, MSE estimates obtained by the EBP or Fay-Herriot approach
#' (see also \code{\link{ebp}} and \code{\link{fh}}).
#' @param label argument that enables to customize title and axis labels. There
#' are three instant options to label the diagnostic plot: (i) original labels ("orig"),
#' (ii) axis lables but no title ("no_title"), (iii) neither axis
#' labels nor title ("blank").
#' (iv) individual labels by a list that needs to
#' have below structure. Six elements can be defined called \code{qq_res, qq_ran,
#' d_res, d_ran, cooks} and \code{opt_lambda} for the six different plots and these
#' list elements need to have three elements each called \code{title, y_lab and
#' x_lab}. Only the labels for the plots that should be different to the original
#' need to be specified. Please see the details section for an example
#' with the default labels.
#' @param color a character vector with two elements. The first element defines
#' the color for the line in the QQ-plots, for the Cook's Distance plot and for
#' the Box-Cox plot. The second element defines the color for the densities.
#' @param gg_theme \code{\link[ggplot2]{theme}} list from package \pkg{ggplot2}.
#' For using this argument, package \pkg{ggplot2} must be loaded via
#' \code{library(ggplot2)}. See also Example 4.
#' @param cooks if \code{TRUE}, a Cook's distance plot is returned when the ebp
#' function is used. The used method \code{mdffits.default} from the package \pkg{HLMdiag}
#' struggles when data sets get large. In these cases, \code{cooks} should be
#' set to \code{FALSE}. It defaults to \code{TRUE}.
#' @param range optional sequence determining the range of the x-axis for plots
#' of the optimal transformation parameter that defaults to \code{NULL}. In that
#' case a range of the default interval is used for the plots of the
#' optimal parameter. This leads in some cases to convergence problems such that
#' it should be changed to e.g. the selected \code{interval}. The default value
#' depends on the chosen data driven transformation and equals the default
#' interval for the estimation of the optimal parameter.
#' @param ... optional arguments passed to generic function.
#' @return Two Q-Q plots in one grid, two density plots, a Cook's distance plot
#' and a likelihood plot for the optimal parameter of transformations with transformation
#' parameter obtained by \code{\link[ggplot2]{ggplot}}. The latter two plots are only provided
#' for ebp object.
#' @details The default settings of the \code{label} argument are as follows (please
#' note that the title for opt_lambda depends on the chosen transformation, for
#' the example Box-Cox is shown):\cr
##' \describe{
##' \item{list(}{}
##' \item{qq_res =}{c(title="Error term", y_lab="Quantiles of pearson residuals",
##'                 x_lab="Theoretical quantiles"),}
##' \item{qq_ran =}{c(title="Random effect",
##'                 y_lab="Quantiles of random effects",
##'                 x_lab="Theoretical quantiles"),}
##' \item{d_res =}{c(title="Density - Pearson residuals",
##'                y_lab="Density",
##'                x_lab="Pearson residuals"),}
##' \item{d_ran =}{c(title="Density - Standardized random effects",
##'                y_lab="Density",
##'                x_lab="Standardized random effects"),}
##' \item{cooks =}{c(title="Cook's Distance Plot",
##'                y_lab="Cook's Distance",
##'                x_lab="Index"),}
##' \item{opt_lambda =}{c(title="Box-Cox - REML",
##'                y_lab="Log-Likelihood",
##'                x_lab="expression(lambda)"))}
##' }
#' @seealso \code{\link{saeTrafoObject}}, \code{\link{NER_Trafo}}
#' @importFrom ggplot2 qplot geom_abline ggtitle ylab xlab ggplot stat_qq
#' @importFrom ggplot2 aes geom_point geom_smooth coord_fixed geom_line
#' @importFrom ggplot2 scale_color_manual scale_fill_manual geom_segment
#' @importFrom ggplot2 scale_linetype_discrete geom_density geom_text
#' @importFrom ggplot2 geom_line geom_vline stat_function geom_qq
#' @importFrom nlme ranef random.effects
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom stats shapiro.test logLik cooks.distance
#' @importFrom HLMdiag mdffits
#' @importFrom stringr str_to_title
#' @name plot.saeTrafo
#' @export

plot.saeTrafo <- function(x,
                      label = "orig",
                      color = c("blue", "lightblue3"),
                      gg_theme = NULL,
                      cooks = TRUE,
                      range = NULL, ...) {

  plot_check(x = x, label = label, color = color, cooks = cooks, range = range)

  # Preparation for plots ----
  residuals <- residuals(x$model, level = 0, type = "pearson")
  rand.eff <- nlme::ranef(x$model)$'(Intercept)'
  srand.eff <- (rand.eff - mean(rand.eff)) / sd(rand.eff)
  tmp <- as.matrix(random.effects(x$model))[,1]
  model <- x$model
  model$call$fixed <- x$fixed
  cook_df <- NULL
  indexer <- NULL
  likelihoods <- NULL
  opt_lambda <- FALSE

  if (cooks == TRUE) {
    cooksdist <- NULL
    # Supress warning is used here due to a small bug in the
    # HLMdiag:::.extract.lmeDesign which is underlying the here used
    # cooks distance. The given warning has no relevance in this case.
    try(cooksdist <- as.vector(suppressWarnings(cooks.distance(model))),
        silent = TRUE)
    if (is.null(cooksdist)) {
      cooks <- FALSE
      warning("Cook's distance could not be calculated, this is usually due to
               exceedence of available memory. Try using cooks = FALSE to avoid
               this message and improve computation time.")
    } else {
      cook_df <- data.frame(index = seq_along(cooksdist), cooksdist)
      indexer <- cook_df[order(cooksdist, decreasing = TRUE),][seq_len(3),]
    }
  }

  if (x$transformation == "log.shift") {
    opt_lambda = TRUE

    if (is.null(range)) {

      range <- seq(x$transform_param$optimal_lambda - (0.5 * x$transform_param$optimal_lambda),
                   x$transform_param$optimal_lambda + (0.5 * x$transform_param$optimal_lambda),
                   length = 50)

    } else {
      range <- range
    }


    likelihoods <- vapply(range,
                          function(lam, fixed , smp_data, smp_domains,
                                   transformation)
                          {
                            result <- NULL
                            try(result <- -as.numeric(
                              generic_opt(lam, fixed, smp_data,
                                          smp_domains, transformation)),
                              silent = TRUE)
                            if (is.null(result)) result <- NA
                            result
                          }, numeric(1),   fixed = x$fixed,
                          smp_data = x$framework$smp_data,
                          smp_domains = x$framework$smp_domains,
                          transformation = x$transformation)

    if (any(is.na(likelihoods))) {
      warning("For some lambda in the chosen range, the likelihood does not
               converge. For these lambdas no likelihood is plotted. Choose a
               different range to avoid this behaviour")
    }
  }

  Residuals <- Random <- index <- lambda <- log_likelihood <- cooksdist <-  NULL

  plotList <- vector(mode = "list", length = 5)
  plotList <- lapply(plotList, function(x) NA)
  names(plotList) <- c("qq_plots", "density_res","density_ran",
                       "cooks_distance", "likelihood")

  label <- define_label(x = x, label = label)

  ## QQ Plots
  # Residuals
  res <- qplot(sample = residuals) +
    geom_abline(colour = color[1]) +
    ggtitle(label$qq_res["title"]) + ylab(label$qq_res["y_lab"]) +
    xlab(label$qq_res["x_lab"]) + gg_theme

  # Random effects
  ran <- ggplot(data.frame(tmp), aes(sample = tmp)) +
    stat_qq(distribution = qnorm, dparams = list(mean = mean(tmp),
                                                 sd = sd(tmp))) +
    geom_abline(intercept = 0, slope = 1, na.rm = TRUE, col = color[1]) +
    ggtitle(label$qq_ran["title"]) + ylab(label$qq_ran["y_lab"]) +
    xlab(label$qq_ran["x_lab"]) + gg_theme

  plotList[[1]] <- arrangeGrob(res, ran ,ncol = 2)
  grid.arrange(plotList[[1]])
  cat("Press [enter] to continue")
  line <- readline()

  print( (plotList[[2]] <- ggplot(data.frame(Residuals = residuals),
                                  aes(x = Residuals),
                                  fill = color[2], color = color[2]) +
            geom_density(fill = color[2], color = color[2],
                         alpha = 0.4) +
            stat_function(fun = dnorm) + ylab(label$d_res["y_lab"]) +
            xlab(label$d_res["x_lab"]) +
            ggtitle(label$d_res["title"]) + gg_theme))
  cat("Press [enter] to continue")
  line <- readline()
  print( (plotList[[3]] <- ggplot(data.frame(Random = srand.eff), aes(x = Random),
                                  fill = color[2], color = color[2]) +
            geom_density(fill = color[2], color = color[2],
                         alpha = 0.4) +
            stat_function(fun = dnorm) + ylab(label$d_ran["y_lab"]) +
            xlab(label$d_ran["x_lab"]) +
            ggtitle(label$d_ran["title"]) +
            gg_theme))

  if (cooks == TRUE) {
    cat("Press [enter] to continue")
    line <- readline()
    print((plotList[[4]] <- ggplot(data = cook_df, aes(x = index, y = cooksdist)) +
             geom_segment(aes(x = index, y = 0, xend = index, yend = cooksdist),
                          colour = color[1]) +
             xlab("Index") + ylab(label$cooks["y_lab"])
           + geom_text(label = indexer[,1], data = indexer) +
             ggtitle(label$cooks["title"]) + gg_theme))
  }

  if (opt_lambda == TRUE) {
    cat("Press [enter] to continue")
    line <- readline()

    if (any(label$opt_lambda["x_lab"] == "expression(lambda)") ||
        any(label$opt_lambda["x_lab"] == "expression(Lambda)")) {

      x_lab <- expression(lambda)
    } else {
      x_lab <- label$opt_lambda["x_lab"]
    }
    if (any(is.na(likelihoods))) {
      warning("For some lambda in the chosen range, the likelihood does not
               converge. For these lambdas no likelihood is plotted. Choose a
              different range to avoid this behaviour")
    }
    print((plotList[[5]] <- ggplot(data.frame(lambda = range,
                                              log_likelihood = likelihoods),
                                   aes(x = lambda, y = log_likelihood)) + geom_line() +
             xlab(x_lab) + ylab(label$opt_lambda["y_lab"]) +
             #geom_vline(xintercept = range[which.max(likelihoods)],
             #            colour = color[1]) + ggtitle(label$opt_lambda["title"]) +
             geom_vline(xintercept = x$transform_param$optimal_lambda,
                        colour = color[1]) + ggtitle(label$opt_lambda["title"]) +
             gg_theme))
  }
  invisible(plotList)
}



# Definition of the labels

define_label <- function(x, label){

   if (!inherits(label, "list")) {
    if (label == "orig") {

      if (inherits(x, "NER")) {
        label <- list(qq_res = c(title = "Error term",
                                 y_lab = "Quantiles of pearson residuals",
                                 x_lab = "Theoretical quantiles"),
                      qq_ran = c(title = "Random effect",
                                 y_lab = "Quantiles of random effects",
                                 x_lab = "Theoretical quantiles"),
                      d_res = c(title = "Density - Pearson residuals",
                                y_lab = "Density",
                                x_lab = "Pearson residuals"),
                      d_ran = c(title = "Density - Standardized random effects",
                                y_lab = "Density",
                                x_lab = "Standardized random effects"),
                      cooks = c(title = "Cook's Distance Plot",
                                y_lab = "Cook's Distance",
                                x_lab = "Index"),
                      opt_lambda = c(title = paste0(str_to_title(gsub("\\.","-", x$transformation)), ' - REML'),
                                  y_lab = "Log-Likelihood",
                                  x_lab = "expression(lambda)"))}


    } else if (label == "blank") {
      label <- list(qq_res = c(title = "",
                               y_lab = "",
                               x_lab = ""),
                    qq_ran = c(title = "",
                               y_lab = "",
                               x_lab = ""),
                    d_res = c(title = "",
                              y_lab = "",
                              x_lab = ""),
                    d_ran = c(title = "",
                              y_lab = "",
                              x_lab = ""),
                    cooks = c(title = "",
                              y_lab = "",
                              x_lab = ""),
                    opt_lambda = c(title = "",
                                y_lab = "",
                                x_lab = ""))
    } else if (label == "no_title") {

      if (inherits(x, "NER")) {
        label <- list(qq_res = c(title = "",
                                 y_lab = "Quantiles of pearson residuals",
                                 x_lab = "Theoretical quantiles"),
                      qq_ran = c(title = "",
                                 y_lab = "Quantiles of random effects",
                                 x_lab = "Theoretical quantiles"),
                      d_res = c(title = "",
                                y_lab = "Density",
                                x_lab = "Pearson residuals"),
                      d_ran = c(title = "",
                                y_lab = "Density",
                                x_lab = "Standardized random effects"),
                      cooks = c(title = "",
                                y_lab = "Cook's Distance",
                                x_lab = "Index"),
                      opt_lambda = c(title = "",
                                  y_lab = "Log-Likelihood",
                                  x_lab = "expression(lambda)"))

      }

    }

  } else if (inherits(label, "list")) {

    if (!any(names(label) %in% c("qq_res", "qq_ran",
                               "d_res", "d_ran",
                               "cooks", "opt_lambda")) ||
        !any(names(label) %in% c("qq_res", "qq_ran",
                                 "d_res", "d_ran",
                                 "cooks", "box_cox"))) {
     stop("List elements must have following names even though not
          all must be included: qq_res, qq_ran, d_res, d_ran, cooks,
          opt_lambda. Every list element must have the elements title,
          y_lab and x_lab. See also help(plot.saeTrafo).")
    }
    for (i in names(label)) {
      if (!all(names(label[[i]]) == c("title", "y_lab", "x_lab"))) {
        stop("Every list element must have the elements title,
              y_lab and x_lab in this order. See also
              help(plot.saeTrafo).")
      }
    }

      orig_label <- list(qq_res = c(title = "Error term",
                                    y_lab = "Quantiles of pearson residuals",
                                    x_lab = "Theoretical quantiles"),
                         qq_ran = c(title = "Random effect",
                                    y_lab = "Quantiles of random effects",
                                    x_lab = "Theoretical quantiles"),
                         d_res = c(title = "Density - Pearson residuals",
                                   y_lab = "Density",
                                   x_lab = "Pearson residuals"),
                         d_ran = c(title = "Density - Standardized random effects",
                                   y_lab = "Density",
                                   x_lab = "Standardized random effects"),
                         cooks = c(title = "Cook's Distance Plot",
                                   y_lab = "Cook's Distance",
                                   x_lab = "Index"),
                         opt_lambda = c(title = paste0(str_to_title(gsub("\\.","-", x$transformation)), ' - REML'),
                                     y_lab = "Log-Likelihood",
                                     x_lab = "expression(lambda)"))

      if (any(names(label) == "qq_res")) {
        label$qq_res <- label$qq_res
      } else {
        label$qq_res <- orig_label$qq_res
      }
      if (any(names(label) == "qq_ran")) {
        label$qq_ran <- label$qq_ran
      } else {
        label$qq_ran <- orig_label$qq_ran
      }
      if (any(names(label) == "d_res")) {
        label$d_res <- label$d_res
      } else {
        label$d_res <- orig_label$d_res
      }
      if (any(names(label) == "d_ran")) {
        label$d_ran <- label$d_ran
      } else {
        label$d_ran <- orig_label$d_ran
      }
      if (any(names(label) == "cooks")) {
        label$cooks <- label$cooks
      } else {
        label$cooks <- orig_label$cooks
      }
      if (any(names(label) == "opt_lambda") || any(names(label) == "box_cox")) {
        if (any(names(label) == "opt_lambda")) {
          label$opt_lambda <- label$opt_lambda
        } else if (any(names(label) == "box_cox")) {
          label$opt_lambda <- label$box_cox
        }
      }
      else {
        label$opt_lambda <- orig_label$opt_lambda
      }
  }

  if (any(!(names(label) %in%  c("qq_res", "qq_ran",
                               "d_res", "d_ran",
                               "cooks", "opt_lambda", "box_cox")))) {
    warning("One or more list elements are not called qq_res, qq_ran, d_res,
             d_ran, cooks or opt_lambda. The changes are for this/these
             element(s) is/are not done. Instead the original labels are used.")
  }

  return(label)
}


