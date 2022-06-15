#' Quantile-quantile plots for an saeTrafo object
#'
#' Normal quantile-quantile plots of the underlying model (see
#' \code{\link{NER_Trafo}}) are obtained. The plots are obtained by
#' \code{\link[ggplot2]{ggplot}}.
#' @param y a model object of type "saeTrafo" (see \code{\link{NER_Trafo}}).
#' @param color a character vector with two elements. The first element defines
#' the color for the line in the Q-Q plots, for the Cook's Distance plot and for
#' the Box-Cox plot. The second element defines the color for the densities.
#' @param gg_theme \code{\link[ggplot2]{theme}} list from package \pkg{ggplot2}.
#' For using this argument, package \pkg{ggplot2} must be loaded via
#' \code{library(ggplot2)}. See also Example 2.
#' @param ... optional arguments passed to generic function.
#' @return Two Q-Q plots in one grid obtained by \code{\link[ggplot2]{ggplot}}.
#' @seealso \code{\link{saeTrafoObject}}, \code{\link{NER_Trafo}}
#' @examples
#'
#' # Examples for Quantile-quantile plots
#'
#' # Load Data
#' data("eusilcA_smp")
#' data("pop_area_size")
#' data("pop_mean")
#' data("pop_cov")
#'
#' # Nested error regression model
#' NER_model <- NER_Trafo(fixed = eqIncome ~ gender + eqsize + cash +
#'                        self_empl + unempl_ben + age_ben + surv_ben +
#'                        sick_ben + dis_ben + rent + fam_allow + house_allow +
#'                        cap_inv + tax_adj,
#'                        smp_domains = "district",
#'                        pop_area_size = pop_area_size,
#'                        pop_mean = pop_mean, pop_cov = pop_cov,
#'                        smp_data = eusilcA_smp)
#'
#' # Example 1: Default Quantile-quantile plots
#' qqnorm(NER_model)
#'
#' # Example 2: Personalized plot using theme
#' require("ggplot2")
#' library(ggplot2)
#' qqnorm(NER_model,
#'        color = c("red", "darkgreen"),
#'        gg_theme = theme(panel.background = element_rect(fill = NA),
#'                         panel.grid.major = element_line(colour = "grey50"),
#'                         panel.ontop = TRUE)
#' )
#'
#' @export
#' @method qqnorm saeTrafo
#' @importFrom ggplot2 qplot geom_abline ggtitle ylab xlab ggplot stat_qq
#' @importFrom ggplot2 aes
#' @importFrom nlme ranef random.effects
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @importFrom stats sd residuals qnorm shapiro.test logLik cooks.distance

qqnorm.saeTrafo <- function(y, color = c("blue", "lightblue3"),
                            gg_theme = NULL, ...) {

  residuals <- residuals(y$model, level = 0, type = "pearson")
  tmp <- as.matrix(random.effects(y$model))[, 1]

  model <- y$model
  model$call$fixed <- y$fixed

  ## QQ Plots
  # Residuals
  res <- qplot(sample = residuals) +
    geom_abline(colour = color[1]) +
    ggtitle("Error term") +
    ylab("Quantiles of pearson residuals") +
    xlab("Theoretical quantiles") +
    gg_theme

  # Random effects
  ran <- ggplot(data.frame(tmp), aes(sample = tmp)) +
    stat_qq(distribution = qnorm, dparams = list(mean = mean(tmp),
                                                 sd = sd(tmp))) +
    geom_abline(intercept = 0, slope = 1, na.rm = TRUE, col = color[1]) +
    ggtitle("Random effect") +
    ylab("Quantiles of random effects") +
    xlab("Theoretical quantiles") +
    gg_theme

  invisible(grid.arrange(arrangeGrob(res, ran, ncol = 2)))
}
