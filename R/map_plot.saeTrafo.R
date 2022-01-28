# External documentation -------------------------------------------------------

#' Visualizes regional disaggregated estimates on a map
#'
#' Function \code{map_plot} creates spatial visualizations of the estimates
#' obtained by small area estimation methods or direct estimation.
#'
#' @param object an object of type emdi, containing the estimates to be
#' visualized.
#' @param MSE optional logical. If \code{TRUE}, the MSE is also visualized.
#' Defaults to \code{FALSE}.
#' @param CV optional logical. If \code{TRUE}, the CV is also visualized.
#' Defaults to \code{FALSE}.
#' @param map_obj an \code{SpatialPolygonsDataFrame} object as defined by the
#' \pkg{sp} package on which the data should be visualized.
#' @param map_dom_id a character string containing the name of a variable in
#' \code{map_obj} that indicates the domains.
#' @param map_tab a \code{data.frame} object with two columns that match the
#' domain variable from the census data set (first column) with the domain
#' variable in the map_obj (second column). This should only be used if the IDs
#' in both objects differ.
#' @param color a \code{vector} of length 2 defining the lowest and highest
#' color in the plots.
#' @param scale_points a structure defining the lowest, the mid and the highest
#' value of the colorscale. If a numeric vector of length two is given, this scale
#' will be used for every plot. Alternatively, a list defining colors for each
#' plot separately may be given.
#' @param guide character passed to
#' \code{scale_colour_gradient} from \pkg{ggplot2}.
#' Possible values are "none", "colourbar", and "legend".
#' @param return_data if set to \code{TRUE}, a fortified data frame including the
#' map data as well as the chosen indicators is returned. Customized maps can
#' easily be obtained from this data frame via the package \pkg{ggplot2}. Defaults
#' to \code{FALSE}.
#' @return Creates the plots demanded, and, if selected, a fortified data.frame
#' containing the mapdata and chosen indicators.
#' @seealso \code{\link{direct}}, \code{\link{ebp}}, \code{\link{fh}}, \code{\link{emdiObject}},
#' \code{\link[maptools]{readShapePoly}}
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes geom_polygon facet_wrap fortify coord_equal labs
#' @importFrom ggplot2 theme element_blank guides scale_fill_gradient
#' @importFrom ggplot2 scale_colour_gradient

map_plot <- function(object,
                     MSE = FALSE,
                     CV = FALSE,
                     map_obj = NULL,
                     map_dom_id = NULL,
                     map_tab = NULL,
                     color = c("white", "red4"),
                     scale_points = NULL,
                     guide = "colourbar",
                     return_data = FALSE
                     ){

  indicator <- c("Mean")

  if (is.null(map_obj)) {
    message("No Map Object has been provided. An artificial polygone is used for
            visualization")
    map_pseudo(object = object , indicator = indicator, panelplot = FALSE, MSE = MSE,
               CV = CV)
  } else if (class(map_obj) != "SpatialPolygonsDataFrame" ||
             attr(class(map_obj), "package") != "sp") {
    stop("map_obj is not of class SpatialPolygonsDataFrame from the sp package")
  } else {
    if (length(color) != 2 || !is.vector(color)) {
      stop("col needs to be a vector of length 2
           defining the starting, mid and upper color of the map-plot")
    }

    plot_real(object,
              indicator = indicator,
              MSE = MSE,
              CV = CV,
              map_obj = map_obj,
              map_dom_id = map_dom_id,
              map_tab = map_tab,
              col = color,
              scale_points = scale_points,
              return_data = return_data,
              guide = guide
    )
  }
}

map_pseudo <- function(object, indicator, panelplot, MSE, CV)
{
  x <- y <- id <- value <- NULL #avoid note due to usage in ggplot
  values <-  estimators(object = object, indicator = indicator,
                                 MSE = MSE, CV = CV)$ind
  indicator <- colnames(values)[-1]

  tplot <- get_polygone(values = values)

  if (panelplot) {
    ggplot2::ggplot(tplot, aes(x = x, y = y)) + ggplot2::geom_polygon(aes(group=id, fill = value))
    + ggplot2::facet_wrap( ~ variable, ncol = ceiling(sqrt(length(unique(tplot$variable)))))
  } else {
    for (ind in indicator) {
      print(ggplot2::ggplot(tplot[tplot$variable == ind,], aes(x = x, y = y)) +
              ggplot2::ggtitle(paste0(ind)) + ggplot2::geom_polygon(aes(group = id, fill = value)) )
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
}

plot_real <- function(object,
                      indicator = "all",
                      MSE = FALSE,
                      CV = FALSE,
                      map_obj = NULL,
                      map_dom_id = NULL,
                      map_tab = NULL,
                      col = col,
                      scale_points = NULL,
                      return_data = FALSE,
                      guide = NULL
) {
  if (!is.null(map_obj) && is.null(map_dom_id)) {
    stop("No Domain ID for the map object is given")
  }
  long <- lat <- group <- NULL


  map_data <- estimators(object = object, indicator = indicator,
                                  MSE = MSE, CV = CV)$ind

  if (!is.null(map_tab)) {
    map_data <- merge(x = map_data, y = map_tab,
                      by.x = "Domain", by.y = names(map_tab)[1])
    matcher <- match(map_obj@data[map_dom_id][,1], map_data[,names(map_tab)[2]])

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop("Domains of map_tab and Map object do not match. Check map_tab")
      } else {
        warnings("Not all Domains of map_tab and Map object could be matched.
                 Check map_tab")
      }
    }
    map_data <- map_data[matcher,]
    map_data <- map_data[,!colnames(map_data) %in% c("Domain",
                                                     map_dom_id,
                                                     names(map_tab)), drop = F]
  } else {
    matcher <- match(map_obj@data[map_dom_id][,1], map_data[,"Domain"])

    if (any(is.na(matcher))) {
      if (all(is.na(matcher))) {
        stop("Domain of EMDI and Map object do not match. Try using map_tab")
      } else {
        warnings("Not all Domains of EMDI and Map object could be matched.
                 Try using map_tab")
      }
    }
    map_data <- map_data[matcher, ]
  }

  map_obj@data[colnames(map_data)] <- map_data

  map_obj.fort <- ggplot2::fortify(map_obj, region = map_dom_id)
  map_obj.fort <- merge(map_obj.fort, map_obj@data,
                        by.x = "id", by.y = map_dom_id)

  indicator <- colnames(map_data)
  indicator <- indicator[!(indicator %in% "Domain")]

  for (ind in indicator) {
    map_obj.fort[ind][,1][!is.finite(map_obj.fort[ind][,1])] <- NA
    scale_point <- get_scale_points(map_obj.fort[ind][,1], ind, scale_points)
    print(ggplot2::ggplot(map_obj.fort, aes(long, lat, group = group,
                                   fill = map_obj.fort[ind][,1])) +
            ggplot2::geom_polygon(color = "azure3") + coord_equal() +
            ggplot2::labs(x = "", y = "", fill = ind) +
            ggplot2::ggtitle(gsub(pattern = "_",replacement = " ",x = ind)) +
            ggplot2::scale_fill_gradient(low = col[1], high = col[2],limits = scale_point,
                                guide = guide) +
            ggplot2::theme(axis.ticks = element_blank(), axis.text = element_blank(),
                  legend.title = element_blank())

    )
    if (!ind == tail(indicator,1)) {
      cat("Press [enter] to continue")
      line <- readline()
    }
  }
  if (return_data) {
    return(map_obj.fort)
  }
}

get_polygone <- function(values) {
  if (is.null(dim(values))) {
    values = as.data.frame(values)
  }
  n <- nrow(values)
  cols <- ceiling(sqrt(n))
  n <- cols^2

  values["id"] <- seq_len(nrow(values))

  poly <- data.frame(id = rep(seq_len(n), each = 4),
                     ordering = seq_len((n*4)),
                     x = c(0,1,1,0) + rep(0:(cols - 1), each = (cols * 4)),
                     y = rep(c(0,0,1,1) + rep(0:(cols - 1), each = 4), cols)
  )

  combo <- merge(poly, values, by = "id", all = TRUE, sort = FALSE)
  melt(combo[order(combo$ordering),], id.vars = c("id","x","y","ordering"))
}

get_scale_points <- function(y, ind, scale_points){
  result <- NULL
  if (!is.null(scale_points)) {
    if (class(scale_points) == "numeric" && length(scale_points) == 2) {
      result <- scale_points
    } else {
      splt <- strsplit(ind, "_\\s*(?=[^_]+$)", perl = TRUE)[[1]]
      indicator_name <- splt[1]
      if (length(splt) == 2) {
        measure <- splt[2]
      } else {
        measure <- "ind"
      }
      if (indicator_name %in% names(scale_points)) {
        pointset <- scale_points[[indicator_name]]
        try(result <- pointset[[measure]])
      }
      if (is.null(result) || length(result) != 2)
      {
        warning("scale_points is of no apropriate form, default values will
                 be used. See the descriptions and examples for details")
        result <- NULL
      }
    }
  }
  if (is.null(result)) {
    rg <- range(y, na.rm = TRUE)
    result <- rg
  }
  return(result)
}

