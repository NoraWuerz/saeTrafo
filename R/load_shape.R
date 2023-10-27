#' Loading the shape file for Austrian districts
#'
#' The function simplifies to load the shape file for Austrian districts.
#'
#' @return A shape file of class \code{sf}.
#' @details The shape file contains the borders of Austrian districts. Thus, it
#' can be used for the visualization of estimation results for Austrian
#' districts.
#' @export

load_shapeaustria <- function() {
  load(file  = system.file("shapes/shape_austria_dis.rda",
                           package = "saeTrafo"),
       envir = .GlobalEnv
  )
}
