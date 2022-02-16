writeexcel_check <- function(object, file, split) {

  if (!inherits(object, "saeTrafo")) {
    stop("First object needs to be of class saeTrafo.")
  }
  if (!(inherits(split, "logical") && length(split) == 1)) {
    stop("split must to be a logical value. Set CV to TRUE or FALSE.")
  }
  if (!inherits(file, "character") || (length(grep(".xlsx", file)) == 0)) {
    stop(paste("file must to be a character string that determines a path and",
               "filename. It should end on .xlsx"))
  }
}

writeods_check <- function(object, file, split) {

  if (!inherits(object, "saeTrafo")) {
    stop("First object needs to be of class saeTrafo.")
  }
  if (!(inherits(split, "logical") && length(split) == 1)) {
    stop("split must to be a logical value. Set CV to TRUE or FALSE.")
  }
  if (!inherits(file, "character") || (length(grep(".ods", file)) == 0)) {
    stop(paste("file must to be a character string that determines a path and",
               "filename. It should end on .ods"))
  }
  testZIP <- 1
  try(testZIP <- shell("zip"), silent = TRUE)
  if (testZIP == 1) {
    stop(paste("No zipping application was found, see details in",
               "help(write.ods) for details."))
  }
}
