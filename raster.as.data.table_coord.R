#' Transform raster to data.table including the coordinates, available at https://gist.github.com/etiennebr/9515738
#' 
#' @param x  Raster* object
#' @param row.names	`NULL` or a character vector giving the row names for the data frame. Missing values are not allowed
#' @param optional	logical. If `TRUE`, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param xy  logical. If `TRUE`, also return the spatial coordinates
#' @param centroids	logical. If TRUE return the centroids instead of all spatial coordinates (only relevant if xy=TRUE)
#' @param sepNA	logical. If TRUE the parts of the spatial objects are separated by lines that are NA (only if xy=TRUE and, for polygons, if centroids=FALSE
#' @param ...	 Additional arguments (none) passed to `raster::as.data.frame`
#' 
#' @value returns a data.table object
#' @examples
#' logo <- brick(system.file("external/rlogo.grd", package="raster"))
#' v <- as.data.table(logo)
#' @import

as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
    coln <- names(x)
    
    if(xy) coln <- c("x", "y", coln)
    setnames(v, coln)
  } else {
    tr <- blockSize(x)
    l <- lapply(1:tr$n, function(i) {
      DT <- as.data.table(as.data.frame(getValues(x, row = tr$row[i], nrows = tr$nrows[i]), ...))  
      if(xy == TRUE) {
        cells <- cellFromRowCol(x, c(tr$row[i], tr$row[i] + tr$nrows[i] - 1), c(1, ncol(x)))
        coords <- xyFromCell(x, cell = cells[1]:cells[2])
        DT[, c("x", "y") := data.frame(xyFromCell(x, cell = cells[1]:cells[2]))]
      } 
      DT
    })
    v <- rbindlist(l)
    coln <- names(x)
    if(xy) {
      coln <- c("x", "y", coln)
      setcolorder(v, coln)
    }
  }
  v
}

if (!isGeneric("as.data.table")) {
  setGeneric("as.data.table", function(x, ...)
    standardGeneric("as.data.table"))
}  

setMethod('as.data.table', signature(x='data.frame'), data.table::as.data.table)
setMethod('as.data.table', signature(x='Raster'), as.data.table.raster)




