## readL3.R
#' @export
readL3 <- function(x, vname) {
  x <- iz2(x)
  vdatalist <- vdatainfo(x)
  if (missing(vname)) {
    vdatalist <- vdatalist[!names(vdatalist) %in% c("SEAGrid", "BinList", "BinIndex")]
  } else {
    vdatalist <- vdatalist[names(vdatalist) == vname]
  }
  
  bindata <- binlist(x, names(vdatalist)[1L])
}

iz2 <- function(x) {
  uncx <- gsub(".bz2$", "", x)
  if (!file.exists(x)) {
    if (!file.exists(uncx)) stop(sprintf("no file: %s", x)) else needsdecompress <- TRUE
  }

  if (needsdecompress) {
    system(sprintf("bunzip2 %s", x))
    return(uncx)
  }
  x
}


#' Bin map
#' 
#' mapping between bins and a given raster
#' @export
binmap <- function(bin, ras, init = NULL) {
  if (is.null(init)) init <- initbin()
  ## TODO do this smarter
  ll <- SpatialPoints(do.call(cbind, bin2lonlat(bin)), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  extract(ras, ll, cellnumbers = TRUE)[,"cells"]
}


#' Initialize values for a particular binning
#' 
#' Set up the basic values for the bin scheme for given number of rows. 
#' @export 
initbin <- function(NUMROWS = 2160) {
## TODO options for lon-lat sub-sets
  latbin <- (((seq(NUMROWS) - 1) + 0.5) * 180 / NUMROWS ) - 90
  numbin <- as.integer(2 * NUMROWS * cos(latbin * pi/180) + 0.5)
  basebin <- cumsum(c(1L, numbin[-length(numbin)]))
  totbins = basebin[NUMROWS] + numbin[NUMROWS] - 1
  list(latbin = latbin, numbin = numbin, basebin = basebin, totbins = totbins)
}
# 
# bin2lonlat <- function(bin) {
#   row = NUMROWS - 1;
#   fint <- findInterval(bin, basebin)
#   clat = latbin[fint];
#   clon = 360.0*(bin - basebin[fint] + 0.5)/numbin[fint] - 180.0;
#   cbind(clon, clat)
# }
#' @export
bin2bounds <- function(bin, NUMROWS = 2160) {
  row = NUMROWS - 1;
  latbin <- (((seq(NUMROWS) - 1) + 0.5) * 180 / NUMROWS ) - 90
  numbin <- as.integer(2 * NUMROWS * cos(latbin * pi/180) + 0.5)
  basebin <- cumsum(c(1L, numbin[-length(numbin)]))
  fint <- findInterval(bin, basebin)
  north <- latbin[fint] + 90.0/NUMROWS
  south <- latbin[fint] - 90.0/NUMROWS
    ##*north = latbin[row] + 90.0/NUMROWS;
    ##*south = latbin[row] - 90.0/NUMROWS;
    lon = 360.0*(bin - basebin[fint] + 0.5)/numbin[fint] - 180.0;
    west = lon - 180.0/numbin[fint];
    east = lon + 180.0/numbin[fint];
list(east = east, south = south, west =   west, north = north)
  }


