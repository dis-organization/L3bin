library(L3bin)
##x <- binlist("/home/mdsumner/Git/L3bin/L3work/S1998001.L3b_DAY_CHL.main")

##x <- binlist("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_CHL.main")

system.time(x <- binlist("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S1998001.L3b_DAY_CHL.main"))

system.time(x <- binlist("S1998001.L3b_DAY_CHL.main"))



load("~/Git/L3bin/L3work/bindex.Rdata")
L3lonlat <- function(bin_num, start_num, hsize) {
  latbase <- seq(-89.958336, 89.958336, length = length(start_num))
  boffset <- findInterval(bin_num, start_num)
  lats <- latbase[boffset]
  lons <- -180 + (bin_num - start_num[boffset]) * hsize[boffset]
  cbind(lons, lats)  
}

ll <- L3lonlat(x$bin_num, bindex$start_num, bindex$hsize)
library(raadtools)
chl <- readchla("1998-01-15", time.resolution = "monthly", product = "oceancolor")
pal <- chl.pal(palette = TRUE)
e <- extent(140, 150, 70, -60)
schl <- crop(chl, e)
plot(schl, col = pal$cols, breaks = pal$breaks, legend = FALSE)

ind <- extract(schl, ll, cellnumbers = TRUE)[,"cells"]
bad <- is.na(ind)
vals <- x$sum[!bad] / x$weights[!bad]

ind <- ind[!bad]

plot(ll[!bad, ], col = chl.pal(vals), pch = 16)



## all monthly files
fs <- list.files("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov", pattern = "L3b_MO_CHL", full.names = TRUE, recursive = TRUE)
files <- data.frame(fullname = fs, stringsAsFactors = FALSE)
  
files$date <- as.POSIXct(as.Date(substr(basename(files$fullname), 2, 8), "%Y%j"), tz = "GMT")
ord <- order(files$date)
files <- files[ord, ]  
  
## unique files
## MODIS trumps SeaWiFS
files <- files[!duplicated(files$date), ]

library(L3bin)
d <- binlist("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/MODISA/L3BIN/2002/182//A20021822002212.L3b_MO_CHL.main")
load("~/Git/L3bin/L3work/bindex.Rdata")
L3lonlat <- function(bin_num, start_num, hsize) {
  latbase <- seq(-89.958336, 89.958336, length = length(start_num))
  boffset <- findInterval(bin_num, start_num)
  lats <- latbase[boffset]
  lons <- -180 + (bin_num - start_num[boffset]) * hsize[boffset]
  cbind(lons, lats)  
}

ll <- L3lonlat(d$bin_num, bindex$start_num, bindex$hsize)

library(raadtools)
isub <- sample(nrow(ll), 1e4)
##plot(ll[isub, ], col = chl.pal(d$sum[isub]/d$weights[isub]), pch = ".")

NUMROWS <- 26
latbin <- (((seq(NUMROWS) - 1) + 0.5) * 180 / NUMROWS ) - 90
numbin <- as.integer(2 * NUMROWS * cos(latbin * pi/180) + 0.5)
basebin <- cumsum(c(1L, numbin[-length(numbin)]))
totbins = basebin[NUMROWS] + numbin[NUMROWS] - 1;

bin2lonlat <- function(bin) {
  row = NUMROWS - 1;
  fint <- findInterval(bin, basebin)
  clat = latbin[fint];
  clon = 360.0*(bin - basebin[fint] + 0.5)/numbin[fint] - 180.0;
  cbind(clon, clat)
}

bin2bounds <- function(bin) {
  
  fint <- findInterval(bin, basebin)
  north = latbin[fint] + 90.0/NUMROWS;
  south = latbin[fint] - 90.0/NUMROWS;
  lon = 360.0*(bin - basebin[fint] + 0.5)/numbin[fint] - 180.0;
  west = lon - 180.0/numbin[fint];
  east = lon + 180.0/numbin[fint];
  c(west, south, east, north)
}

l <- lapply(1:totbins, bin2bounds)
rect2P <- function(x) Polygon(as.matrix(expand.grid(x[c(1, 3)], x[c(2, 4)]))[c(1, 3, 4, 2, 1), ])


##zap <- lapply(l, function(x) rect(x[1], x[2], x[3], x[4]))
zap <- lapply(l, function(x) rect2P(x))
ps <- SpatialPolygonsDataFrame(SpatialPolygons(lapply(seq(along = zap), function(x) Polygons(zap[x], as.character(x))), proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")), data.frame(x = 1:length(zap)))


       


