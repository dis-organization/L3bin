library(L3bin)
##x <- binlist("/home/mdsumner/Git/L3bin/L3work/S1998001.L3b_DAY_CHL.main")

x <- binlist("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_CHL.main")

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



