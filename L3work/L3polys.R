load("D:\\Toolbox\\data_candidates\\chlor\\L3BIN\\A20100012010031.L3b_MO_CHL.main.Rdata")

library(raster)
latstart <- -90

latstep <- 180 / 4320


lastID <- 0


l <- vector("list", max(bin$start_num))
for (i in seq_along(bin$start_num)[-1]) {
  nbins <- bin$start_num[i] - bin$start_num[i-1]
  p <- as(raster(xmn = -180, xmx = 180, ymn = -90 + (i-1) * latstep, ymx = -90 + i * latstep, nrows = 1, ncols = nbins)	, "SpatialPolygons")
  
  ind <- lastID + seq_len(nbins)
  p <- spChFIDs(p, as.character(ind))
  lastID <- lastID + nbins
  l[ind] <- p@polygons
  if (i %% 100 == 0) print(i)
  
}