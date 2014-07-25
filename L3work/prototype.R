## 
##hdir <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001"
##f1 <- "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S19980011998008.L3b_8D_CHL.main.bz2"
##download.file(f1, file.path(hdir, basename(f1)), mode = "wb")

## S1998003.L3b_DAY_CHL.main.bz2"
## prototype L3 re-binning

##f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998008.L3b_8D_CHL.main"
##f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_CHL.main"

f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998008.L3b_8D_CHL.main"
library(raadtools)
##ochl <- readchla("1998-01-01", product = "oceancolor", time.resolution = "monthly")
##ochl <- raster("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/9km/chlor/S19980011998031.L3m_MO_CHL_chlor_a_9km")
##ochl <- raster("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/9km/chlor/S19980011998008.L3m_8D_CHL_chlor_a_9km")
ochl <- readchla("1998-01-01", product = "oceancolor", time.resolution = "weekly")
extent(ochl) <- extent(-180, 180, -90, 90)
projection(ochl) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
ochl[ochl <= 0] <- NA
##f <- file.path(hdir, basename(f1))
##f <- L3bin:::iz2(f)

rbin <- raster(ochl)
rbin[] <- 0

rawbins <- readL3(f)
ll <- do.call(cbind, bin2lonlat(rawbins$bin_num))
cells <- extract(ochl, ll, cellnumbers = TRUE)
rbin[cells[,1]] <- rbin[cells[,1]] + rawbins$sum / rawbins$weights

e1 <- extent(9.75, 19.1, -42, -33) - c(7, 7.2)
e2 <- extent(-59, -44, -77, -73) - c(8.5, 2.1)
asub1 <- ll[,1] >= xmin(e1) & ll[,1] <= xmax(e1) & ll[,2] >= ymin(e1) & ll[,2] <= ymax(e1)

asub2 <- ll[,1] >= xmin(e2) & ll[,1] <= xmax(e2) & ll[,2] >= ymin(e2) & ll[,2] <= ymax(e2)
bb1 <- bin2bounds(rawbins$bin_num[asub1])
bb2 <- bin2bounds(rawbins$bin_num[asub2])


pal <- chl.pal(palette = TRUE)
pal$cols[1] <- "transparent"

cr1_ochl <- crop(ochl, e1, snap = "out")
cr1_rbin <- crop(rbin, e1, snap = "out")
cr2_ochl <- crop(ochl, e2, snap = "out")
cr2_rbin <- crop(rbin, e2, snap = "out")


latlon2bin <- function(x) {
  row <- findInterval(x[,2], init$latbin)
  col = trunc((x[,1] + 180.0)*init$numbin[row]/360.0);
  ##if(col >= numbin[row]) col = numbin[row] - 1;
  col <- ifelse(col > init$numbin[row], init$numbin[row] - 1, col)
  ##col[col >= numbin[row] <- numbin[row] - 1;
  init$basebin[row] + col;
  
}


rebin <- cr2_rbin
rebin[] <- 0
binid <- latlon2bin(coordinates(rebin))
rebin[] <- rawbins$sum[match(binid, rawbins$bin_num)] / rawbins$weights[match(binid, rawbins$bin_num)]

recells <- extract(cr2_rbin, ll[asub2,], cellnumbers  = TRUE)
notbad <- !is.na(recells[,1])
rebin[unique(recells[notbad,1])] <- tapply(rawbins$sum[asub2][notbad] / rawbins$weights[asub2][notbad], recells[notbad,1], FUN = function(x) exp(1)^(mean(log(x), na.rm = TRUE)))

rr <- rasterize(data.frame)
png("afile.png", 1200, 2000)
#pdf("afile.pdf")
op <- par(mfcol = c(3,2))
# plot(cr1_ochl, col = chl.pal(values(cr1_ochl)))
# plot(cr1_rbin, col = chl.pal(values(cr1_rbin)))
# plot(cr2_ochl, col = chl.pal(values(cr2_ochl)))
# plot(cr2_rbin, col = chl.pal(values(cr2_rbin)))
plot(cr1_ochl, col = pal$cols, breaks = pal$breaks, legend = FALSE)
plot(cr1_rbin, col = pal$cols, breaks = pal$breaks, legend = FALSE)
plot(cr2_ochl, col = pal$cols, breaks = pal$breaks, legend = FALSE)
plot(cr2_rbin, col = pal$cols, breaks = pal$breaks, legend = FALSE)

plot(cr1_ochl)
rect(bb1$east, bb1$south, bb1$west, bb1$north, col = chl.pal((rawbins$sum/rawbins$weight)[asub1]), border = "transparent")
plot(cr2_rbin)
rect(bb2$east, bb2$south, bb2$west, bb2$north, col = chl.pal((rawbins$sum/rawbins$weight)[asub2]), border = "transparent")
##par(op)
dev.off()

rebin <- readsst()
rebin[] <- 0
library(rgdal)
##recells <- extract(rebin, project(ll, projection(rebin)), cellnumbers = TRUE)
recells <- extract(rebin, ll, cellnumbers = TRUE)
notbad <- !is.na(recells[,1])

rebin[unique(recells[notbad,1])] <- tapply(rawbins$sum[notbad] / rawbins$weights[notbad], recells[notbad,1], FUN = function(x) exp(1)^(mean(log(x), na.rm = TRUE)))
plot(rebin, col = pal$cols, breaks = pal$breaks)

##cnt[cells[,1]] <- cnt[cells[,1]] + 1
##print(i)



## FINAL REBIN to match original
rebin <- ochl
rebin[] <- 0
binid <- latlon2bin(coordinates(rebin))
rebin[] <- rawbins$sum[match(binid, rawbins$bin_num)] / rawbins$weights[match(binid, rawbins$bin_num)]





fs <- list.files("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998", full.names = TRUE, recursive = TRUE, pattern = "DAY_CHL")[1:8]
library(raadtools)
chl <- readchla("1998-01-01", product = "oceancolor")


r <- raster(chl)
r[1:ncell(r)] <- 0
cnt <- r

for (i in seq_along(fs)) {
  x <- readL3(fs[i])
  ll <- bin2lonlat(x$bin_num)
  
  cells <- extract(chl, do.call(cbind, ll), cellnumbers = TRUE)

  r[cells[,1]] <- r[cells[,1]] + x$sum / x$weights
  cnt[cells[,1]] <- cnt[cells[,1]] + 1
  print(i)
  plot(log(r))
}

r <- r/cnt

pal <- chl.pal(palette = TRUE)

plot(chl, col = pal$cols, breaks = pal$breaks, legend = FALSE)