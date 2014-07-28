bchlfile <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_CHL.main"
brrsfile <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_RRS.main"
ocfile <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/9km/chlor/S19980011998031.L3m_MO_CHL_chlor_a_9km"
jfile <- "/rdsi/PRIVATE/chl/johnson/seawifs/monthly/S19980011998031.L3m_MO_SO_Chl_9km.Johnson_SO_Chl.nc"


library(L3bin)

##file.copy(bchlfile, basename(bchlfile))
##file.copy(brrsfile, basename(brrsfile))

bchlfile <- "S19980011998031.L3b_MO_CHL.main.bz2"
brrsfile <- "S19980011998031.L3b_MO_RRS.main.bz2"

b <- "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile"

download.file(file.path(b, bchlfile), bchlfile, mode = "wb")
download.file(file.path(b, brrsfile), brrsfile, mode = "wb")

system(sprintf("bunzip2 %s", bchlfile))
system(sprintf("bunzip2 %s", brrsfile))

bchlf <- gsub(".bz2", "", bchlfile)
brrsf <- gsub(".bz2", "", brrsfile)


bchl <- readL3(basename(bchlf))
brrs <- readL3(basename(brrsf)
               
               
               save(bchl, brrs, file = "rawMO.Rdata")
               
               
               library(raadtools)
               oc <- readchla("1998-01-15", time.resolution = "monthly", product = "oceancolor")
               jc <- readchla("1998-01-15", time.resolution = "monthly", product = "johnson")
               
               # 
               bin2lonlat <- function(bin) {
                 row = NUMROWS - 1;
                 fint <- findInterval(bin, basebin)
                 clat = latbin[fint];
                 clon = 360.0*(bin - basebin[fint] + 0.5)/numbin[fint] - 180.0;
                 cbind(clon, clat)
               }
               
               ll <- bin2lonlat(bchl$bin_num)
               ##ll <- cbind(ll$x, ll$y)
               
               
               
               
               e <- new("Extent"
                        , xmin = -80.1248170756666
                        , xmax = -58.1728095345417
                        , ymin = -61.5656963058103
                        , ymax = -40.5821596856174
               )
               new("Extent"
                   , xmin = 164.521045068331
                   , xmax = 175.748592647521
                   , ymin = -79.9102234492261
                   , ymax = -69.4311790419819
               )
               
               
               
     
               asub <- ll[,1] >= xmin(e) & ll[,1] <= xmax(e) & ll[,2] >= ymin(e) & ll[,2] <= ymax(e)
               bb <- bin2bounds(brrs$bin_num[asub])
               
               xy <- ll[asub, ]
               ocbin <- swchl(brrs)[asub]
               jbin <- swchl(brrs, johnson = TRUE)[asub]
               
               
               par(mfrow = c(2,2))
               plot(crop(oc, e), col = pal$cols, breaks = pal$breaks, asp = cos(-50/pi/180), legend = FALSE)
               plot(crop(oc, e), col = "white", breaks = pal$breaks, asp = cos(-50/pi/180), legend = FALSE)
               ##plot(xy, type = "n", asp = cos(-50/pi/180))
               rect(bb$east, bb$south, bb$west, bb$north, col = chl.pal(ocbin), border = NA)
               
               plot(crop(jc, e), col = pal$cols, breaks = pal$breaks, asp = cos(-50/pi/180), legend = FALSE)
               plot(crop(oc, e), col = "white", breaks = pal$breaks, asp = cos(-50/pi/180), legend = FALSE)
               
               ##plot(xy, type = "n", asp = cos(-50/pi/180))
               rect(bb$east, bb$south, bb$west, bb$north, col = chl.pal(jbin), border = NA)
               
               
               
               
               
               
               
               




















## 
##hdir <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001"
##f1 <- "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S19980011998008.L3b_8D_CHL.main.bz2"
##download.file(f1, file.path(hdir, basename(f1)), mode = "wb")

## S1998003.L3b_DAY_CHL.main.bz2"
## prototype L3 re-binning

##f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998008.L3b_8D_CHL.main"
f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998031.L3b_MO_CHL.main"

##f = "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S19980011998008.L3b_8D_CHL.main"
library(raadtools)
##ochl <- readchla("1998-01-01", product = "oceancolor", time.resolution = "monthly")
##ochl <- raster("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/9km/chlor/S19980011998031.L3m_MO_CHL_chlor_a_9km")
##ochl <- raster("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/Mapped/Monthly/9km/chlor/S19980011998008.L3m_8D_CHL_chlor_a_9km")
ochl <- readchla("1998-01-01", product = "oceancolor", time.resolution = "monthly")
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
  row[row < 1] <- 1
  col = trunc((x[,1] + 180.0)*init$numbin[row]/360.0);
  ##if(col >= numbin[row]) col = numbin[row] - 1;
  col <- ifelse(col > init$numbin[row], init$numbin[row] - 1, col)
  col[col < 1] <- 1
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



## wrong, this must aggregate the bins first 
rebin <- readice()
res(rebin) <- res(rebin) * 12
rebin[] <- 0
recells <- extract(rebin, ll, cellnumbers = TRUE)

rebin[unique(recells[,1])] <- tapply((rawbins$sum/rawbins$weight), recells[,1], FUN = mean)

plot(rebin, col = pal$cols, breaks = pal$breaks, legend = FALSE)



rebin <- readsst()
res(rebin) <- 2.5
rebin[] <- 0
library(rgdal)
##recells <- extract(rebin, project(ll, projection(rebin)), cellnumbers = TRUE)
recells <- extract(rebin, ll, cellnumbers = TRUE)
notbad <- !is.na(recells[,1])

rebin[unique(recells[notbad,1])] <- tapply(rawbins$sum[notbad] / rawbins$weights[notbad], recells[notbad,1], FUN = function(x) exp(1)^(mean(log(x), na.rm = TRUE)))
plot(rebin, col = pal$cols, breaks = pal$breaks)

##cnt[cells[,1]] <- cnt[cells[,1]] + 1
##print(i)




rebin <- readice()
res(rebin) <- res(rebin) * 2.8
rebin[] <- 0

rebin <- projectExtent(crop(ochl, extent(-180, 180, -90, 30)), projection(readice()))

 res(rebin) <- c(1e5)
rebin[] <- 0

recells <- extract(rebin, project(ll, projection(rebin)), cellnumbers = TRUE)

rebin[unique(recells[,1])] <- tapply((rawbins$sum/rawbins$weight), recells[,1], FUN = mean)

plot(rebin, col = pal$cols, breaks = pal$breaks, legend = FALSE)

ll <-  coordinates(rebin)
ll <- project(ll, projection(rebin), inv = TRUE)
binid <- latlon2bin(ll)

rebin[] <- rawbins$sum[match(binid, rawbins$bin_num)] / rawbins$weights[match(binid, rawbins$bin_num)]

plot(rebin, col = pal$cols, breaks = pal$breaks)





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