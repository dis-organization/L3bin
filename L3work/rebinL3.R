source("L3read.R")
dp <- "oceandata.sci.gsfc.nasa.gov/MODISA/L3BIN/2010"
datadir <- getOption("default.datadir")
f <- "A20100012010031.L3b_MO_CHL.main"
fs <- list.files(file.path(datadir, dp), recursive = TRUE, pattern = "MO_CHL")


list.files(file.path(datadir, "oceandata.sci.gsfc.nasa.gov/MODISA/L3BIN"))
bin <- readL3bin(file.path(datadir, dp, f))

library(raadtools)
x <- readchla("2010-01-01", time.resolution = "monthly", product = "oceancolor")
x1 <- aggregate(x, fact = 2.5/mean(res(x)), fun = mean, na.rm = TRUE)
chl <- chl.pal(palette = TRUE)
cols <- chl$cols
brks <- chl$breaks



cn <- extract(x1, cbind(lons, lats), cellnumbers = TRUE)

## mean of bin sums
aggvals <- tapply(bin$data[!is.na(cn[,1]), 1] / bin$weights[!is.na(cn[,1])], cn[!is.na(cn[,1]), 1] , mean)


y1 <- x1 * 0
y1[as.numeric(names(aggvals))] <- aggvals
par(mfrow = c(3, 1))
plot(x1, col = "white", legend = FALSE);plot(x, col = cols, breaks = brks, legend = FALSE, add = TRUE)
plot(x1, col= "white", legend = FALSE);plot(x1, col = cols, breaks = brks, legend = FALSE, add = TRUE)
plot(x1, col = "white", legend = FALSE);plot(y1, col = cols, breaks = brks, legend = FALSE, add = TRUE)
plot(values(y1), values(x1))
plot(crop(abs(x1 - y1), extent(-180, 180, -90, -30)), col = sst.pal(100))


binll <- L3lonlat(bin$bin_num, bin$start_num, bin$hsize)

px <- projectRaster(crop(x, extent(-180, 180, -90, 0)), crs = "+proj=stere +lat_ts=-71", res = c(5e4, 5e4))

asub <- lats <= 0
cn <- extract(px, project(cbind(lons[asub], lats[asub]), projection(px)), cellnumbers = TRUE)

## mean of bin sums
aggvals <- tapply(bin$data[asub, ][!is.na(cn[,1]), 1] / bin$weights[asub][!is.na(cn[,1])], cn[!is.na(cn[,1]), 1] , mean)
y1 <- px * 0
y1[as.numeric(names(aggvals))] <- aggvals

plot(y1, col = cols, breaks = brks)




asub <- lats <= -30
px <- projectRaster(crop(x, extent(-180, 180, -90, -30)), crs = "+proj=laea +lat_0=-90", res = c(1e5, 1e5))

cn <- extract(px, project(cbind(lons[asub], lats[asub]), projection(px)), cellnumbers = TRUE)


## bin sums
wtmu <- bin$data[asub, ][!is.na(cn[,1]), 1] / bin$weights[asub][!is.na(cn[,1])]
wtsxx <- (bin$data[asub, ][!is.na(cn[,1]), 2] / bin$weights[asub][!is.na(cn[,1])]) - (wtmu^2)


wtmu <- log(bin$data[asub, ][!is.na(cn[,1]), 1]) / bin$weights[asub][!is.na(cn[,1])]
wtsxx <- (log(bin$data[asub, ][!is.na(cn[,1]), 2]) / bin$weights[asub][!is.na(cn[,1])]) - (wtmu^2)

wtmle <- exp(1) ^ (wtmu + 0.5 * wtsxx)
wtsd <- wtmle * sqrt(exp(1)^wtsxx - 1)


aggmu <- tapply(wtmle, cn[!is.na(cn[,1]), 1], mean)
aggsxx <- tapply(wtsd, cn[!is.na(cn[,1]), 1], mean)

y2 <- y1 <- px * 0
y1[as.numeric(names(aggmu))] <- aggmu
y2[as.numeric(names(aggsxx))] <- aggsxx

par(mfrow = c(1, 2))
plot(y1, col = cols, breaks = brks, legend = FALSE)
plot(y2, col = sst.pal(100), legend = FALSE)





## can we reconstruct the L3mapped grid?
cn <- extract(x, cbind(lons, lats), cellnumbers = TRUE)

## mean of bin sums
aggvals <- tapply(bin$data[!is.na(cn[,1]), 1] / bin$weights[!is.na(cn[,1])], cn[!is.na(cn[,1]), 1] , mean)

y1 <- x * 0
y1[as.numeric(names(aggvals))] <- aggvals
