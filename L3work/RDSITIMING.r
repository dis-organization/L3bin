##rdsi timing

#system.time(x <- binlist("/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S1998001.L3b_DAY_CHL.main"))
##user  system elapsed 
##3.123 134.793 137.921 
#system.time(x <- binlist("L3work/S1998001.L3b_DAY_CHL.main"))
##user  system elapsed 
##1.835   1.102   2.937

dp <- "/tmp"
library(raadtools)
fs <- chlafiles()
ss <- sample(nrow(fs), 1L)
f <- fs$fullname[ss]

f <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov/SeaWiFS/L3BIN/1998/001/S1998001.L3b_DAY_RRS.main"
file.copy(f, file.path(dp, basename(f)))


system.time(x <- binlist(file.path(dp, basename(f))))
system.time(x2 <- binlist(f))

system.time(readAll(raster(f)))
system.time(readAll(raster(file.path(dp, basename(f)))))
            
