readL3bin <- function(x) {
   if (!file.exists(x)) stop(sprintf("file does not exist: %s", x))
   binfile <- tempfile()
   datfile <- tempfile()
   stfile <- tempfile()
   hfile <- tempfile()
   wfile <- tempfile()
   call1 <- sprintf("hdp dumpvd -d -f bin_num -n BinList -b -o %s %s", binfile, x)
   call2 <- sprintf("hdp dumpvd -d -n chlor_a -b -o %s %s", datfile, x)
   call3 <- sprintf("hdp dumpvd  -f start_num -d -n BinIndex -b -o %s %s", stfile, x)
   call4 <- sprintf("hdp dumpvd  -f hsize -d -n BinIndex -b -o %s %s", hfile, x)
   call5 <- sprintf("hdp dumpvd -d -f weights -n BinList -b -o %s %s", wfile, x)
   system(call1)
   system(call2)
   system(call3)
   system(call4)
   system(call5)
   bsize <- file.info(binfile)$size
   dsize <- file.info(datfile)$size
   ssize <- file.info(stfile)$size
   hsize <- file.info(hfile)$size
   wsize <- file.info(wfile)$size
   binlist <- readBin(binfile, "integer", n = bsize /4, size = 4)
   datlist <- matrix(readBin(datfile, "numeric", n = dsize / 4, size = 4 ), ncol = 2, byrow = TRUE)
   stlist <- readBin(stfile, "integer", n = ssize/4, size = 4)
   hlist <- readBin(hfile, "numeric", n = hsize/8, size = 8)
   wlist <- readBin(wfile, "numeric", n = wsize/4, size = 4)
   unlink(binfile)
   unlink(datfile)
   unlink(stfile)
   unlink(hfile)
   unlink(wfile)
   list(bin_num = binlist, data = datlist, start_num = stlist, hsize = hlist, weights = wlist)
 }

bindex <- readL3bin("/home/mdsumner/Git/L3bin/L3work/S1998001.L3b_DAY_CHL.main")
L3lonlat <- function(bin_num, start_num, hsize) {
  latbase <- seq(-89.958336, 89.958336, length = length(start_num))
  boffset <- findInterval(bin_num, start_num)
  lats <- latbase[boffset]
  lons <- -180 + (bin_num - start_num[boffset]) * hsize[boffset]
  cbind(lons, lats)  
}
