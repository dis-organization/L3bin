baseurl <- "http://oceandata.sci.gsfc.nasa.gov"
sensor <- c("SeaWiFS", "MODISA", "VIIRS")
level <- "L3BIN"
years <- sprintf("%i", seq(1997, as.numeric(format(Sys.Date(), "%Y"))))



days <- sprintf( "%.03i", seq(1, 366))
## YR only so 
##days <- "001"
## MO only so 
##days <- c("001", "032", "060", "091", "121", "152", "182", "213", "244", "274", "305", "305")

## VIIRS needs more thought "MO_*.CHL" etc. doesn't match properly local/remote
filter <- "DAY_CHL"
##filter <- "YR_CHL"
baselocal <- "/rdsi/PRIVATE/oceandata.sci.gsfc.nasa.gov"

downloadyesno <- function(x) {
  ## return TRUE if file or file without extension does not exist
  test <- TRUE
  if (file.exists(x)) test <- FALSE
  if (file.exists(gsub(".bz2", "", x))) test <- FALSE
  test
}

## TODO:
##  this urltxt scan should be done once, comprehensively to save all file names and structure
## use it to build a filenames database
## query is filename where [time-bin], [within-date], [within-bin-index], [product-var]/[sensor]

errcon <- file(sprintf("%s_err.txt", format(Sys.Date())), open = "wt")
for (isensor in seq_along(sensor)) {
  for (iyear in seq_along(years)) {
    for (iday in seq_along(days)) {
     sourceurl <- file.path(baseurl, sensor[isensor], level[1L], years[iyear], days[iday])
        localfolder <- file.path(baselocal, sensor[isensor], level[1L], years[iyear], days[iday])
        urltxt <- try(readLines(sourceurl, warn = FALSE), silent = TRUE)
        if (!inherits(urltxt, "try-error")) {
          remotefiles <- sapply(strsplit(grep("getfile", urltxt, value = TRUE), "\""), "[", 2L)
          remotefiles <- grep(filter, remotefiles, value = TRUE)
         ## print(localfolder)
        ##  print(grep(filter, basename(remotefiles), value = TRUE))
          for (igetfile in seq_along(remotefiles)) {
            localfilepath  <- file.path(localfolder, basename(remotefiles[igetfile]))    
              if (downloadyesno(localfilepath)) {
                if (!file.exists(dirname(localfilepath))) dir.create(dirname(localfilepath), recursive = TRUE)
                
                print(remotefiles[igetfile])
                print(localfilepath)
                test <- try(download.file(remotefiles[igetfile], localfilepath, mode = "wb"))
                if (inherits(test, "try-error")) writeLines(sprintf("problem with: %s", remotefiles[igetfile]), errcon)
             
                
              }
          }
       
        
        }
      }
  }
}

close(errcon)



