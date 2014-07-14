baseurl <- "http://oceandata.sci.gsfc.nasa.gov"
sensor <- c("SeaWiFS", "MODISA")
level <- "L3BIN"

parseFolders <- function(x, ndigit = 4) {
    tx <- readLines(x, warn = FALSE)
    sapply(strsplit(grep(sprintf("href=\"[0-9]{%i}", ndigit), tx, value = TRUE), "\""), "[", 2)
}

# 
# ## base folder for this sensor
# x <- file.path(baseurl, sensor[isensor], level)
# ## year folders
# years <- parseFolders(x)
# xx <- file.path(x, years)
# iday <- 1
# dayfolders <- parseFolders(xx[iday], ndigit = 3)
# xxx <- file.path(xx, dayfolders)


## VIIRS needs more thought "MO_*.CHL" etc. doesn't match properly local/remote
##filter <- "DAY_CHL"
filter <- "YR_"
daysfilter <- "001"
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
## NEED to filter out the multiple RRS xNN variants

errcon <- file(sprintf("%s_err.txt", format(Sys.Date())), open = "wt")
for (isensor in seq_along(sensor)) {
  ## base folder for this sensor
  x <- file.path(baseurl, sensor[isensor], level)
  ## year folders
  years <- parseFolders(x)
 for (iyear in seq_along(years)) {
   
   xx <- file.path(x, years[iyear])
   days <- parseFolders(xx, ndigit = 3)
  
   days <- grep(daysfilter, days, value = TRUE)
    for (iday in seq_along(days)) {
      ## this is years so kick things along after first day
  
      sourceurl <- file.path(xx, days[iday])
      
##     sourceurl <- file.path(baseurl, sensor[isensor], level[1L], years[iyear], days[iday])
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
                ##scan("", 1)
                test <- try(download.file(remotefiles[igetfile], localfilepath, mode = "wb"))
                if (inherits(test, "try-error")) writeLines(sprintf("problem with: %s", remotefiles[igetfile]), errcon)
             
                
              }
          }
       
        
        }
      }
  }
}

close(errcon)



