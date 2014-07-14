library(raadtools)
        
datadir <- getOption("default.datadir")

fs <- list.files(file.path(datadir, "oceandata.sci.gsfc.nasa.gov"), recursive = TRUE, pattern = "DAY_CHL",
           full.names = TRUE)

sum(file.info(fs)$size / 1e6) / length(fs)


length(fs)