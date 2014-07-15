## http://oceancolor.gsfc.nasa.gov/DOCS/FormatChange.html

## FormatChange

fpath <- "/rdsi/PRIVATE/scratch/mdsumner/L3bin"
dir.create(fpath)

## L3bin
bin <- c("http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S20030602003090.L3b_MO_ST92_CHL.nc", 
"http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20030602003090.L3b_MO_AT108_CHL.nc")

for (i in seq(along = bin)) download.file(bin[i], file.path(fpath, basename(bin[i])), mode = "wb")

## L3smi
##http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/S20030602003090.L3m_MO_ST92_CHL_chlor_a_9km.nc
##htp://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A20030602003090.L3m_MO_AT108_CHL_chlor_a_4km.nc

fs <- file.path(fpath, basename(bin))

library(RNetCDF)
