L3bin
=====

R package to read L3 bin files of ocean colour data. Read all variables from any Level-3 binned file (only SeaWiFS and 
MODISA extensively tested). 

There are helper functions to generate longitude / latitude coordinates for the bins (centre or corners). These bins
are a sparse grid on Sinusoidal projection, see 

http://oceandata.sci.gsfc.nasa.gov/

http://oceancolor.gsfc.nasa.gov/SeaWiFS/TECH_REPORTS/PreLPDF/PreLVol32.pdf

Package is built using roxygen2 and Rcpp. Only tested on Linux for now - help welcome to port to Windows. 

Limitations
====
The SeaWiFS and MODIS bins of 2160 rows are assumed, extracting lonlat from CZCS does not yet work. This should be as simple as detecting the value of NUMROWS from the Bin Index. 

Basic usage
====

Get a Level-3 bin file  (this one is 15 Mb) of ocean colour, and read all bins. 

```{r}
f <- "http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/A2014208.L3b_DAY_KD490.main.bz2"
download.file(f, basename(f), mode = "wb")
library(L3bin)
system(sprintf("bunzip2 %s", basename(f)))
x <- readL3(gsub(".bz2", "", basename(f)))
List of 5
 $ bin_num: int [1:1506363] 3000338 3006072 3006080 3006081 3011818 3011827 3011828 3011829 3011830 3017567 ...
 $ nobs   : int [1:1506363] 1 1 2 4 2 6 9 6 3 1 ...
 $ nscenes: int [1:1506363] 1 1 1 1 1 1 1 1 1 1 ...
 $ weights: num [1:1506363] 1 1 1.41 2 1.41 ...
 $ Kd_490 :List of 2
  ..$ sum: num [1:1506363] 0.0414 0.0464 0.0595 0.0809 0.0564 ...
  ..$ ssq: num [1:1506363] 0.00171 0.00215 0.00251 0.00331 0.00225 ...

xy <- bin2lonlat(x$bin_num)
```

Default behaviour is to read both sum/ssq of *all* variables. You can limit to just some variables with the *vname* argument. (This can be important for the RRS files which contain many Remote Sensing Reflectance variables. 

TODO
====
- Read of multiple variables is inefficient, need to set up output vectors dynamically and only loop over bins once. 
- Extend capacity to deal with bins spatially. 
- Generalize bin definition for other sensors. 
- 
Dependencies
====

Bindings to the HDF4 library are built with Rcpp: 

You'll need at least 

- command line access, e.g. 

```{r} 
system(" hdp dumpvd")
```

- access to the HDF4 source, e.g. 
```{bash}
#include "hdf.h"
```


Build notes for HDF4
=====
I originally used these notes as a guide, building a chain of tools ultimate for use in R - but for L3bin 
you only need HDF4 in your system. 

http://scigeo.org/articles/howto-install-latest-geospatial-software-on-linux.html

```{bash}
## dependencies (on CentOS)
## sudo yum install libjpeg-devel bison flex ## ~3Mb

## HDF4
# download latest release
# check here: http://www.hdfgroup.org/ftp/HDF/HDF_Current/src/
wget http://www.hdfgroup.org/ftp/HDF/HDF_Current/src/[hdf4].tar.gz
tar xvfz [hdf4].tar.gz

## I edited the hlimits.h before compiling, see 
## https://trac.osgeo.org/gdal/wiki/HDF
##cd [hdf4]
## grep MAXFILE hdf/src/hlimits.h


# set compile options
# --disable-netcdf and --disable-fortran are necessary
# when compiling netcdf with hdf4 support
./configure \
  --prefix=/opt/source/$hdf4/build \
  --with-zlib \
  --with-jpeg \
  --enable-shared \
  --disable-netcdf \
  --disable-fortran

# compile
make 
# check build (should all pass)
make check
# install into build dir
make install
make install-examples
# check install
make installcheck

## I added this script to /etc/profile.d to get access to local/lib
##cat /etc/profile.d/custom_exports.sh
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
sudo ldconfig

```

