
## Process this file with autoconf to produce a configure script.
##
## Configure.ac for L3bin
##
## Copyright (C) 2014 Michael Sumner
##
## Licensed under GNU GPL 2 or later

# The version set here will propagate to other files from here
AC_INIT([L3bin], 0.0-1)

# Checks for common programs using default macros
AC_PROG_CC

HDF_CFLAGS='/usr/local/include/'
HDF_LIBS='-L/usr/local/lib -lmfhdf -ldf -ljpeg -lz'


## Use Rscript to query Rcpp for compiler and linker flags
##
## look for Rscript, but use the one found via R_HOME to allow for multiple installations
dnl AC_DEFUN(AC_PROG_RSCRIPT, [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
dnl AC_PROG_RSCRIPT
dnl if test x"${RSCRIPT}" == x"yes" ; then
dnl      : ${R_HOME=`R RHOME`}
dnl      if test -z "${R_HOME}"; then
dnl          AC_MSG_ERROR([Could not determine R_HOME.])
dnl      fi
dnl      ## Rcpp compile flag providing header directory containing Rcpp.h
dnl      RCPP_CXXFLAGS=`${R_HOME}/bin/Rscript -e 'Rcpp:::CxxFlags()'`
dnl      ## link flag providing libary as well as path to library, and optionally rpath
dnl      RCPP_LDFLAGS=`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`
dnl else
dnl     AC_MSG_ERROR([Rscript not found, is R installed?])
dnl fi

# Now substitute these variables in src/Makevars.in to create src/Makevars
AC_SUBST(HDF_CFLAGS)
AC_SUBST(HDF_LIBS)
dnl AC_SUBST(RCPP_CXXFLAGS)
dnl AC_SUBST(RCPP_LDFLAGS)

AC_OUTPUT(src/Makevars)

