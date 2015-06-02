# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author: Aykut Bulut, Lehigh University, 05/22/2015

#############################################################################
#                       Names and other basic things                        #
#############################################################################
AC_PREREQ(2.59)
AC_INIT([Cola],[0.1],[https://github.com/aykutbulut/COLA/issues/new])

AC_COPYRIGHT([
/*===========================================================================*/
/*                                                                           */
/* Cola is developed by Aykut Bulut (aykutblt@gmail.com)                     */
/*                                                                           */
/* (c) Copyright 2015 Lehigh University. All Rights Reserved.                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/ColaModel.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################
# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(Cola)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Cola)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
if test "x$F77" != xunavailable
then
  AC_COIN_F77_WRAPPERS
fi

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN-OR components                           #
#############################################################################
AC_COIN_CHECK_PACKAGE(CoinUtils, [coinutils], [ColaLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHECK_PACKAGE(Osi,       [osi],       [ColaLib])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not available.])
fi
AC_COIN_CHECK_PACKAGE(OsiClp,    [osi-clp],   [ColaLib])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package OsiClp not available.])
fi
AC_COIN_CHECK_PACKAGE(OsiConic,  [osiconic],  [ColaLib])
if test $coin_has_osiconic != yes ; then
  AC_MSG_ERROR([Required package OsiConic not available.])
fi

AC_COIN_CHECK_PACKAGE_BLAS(ColaLib)
AC_COIN_CHECK_PACKAGE_LAPACK(ColaLib)
AC_COIN_CHECK_CXX_CHEADER(math)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################
AC_COIN_DOXYGEN(CoinUtils Osi OsiConic)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################
# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 cola.pc
                 cola-uninstalled.pc])

AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cola.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
