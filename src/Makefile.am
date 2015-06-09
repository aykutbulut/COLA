# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University 2015-05-22

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCola                                #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCola.la

# List all source files for this library, including headers
libCola_la_SOURCES=\
	config_cola.h\
	ColaConfig.h\
	Cone.hpp\
	Cone.cpp\
	LorentzCone.hpp\
	LorentzCone.cpp\
	ScaledCone.hpp\
	ScaledCone.cpp\
	Separate.hpp\
	Separate.cpp\
	Options.hpp\
	Options.cpp\
	ColaModel.hpp\
	ColaModel.cpp

# List all additionally required libraries
if DEPENDENCY_LINKING
libCola_la_LIBADD = $(COLALIB_LIBS)
endif

# This is for libtool
libCola_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags.
AM_CPPFLAGS = $(COINUTILS_CFLAGS)\
	 $(OSI_CFLAGS)\
	 $(OSICONIC_CFLAGS)\
	 $(OSICLP_CFLAGS)

########################################################################
#                            cola program                              #
########################################################################

# Name of the executable compiled in this directory.  We want it to be
# installed in the 'bin' directory
bin_PROGRAMS = cola

# List all source files for this executable, including headers
cola_SOURCES = Main.cpp

# List all additionally required libraries
cola_LDADD = libCola.la $(COLALIB_LIBS)

# List all dependency libraries (similar to LDADD, but without -l, -L flags)
cola_DEPENDENCIES = libCola.la $(COLALIB_DEPENDENCIES)

########################################################################
#                            Additional flags                          #
########################################################################

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	Cone.hpp\
	LorentzCone.hpp\
	ScaledCone.hpp\
	ColaModel.hpp\
	Separate.hpp\
	Options.hpp

install-exec-local:
	$(install_sh_DATA) config_cola.h $(DESTDIR)$(includecoindir)/ColaConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/ColaConfig.h
