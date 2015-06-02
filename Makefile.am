# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University 2015-05-22

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

# subdirs is set by configure as the list of all subdirectories to recurse
# into
SUBDIRS = src

########################################################################
#                    Creation of the addlibs file                      #
########################################################################

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = cola.pc

addlibsdir = $(DESTDIR)$(datadir)/coin/doc/Cola

install-data-hook:
	@$(mkdir_p) "$(addlibsdir)"
if COIN_HAS_PKGCONFIG
	PKG_CONFIG_PATH=@COIN_PKG_CONFIG_PATH@ \
	$(PKG_CONFIG) --libs cola > $(addlibsdir)/cola_addlibs.txt
else
if COIN_CXX_IS_CL
	echo "-libpath:`$(CYGPATH_W) @abs_lib_dir@` libCola.lib @COLALIB_LIBS_INSTALLED@" > $(addlibsdir)/cola_addlibs.txt
else
	echo -L@abs_lib_dir@ -lCola @COLALIB_LIBS_INSTALLED@ > $(addlibsdir)/cola_addlibs.txt
endif
endif

uninstall-hook:
	rm -f $(addlibsdir)/cola_addlibs.txt

########################################################################
#                         Maintainer Stuff                             #
########################################################################

CLEANFILES =
# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES =
EXTRA_DIST =
include BuildTools/Makemain.inc
