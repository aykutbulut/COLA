AC_PREREQ([2.67])
AC_INIT([cola], [0.0], [aykut@lehigh.edu])
AC_CONFIG_SRCDIR([config.h.in])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_HEADERS([config.h])

# for creating library
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX

# Installs to ./bin, changes default behavior
AC_PREFIX_DEFAULT([`pwd`])

# Checks for header files.
AC_CHECK_HEADERS([math.h])

# Checks for typedefs, str
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Check dependent libraries
AC_CHECK_LIB([m], [sqrt])
# Check dependent libraries through pkg-config
PKG_CHECK_MODULES([OSICLP], [osi-clp >= 0.105])
PKG_CHECK_MODULES([OSICONIC], [osiconic >= 0.0])
PKG_CHECK_MODULES([OSI], [osi >= 0.105])
PKG_CHECK_MODULES([COINUTILS], [coinutils >= 2.8])

AC_CONFIG_FILES([Makefile
	src/Makefile
	cola.pc
	cola-uninstalled.pc])
AC_OUTPUT
