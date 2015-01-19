lib_LTLIBRARIES = libCola.la

libCola_la_SOURCES=Cone.hpp\
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

libCola_la_LIBADD=\
	$(OSICLP_LIBS)\
	$(OSICONIC_LIBS)\
	$(OSI_LIBS)\
	$(COINUTILS_LIBS)\
	$(LIBS)

bin_PROGRAMS = cola

cola_SOURCES=Main.cpp

cola_LDADD=\
	libCola.la

AM_CPPFLAGS=\
	-I../../include\
	$(OSICLP_CFLAGS)\
	$(OSICONIC_CFLAGS)\
	$(OSI_CFLAGS)\
	$(CLP_CFLAGS)\
	$(COINUTILS_CFLAGS)


# header files to be installed to $(includedir)
include_HEADERS=Cone.hpp\
	LorentzCone.hpp\
	ScaledCone.hpp\
	ColaModel.hpp\
	Separate.hpp\
	Options.hpp

# this are in the disctribution but not in include
#noinst_HEADERS=
