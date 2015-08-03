#!/bin/bash
mkdir build
cd build
build_dir=$PWD
inc_dir=${build_dir%%/}/include
lib_dir=${build_dir%%/}/lib
pkg_dir=${lib_dir%%/}/pkgconfig
PKG_CONFIG_PATH=${pkg_dir}:$PKG_CONFIG_PATH
export CXXFLAGS="-std=c++11 -g"
# configure and install ThirdParty
mkdir -p ThirdParty/Blas
cd ThirdParty/Blas
../../../ThirdParty/Blas/configure --prefix=$build_dir
make -j 10 install
cd ..
mkdir Lapack
cd Lapack
../../../ThirdParty/Lapack/configure --prefix=$build_dir
make -j 10 install
cd ..
mkdir HSL
cd HSL
../../../ThirdParty/HSL/configure --prefix=$build_dir
make -j 10 install
cd ..
mkdir Metis
cd Metis
../../../ThirdParty/Metis/configure --prefix=$build_dir
make -j 10 install
cd ..
mkdir Mumps
cd Mumps
../../../ThirdParty/Mumps/configure --prefix=$build_dir
make -j 10 install
cd ../..
# configure and install CoinUtils
mkdir CoinUtils
cd CoinUtils
../../CoinUtils/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Osi
mkdir Osi
cd Osi
../../Osi/configure --prefix=$build_dir
make -j 10 install
cd ..
# configure and install Clp
mkdir Clp
cd Clp
../../Clp/configure --prefix=$build_dir
make -j 10 install
cd ..
#configure and install OsiConic
mkdir OsiConic
cd OsiConic
../../OsiConic/configure --prefix=$build_dir
make -j 10 install
cd ..
#configure and install Ipopt
mkdir Ipopt
cd Ipopt
../../Ipopt/configure --prefix=$build_dir
make -j 10 install
cd ..
#configure and install Ipopt
mkdir OsiIpopt
cd OsiIpopt
../../OsiIpopt/configure --prefix=$build_dir
make -j 10 install
cd ..
### configure and install OsiMosek, uncomment if you have mosek
# mkdir OsiMosek
# cd OsiMosek
# ../../OsiMosek/configure --prefix=$build_dir
# make -j 10 install
# cd ..
### configure and install OsiCplex, uncomment if you have cplex
# mkdir OsiCplex
# cd OsiCplex
# ../../OsiCplex/configure --prefix=$build_dir
# make -j 10 install
# cd ..
#configure and install Cola
mkdir Cola
cd Cola
../../configure --prefix=$build_dir
make -j 10 install
cd ..
