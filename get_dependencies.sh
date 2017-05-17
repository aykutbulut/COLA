#!/bin/bash

# get blas
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Blas/stable/1.4 ThirdParty/Blas
cd ThirdParty/Blas && ./get.Blas && cd ../..
# get lapack
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Lapack/stable/1.5 ThirdParty/Lapack
cd ThirdParty/Lapack && ./get.Lapack && cd ../..
# coin glpk
# svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Glpk/stable/1.10 ThirdParty/Glpk
# get coin-or coinutils
svn co https://projects.coin-or.org/svn/CoinUtils/stable/2.10/CoinUtils
# get coin-or osi
svn co https://projects.coin-or.org/svn/Osi/stable/0.107/Osi
# get coin-or clp
svn co https://projects.coin-or.org/svn/Clp/stable/1.16/Clp
# get osiconic
git clone https://github.com/aykutbulut/OSI-CONIC OsiConic
# get Ipopt
svn co https://projects.coin-or.org/svn/Ipopt/releases/3.12.3/Ipopt Ipopt
# IPOPT dependencies
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/HSL/stable/1.5    ThirdParty/HSL
##cd ThirdParty/HSL && ./get.HSL && cd ../..
wait
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Metis/stable/1.3  ThirdParty/Metis
cd ThirdParty/Metis && ./get.Metis && cd ../..
wait
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Mumps/stable/1.5  ThirdParty/Mumps
cd ThirdParty/Mumps && ./get.Mumps && cd ../..
wait
# get OsiIpopt
git clone https://github.com/aykutbulut/OsiIpopt OsiIpopt
# get coin-or OsiMosek, uncomment if you have mosek, optional
#git clone https://github.com/aykutbulut/OSI-MOSEK OsiMosek
# get coin-or OsiCplex, uncomment if you have cplex, optional
#git clone https://github.com/aykutbulut/OsiCplex OsiCplex
