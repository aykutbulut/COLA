#!/bin/bash

# get blas
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Blas/stable/1.4 ThirdParty/Blas
cd ThirdParty/Blas && ./get.Blas && cd ../..
# get lapack
svn co https://projects.coin-or.org/svn/BuildTools/ThirdParty/Lapack/stable/1.5 ThirdParty/Lapack
cd ThirdParty/Blas && ./get.Lapack && cd ../..
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
