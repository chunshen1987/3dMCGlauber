#!/usr/bin/env bash

wget http://www.hepforge.org/archive/lhapdf/LHAPDF-6.2.0.tar.gz
tar -xf LHAPDF-6.2.0.tar.gz
rm -fr LHAPDF-6.2.0.tar.gz
LIBPATH=`echo $PWD/LHAPDF_Lib`
(
cd LHAPDF-6.2.0/
./configure --prefix=`echo $LIBPATH`
make
make install
)
rm -fr LHAPDF-6.2.0
