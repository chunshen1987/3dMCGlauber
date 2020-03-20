#!/usr/bin/env bash

version="6.2.1"
wget http://www.hepforge.org/archive/lhapdf/LHAPDF-$version.tar.gz
tar -xf LHAPDF-$version.tar.gz
rm -fr LHAPDF-$version.tar.gz
LIBPATH=`echo $PWD/LHAPDF_Lib`
(
cd LHAPDF-$version/
./configure --prefix=`echo $LIBPATH` --disable-python
make
make install
)
rm -fr LHAPDF-$version

echo "downloading pdfsets ... "
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT10nnlo.tar.gz -O- | tar xz -C $PWD/LHAPDF_Lib/share/LHAPDF
