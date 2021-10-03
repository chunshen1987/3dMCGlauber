#!/usr/bin/env bash

num_of_cores=$1
if [ -z "$num_of_cores" ]
then
    num_of_cores=4
fi

version="6.2.1"
echo "Building LHAPDF ${version} ... "
tar -xf utilities/LHAPDF-${version}.tar.gz
LIBPATH=`echo $PWD/LHAPDF_Lib`
(
    cd LHAPDF-${version}/
    ./configure --prefix=`echo $LIBPATH` --disable-python
    make -j${num_of_cores}
    make install
)
rm -fr LHAPDF-${version}

echo "downloading pdfsets ... "
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/CT10nnlo.tar.gz -O- | tar xz -C $PWD/LHAPDF_Lib/share/LHAPDF
