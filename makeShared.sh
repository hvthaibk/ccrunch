#!/bin/bash
#
#  To make binaries for gEMfitter/gEMpicker/gEMcrystal using shared libraries
#  ==========================================================================
#

numCore=$(nproc)
PCNAME=$(hostname)
dirProj=$(pwd)

if [ "$PCNAME" = "hoang-w520" ]; then
    numCore=1
fi

DIR_build64=$dirProj/build64
FILE_gem=$dirProj/bin/gEM*
FILE_log=$dirProj/test/log/makeShared.log

# bin
mkdir -p $dirProj/test/log
mkdir -p $dirProj/bin
mkdir -p $dirProj/build64

cd $dirProj/bin
rm -f $FILE_gem
cd ..

# build
rm -rf $DIR_build64
mkdir  $DIR_build64
cd     $DIR_build64
printf "Compiling shared gEM binaries...\n"
cmake .. &> $FILE_log
make -j $numCore VERBOSE=1 &>> $FILE_log
make install &>> $FILE_log
