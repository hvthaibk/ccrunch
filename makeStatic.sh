#!/bin/bash

# NOTE
#   - Ubuntu 12.04: this script works well
#   - Ubuntu 14.04: libtiff is 4.0 (BigTIFF) and it requires libjbig
#     in linking. But since libjbig does not yet has static lib, we
#     cannot static link gEMpicker and gEMfitter in Ubuntu 14.04
#       https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=741379

numCore=$(nproc)
PCNAME=$(hostname)
dirProj=$(pwd)

if [ "$PCNAME" = "hoang-w520" ]; then
    numCore=1
fi

CUDA_ver=$(nvcc --version | grep -oP "release \K([0-9]+)\\.([0-9]+)" | tr -dc "[:alnum:]\n\r")
OS_name=$(lsb_release -si)_$(lsb_release -sr)

DIR_nocuda=$dirProj/build64_nocuda
DIR_cuda=$dirProj/build64_cuda
FILE_gem=$dirProj/bin/gEM*

# CCRUNCH
mkdir -p $dirProj/test/log
mkdir -p $dirProj/bin
mkdir -p $dirProj/bin/exe/${OS_name}

cd $dirProj/bin
rm -f $FILE_gem
cd ..

# Without CUDA
FILE_log=$dirProj/test/log/makeStatic.log

rm -rf $DIR_nocuda
mkdir  $DIR_nocuda
cd     $DIR_nocuda
printf "Compiling static gEM binaries...\n"
cmake .. -DGEM_LINK_STATIC=ON -DGEM_FIND_CUDA=OFF &> $FILE_log
make -j $numCore VERBOSE=1 &>> $FILE_log
make install &>> $FILE_log

cd $dirProj/bin

for progName in gEMpicker gEMfitter # gEMaligner
do
    fileRoot=${progName}_nocuda
    mv ${progName} ${fileRoot}.x64
    tar -pcvzf ${fileRoot}.tar.gz ${fileRoot}.x64
    mv ${fileRoot}.tar.gz $dirProj/bin/exe/${OS_name}
done

rm -f  $FILE_gem
rm -rf $DIR_nocuda

# With CUDA
if [ "$CUDA_ver" != "" ]; then
	FILE_log=$dirProj/test/log/makeStaticCUDA.log

    rm -rf $DIR_cuda
    mkdir  $DIR_cuda
    cd     $DIR_cuda
    printf "\nCompiling CUDA-enabled static gEM binaries...\n"
    cmake .. -DGEM_LINK_STATIC=ON -DGEM_FIND_CUDA=ON &> $FILE_log
    make -j $numCore VERBOSE=1 &>> $FILE_log
    make install &>> $FILE_log

    cd $dirProj/bin

    for progName in gEMpicker gEMfitter # gEMaligner
    do
        fileRoot=${progName}_cuda_${CUDA_ver}
        mv ${progName} ${fileRoot}.x64
        tar -pcvzf ${fileRoot}.tar.gz ${fileRoot}.x64
        mv ${fileRoot}.tar.gz $dirProj/bin/exe/${OS_name}
    done

    rm -f  $FILE_gem
    rm -rf $DIR_cuda
fi
