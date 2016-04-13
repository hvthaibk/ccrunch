#!/bin/bash

# lastest boost: https://launchpad.net/~boost-latest/+archive/ubuntu/ppa

dirCur=$(pwd)
numCore=$(nproc)

#COMPILE_VTK=yes
#COMPILE_ITK=yes
COMPILE_OPENCV=yes
COMPILE_TIFF=yes
COMPILE_GPP4=yes
COMPILE_TCLAP=yes
COMPILE_ARMADILLO=yes

function build_cmake()
{
  srcdir="$1"
  args="$2"
  builddir="$srcdir/build"

  echo 'Compiling '$(basename "$srcdir")'...'

  mkdir -p "$builddir"
  pushd "$builddir"
  cmake -D CMAKE_BUILD_TYPE=RELEASE -D BUILD_SHARED_LIBS=OFF -D BUILD_TESTING=OFF -D BUILD_EXAMPLES=OFF $args -D CMAKE_INSTALL_PREFIX=$dirCur/3rdParty/install64 ..
  make -j $numCore install
  cmake -D CMAKE_BUILD_TYPE=RELEASE -D BUILD_SHARED_LIBS=ON -D BUILD_TESTING=OFF -D BUILD_EXAMPLES=OFF $args -D CMAKE_INSTALL_PREFIX=$dirCur/3rdParty/install64 ..
  make -j $numCore install
  popd
  rm -r "$builddir"
}

echo 'We assume that the following dependencies are available:'
echo '
  - g++
  - cmake
  - pkg-config
  - libgtk2.0-dev (for OpenCV windows)
  - automake/autoconf (for tclap from git)
  - python-dev (for python wrapper)
  - python-numpy (for OpenCV python wrapper)
  - python-matplotlib

  - libgsl0-dev
  - libtiff4-dev
  - libfftw3-dev
  - libboost1.55-all-dev or
  -     libboost-thread1.55-dev
  -     libboost-system1.55-dev
  -     libboost-filesystem1.55-dev
  -     libboost-random1.55-dev
  -     libboost-python1.55-dev
  -     libboost-math1.55-dev (for multiprecision)

  - libblas-dev
  - liblapack-dev
  - libarpack2-dev
  - libhdf5-serial-dev
'
echo 'They are all often available with your distribution package manager.'

echo 'Other dependencies are handled directly by this script or with git submodules.'

echo -n 'Finished reading? Then type any key...' ; read -n 1 ; echo

echo 'Updating git submodules...'
# esbtl / tclap / itk / vtk / opencv

git submodule init
# checkout | rebase | merge | none
git config submodule.3rdParty/itk.update none
git config submodule.3rdParty/vtk.update none
git submodule update

[ -n "$COMPILE_VTK" -a ! -d $dirCur/3rdParty/install64/include/vtk-* ] && \
  build_cmake 3rdParty/vtk "-D VTK_WRAP_JAVA=OFF -D VTK_WRAP_PYTHON=ON -D VTK_WRAP_TCL=OFF -D"
[ -n "$COMPILE_ITK" -a ! -d $dirCur/3rdParty/install64/include/ITK-* ] && \
  build_cmake 3rdParty/itk "-D ITK_BUILD_DEFAULT_MODULES=ON -D USE_WRAP_ITK=ON -D WRAP_ITK_PYTHON=ON -D VTK_DIR=$dirCur/3rdParty/install64/lib/cmake/vtk-6.2 -D Module_ITKVtkGlue=ON"
[ -n "$COMPILE_OPENCV" -a ! -d $dirCur/3rdParty/install64/include/opencv ] && \
  build_cmake 3rdParty/opencv "-D BUILD_TESTS=OFF -D WITH_CUDA=OFF -D WITH_OPENCL=OFF -D BUILD_TIFF=ON"

# Armadillo 6.200.2
if [ -n "$COMPILE_ARMADILLO" -a ! -d $dirCur/3rdParty/install64/include/armadillo_bits ]; then
  libName=armadillo-6.200.2
  echo 'Compiling Armadillo...'
  cd $dirCur/3rdParty
  [ -r $libName.tar.gz ] || wget http://sourceforge.net/projects/arma/files/$libName.tar.gz
  tar xf $libName.tar.gz
  build_cmake $libName ""
  rm -rf $libName
  cd $dirCur
fi

# TIFF 4.0.3
if [ -n "$COMPILE_TIFF" -a ! -f $dirCur/3rdParty/install64/include/tiff.h ]; then
  libName=tiff-4.0.3
  echo 'Compiling' ${libName} '...'
  cd $dirCur/3rdParty
  [ -r $libName.tar.gz ] || wget ftp://ftp.remotesensing.org/pub/libtiff/$libName.tar.gz
  tar xf $libName.tar.gz
  cd $dirCur/3rdParty/$libName
  $dirCur/3rdParty/$libName/configure --prefix=$dirCur/3rdParty/install64/ --enable-static=yes --enable-shared=yes --disable-jbig --disable-lzma
  make -j $numCore
  make install
  cd $dirCur/3rdParty
  rm -rf $dirCur/3rdParty/$libName
  cd $dirCur
fi

# GPP4 1.3.1
if [ -n "$COMPILE_GPP4" -a ! -d $dirCur/3rdParty/install64/include/gpp4 ]; then
  libName=gpp4-1.3.1
  echo 'Compiling' ${libName} '...'
  cd $dirCur/3rdParty
  [ -r ${libName}.tar.gz ] || wget https://launchpad.net/gpp4/1.3/1.3.1/+download/${libName}.tar.gz
  tar xf ${libName}.tar.gz
  cd $dirCur/3rdParty/$libName
  export MMDB_CFLAGS="nimp"
  export MMDB_LIBS="nimp"
  $dirCur/3rdParty/$libName/configure --without-fortran-api --prefix=$dirCur/3rdParty/install64
  make -j $numCore
  make install
  cd $dirCur/3rdParty
  rm -rf $dirCur/3rdParty/$libName
  cd $dirCur
fi

# tclap 1.2.1
if [ -n "$COMPILE_TCLAP" -a ! -d $dirCur/3rdParty/install64/include/tclap ]; then
  echo 'Compiling tclap...'
  cd $dirCur/3rdParty/tclap
  ./autotools.sh
  ./configure --prefix=$dirCur/3rdParty/install64/
  make -j $numCore
  # workaround BS bug without doxygen
  mkdir -p docs/html
  make install
  cd $dirCur
fi

# NR300 / ESBTL / Eigen
cd $dirCur/3rdParty
cp -r NR300 ./install64/include

cd $dirCur/3rdParty
mkdir -p ./install64/include/ESBTL
cp -r ./esbtl/include/ESBTL/ ./install64/include/

cd $dirCur/3rdParty
eigen=eigen-3.2.2.tar.bz2
[ -r $eigen ] || wget http://bitbucket.org/eigen/eigen/get/3.2.2.tar.bz2 -O $eigen
tar xf $eigen
mkdir -p ./install64/include/Eigen
cp -r ./eigen-eigen-*/Eigen/ ./install64/include/
rm -r ./eigen-eigen-*

cd $dirCur

# test directory
cd $dirCur
mkdir -p test
[ -r test/CMakeLists.txt ] || touch test/CMakeLists.txt

# CCRUNCH - development
echo 'Compiling ccrunch...'
cd $dirCur
mkdir -p bin
mkdir -p build64
cd build64
[ -r Makefile ] || cmake -D BUILD_TESTS=OFF ..
make -j $numCore
make install
#make clean
