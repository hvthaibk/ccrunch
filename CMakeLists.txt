cmake_minimum_required(VERSION 2.8.2)

project(ccrunch)

#-----------------------------------------------------------------------
# Version number
set(ccrunch_VERSION_MAJOR 0)
set(ccrunch_VERSION_MINOR 1)
set(ccrunch_VERSION_PATCH 0)

#-----------------------------------------------------------------------
# Shared libs directory
set(3RDPARTY_ROOT_DIR       ${PROJECT_SOURCE_DIR}/3rdParty/install64/)
set(3RDPARTY_INCLUDE_DIR    ${PROJECT_SOURCE_DIR}/3rdParty/install64/include/)
set(3RDPARTY_LIB_DIR        ${PROJECT_SOURCE_DIR}/3rdParty/install64/lib/)
set(ENV{PKG_CONFIG_PATH}    ${3RDPARTY_LIB_DIR}/pkgconfig)

include_directories(${3RDPARTY_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/3rdParty/)

link_directories(${3RDPARTY_LIB_DIR})
include(CheckLibraryExists)
include(FindPkgConfig)

set(PROJECT_INCLUDE_DIR     ${PROJECT_SOURCE_DIR}/bin/include)
set(PROJECT_LIB_DIR         ${PROJECT_SOURCE_DIR}/bin/lib)
set(GEMLIB_INSTALL_RPATH    "$ORIGIN/../../3rdParty/install64/lib/")
set(GEMMEX_INSTALL_RPATH    "$ORIGIN/../../3rdParty/install64/lib/;$ORIGIN/../lib/")
set(GEMEXE_INSTALL_RPATH    "$ORIGIN/../3rdParty/install64/lib/;$ORIGIN/lib/")
set(CMAKE_INSTALL_PREFIX    ${PROJECT_SOURCE_DIR})

#-----------------------------------------------------------------------
# GEM directives
option(BUILD_TESTS          "CTest using CDash" OFF)

option(GEM_LINK_STATIC      "Static linking"    OFF)

option(GEM_FIND_HDF5        "Find HDF5"         ON)
option(GEM_FIND_MKL         "Find MKL"          OFF)
option(GEM_FIND_CUDA        "Find CUDA"         ON)
option(GEM_FIND_MPI         "Find MPI"          OFF)
option(GEM_FIND_OpenMP      "Find OpenMP"       ON)
option(GEM_FIND_OpenCV      "Find OpenCV"       ON)
option(GEM_FIND_GSL         "Find GSL"          ON)
option(GEM_FIND_QT4         "Find QT4"          ON)
option(GEM_FIND_Armadillo   "Find Armadillo"    ON)
option(GEM_FIND_HEX         "Find HEX"          OFF)
option(GEM_FIND_PYTHON      "Find Python"       ON)

if(GEM_LINK_STATIC)
    set(GEM_FIND_HDF5        ON)
    set(GEM_FIND_MKL         OFF)
    set(GEM_FIND_MPI         OFF)
    set(GEM_FIND_OpenMP      OFF)
    set(GEM_FIND_OpenCV      OFF)
    set(GEM_FIND_GSL         ON)
    set(GEM_FIND_QT4         OFF)
    set(GEM_FIND_Armadillo   ON)
    set(GEM_FIND_HEX         OFF)
    set(GEM_FIND_PYTHON      OFF)

    option(__GEM_LIMIT_FUNC__   "Limit functionality"   ON)

    message(STATUS "Linking against static libraries")
else(GEM_LINK_STATIC)
    message(STATUS "Linking against dynamic libraries")
endif(GEM_LINK_STATIC)

if(GEM_FIND_VTK AND GEM_FIND_CUDA)
    message(FATAL_ERROR "ERROR: cannot compile CUDA code with VTK header files")
endif()

#-----------------------------------------------------------------------
# Macros
macro(gem_find_dir return_list root_path)
    file(GLOB_RECURSE new_list ${root_path})
    set(dir_list "")
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()
    list(REMOVE_DUPLICATES dir_list)
    set(${return_list} ${dir_list})
endmacro(gem_find_dir)

macro(gem_use_notuse key str)
    if(${key})
        set(MSG_USE    ${MSG_USE}    "${str}, ")
    else()
        set(MSG_NOTUSE ${MSG_NOTUSE} "${str}, ")
    endif()
endmacro(gem_use_notuse)

macro(gem_install_header root_path)
    file(GLOB_RECURSE new_list ${root_path})
    foreach(file_path ${new_list})
        install(FILES      ${file_path}    DESTINATION ${PROJECT_INCLUDE_DIR})
    endforeach()
endmacro(gem_install_header)

macro(gem_find_source return_list root_path)
    file(GLOB_RECURSE new_list ${root_path})
    set(file_list "")
    foreach(file_path ${new_list})
        set(file_list ${file_list} ${file_path})
    endforeach()
    set(${return_list} ${file_list})
endmacro(gem_find_source)

#-----------------------------------------------------------------------
# Python
if(GEM_FIND_PYTHON)
    find_package(PythonLibs REQUIRED)
    if(PYTHONLIBS_FOUND)
        message(STATUS "Found Python version " ${PYTHONLIBS_VERSION_STRING})
        include_directories("${PYTHON_INCLUDE_DIRS}")
    endif()
endif()

#-----------------------------------------------------------------------
# Boost
set(BOOST_ROOT ${3RDPARTY_ROOT_DIR})
set(Boost_ADDITIONAL_VERSIONS "1.53" "1.53.0")
set(Boost_USE_MULTITHREADED ON)
if(GEM_LINK_STATIC)
    set(Boost_USE_STATIC_LIBS ON)
else()
    set(Boost_USE_STATIC_LIBS OFF)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_RUNTIME OFF)
endif()

find_package(Boost 1.53.0 REQUIRED COMPONENTS thread system filesystem random python math_c99 math_tr1)
if(Boost_FOUND)
    message(STATUS "Found Boost")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    if(GEM_LINK_STATIC)
        set(BOOST_LIBRARY boost_thread.a boost_system.a boost_filesystem.a boost_random.a boost_python.a boost_math_c99.a boost_math_tr1.a pthread rt)
        set(EXTRA_LIBS ${EXTRA_LIBS} ${BOOST_LIBRARY})
    else(GEM_LINK_STATIC)
        set(BOOST_LIBRARY boost_thread boost_system boost_filesystem boost_random boost_python boost_math_c99 boost_math_tr1)
        set(EXTRA_LIBS ${EXTRA_LIBS} ${BOOST_LIBRARY})
    endif(GEM_LINK_STATIC)
endif()

#-----------------------------------------------------------------------
# Intel MKL & FFTW 3.3.2
if(GEM_FIND_MKL)
    include(CMakeForceCompiler)
    CMAKE_FORCE_C_COMPILER(icc "Intel C Compiler")
    CMAKE_FORCE_CXX_COMPILER(icpc "Intel C++ Compiler")
    set(MKL_ROOT "/opt/intel/mkl")
    include_directories(${MKL_ROOT}/include/fftw/)
    link_directories(${MKL_ROOT}/lib/intel64)
    if(GEM_LINK_STATIC)
        set(MKL_LIBS mkl_cdft_core.a mkl_intel_ilp64.a mkl_intel_thread.a mkl_core.a mkl_blacs_openmpi_ilp64.a mkl_gnu_thread.a pthread m)
    else(GEM_LINK_STATIC)
        set(MKL_LIBS mkl_intel_ilp64 mkl_core mkl_gnu_thread pthread m)
    endif(GEM_LINK_STATIC)
    set(EXTRA_LIBS ${EXTRA_LIBS} ${MKL_LIBS})
else()
    check_library_exists("fftw3"    fftw_execute     ${3RDPARTY_LIB_DIR} HAVE_FFTW)
    check_library_exists("fftw3f"   fftwf_execute    ${3RDPARTY_LIB_DIR} HAVE_FFTWF)
    if(HAVE_FFTW AND HAVE_FFTWF)
        message(STATUS "Found FFTW")
        if(GEM_LINK_STATIC)
            set(FFTW_LIBRARY fftw3.a fftw3f.a fftw3_threads.a fftw3f_threads.a)
        else(GEM_LINK_STATIC)
            set(FFTW_LIBRARY fftw3 fftw3f fftw3_threads fftw3f_threads)
        endif(GEM_LINK_STATIC)
        set(EXTRA_LIBS ${EXTRA_LIBS} ${FFTW_LIBRARY})
    else()
        message(FATAL_ERROR "ERROR: unable to find FFTW library")
    endif()
endif()

#-----------------------------------------------------------------------
# MPI and OpenMP
if(GEM_FIND_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        message(STATUS "Found MPI")
        option(__GEM_USE_MPI__ "Use MPI" ON)
        include_directories(${MPI_INCLUDE_PATH})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY})
    endif()
endif()

if(GEM_FIND_OpenMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        message(STATUS "Found OpenMP")
        option(__GEM_USE_OPENMP__ "Use OpenMP" ON)
        set(CMAKE_C_FLAGS           ${CMAKE_C_FLAGS}           ${OpenMP_C_FLAGS})
        set(CMAKE_CXX_FLAGS         ${CMAKE_CXX_FLAGS}         ${OpenMP_CXX_FLAGS})
        set(CMAKE_EXE_LINKER_FLAGS  ${CMAKE_EXE_LINKER_FLAGS}  ${OpenMP_EXE_LINKER_FLAGS})
    endif()
endif()

#-----------------------------------------------------------------------
# CUDA (at least 5.5)
if(GEM_FIND_CUDA)
    find_package(CUDA 5.5 REQUIRED)
endif()
if(CUDA_FOUND)
    message(STATUS "Found CUDA")
    include_directories(${CUDA_INCLUDE_DIRS})
    message(STATUS "  toolkit version: " ${CUDA_VERSION_MAJOR} "." ${CUDA_VERSION_MINOR})

    if(__GEM_USE_STD11__)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)  # to avoid -std=c++0x
    else()
        set(CUDA_PROPAGATE_HOST_FLAGS ON)
        string(REPLACE " -ftemplate-depth=50" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        string(REPLACE " -ftemplate-depth-50" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
        string(REPLACE " -ftemplate-depth=50" "" CUDA_HOST_FLAGS "${CUDA_HOST_FLAGS}")
        string(REPLACE " -ftemplate-depth-50" "" CUDA_HOST_FLAGS "${CUDA_HOST_FLAGS}")
    endif()

    # Determine CUDA compute capability
    try_run(CUDA_QUERY_RUN CUDA_QUERY_COMPILE
            ${CMAKE_BINARY_DIR}
            ${PROJECT_SOURCE_DIR}/programs/infoCUDA.c
            CMAKE_FLAGS
            -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
            -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
            COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
            RUN_OUTPUT_VARIABLE     RUN_OUTPUT_VAR)

    # CUDA_QUERY_COMPILE is TRUE when compile succeeds
    # CUDA_QUERY_RUN is zero when a GPU is found
    if(CUDA_QUERY_COMPILE AND NOT CUDA_QUERY_RUN)
        set(CUDA_HAVE_GPU TRUE CACHE BOOL "Whether CUDA-capable GPU is present")
        set(CUDA_COMPUTE_CAPABILITY ${RUN_OUTPUT_VAR} CACHE STRING "Compute capability of CUDA-capable GPU present")
        set(CUDA_GENERATE_CODE "-gencode arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}" CACHE STRING "Which GPU architectures to generate code for (each arch/code pair will be passed as --generate-code option to nvcc, separate multiple pairs by ;)")
        mark_as_advanced(CUDA_COMPUTE_CAPABILITY CUDA_GENERATE_CODE)

        # Compilation directives
        if((${CUDA_VERSION_MAJOR} LESS 4) OR ((${CUDA_VERSION_MAJOR} EQUAL 4) AND (${CUDA_VERSION_MINOR} EQUAL 0)))
            option(__GEM_CUDA_NVCC_VERSION_40__ "CUDA version 4.0" ON)
        endif()
        if(${CUDA_COMPUTE_CAPABILITY} LESS 20)
            #option(__GEM_CUDA_ARCH_HOST_130__ "CUDA architecture 1.3" ON)
        endif()

        # Set nvcc flags
        set(CUDA_ARCH_FLAGS "")
        if(GEM_LINK_STATIC)
            if(${CUDA_VERSION_MAJOR} EQUAL 4)
                set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_20")
                set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_21")
            endif()
            if(${CUDA_VERSION_MAJOR} EQUAL 5)
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_20")
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_21")
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_30,code=sm_30")
                set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_35,code=sm_35")
            endif()
        else()
            if(${CUDA_COMPUTE_CAPABILITY} EQUAL 20 OR ${CUDA_COMPUTE_CAPABILITY} EQUAL 21)
                set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_21")
            endif()
            if(${CUDA_COMPUTE_CAPABILITY} EQUAL 35)
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_20")
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_20,code=sm_21")
                #set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_30,code=sm_30")
                set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS} " -gencode arch=compute_35,code=sm_35")
            endif()

            # Keep this type of flag setting since -gencode does not work on IGBMC's workstation, adonis cluster
            #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch;compute_13;-code;sm_13")
            #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch;compute_20;-code;sm_20,sm_21")
            #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch;compute_20;-code;sm_20,sm_30")
            #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} "-arch;compute_35;-code;sm_35")
        endif()

        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${CUDA_ARCH_FLAGS})
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -Xcudafe --diag_suppress=unrecognized_gcc_pragma")
        set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} " -lineinfo")                 # for Visual Profiler

        message(STATUS "  hardware compute capability: " ${CUDA_COMPUTE_CAPABILITY})
        message(STATUS "  compiler option: " ${CUDA_NVCC_FLAGS})

        option(__GEM_USE_CUDA__   "Use CUDA"   ON)
        option(__GEM_USE_THRUST__ "Use Thrust" ON)
    else()
        set(CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present")
        message(STATUS "  ----------------------------------------")
        message(STATUS "  cannot determine CUDA_COMPUTE_CAPABILITY")
        message(STATUS "  ----------------------------------------")
        message(STATUS "  CUDA_QUERY_COMPILE = " ${CUDA_QUERY_COMPILE})
        message(STATUS "  CUDA_QUERY_RUN     = " ${CUDA_QUERY_RUN})

        if(NOT CUDA_QUERY_COMPILE)
            message(STATUS "  cannot compile CUDA query code")
        endif()

        if(CUDA_QUERY_RUN)
            message(STATUS "  cannot detect CUDA-capable gpus")
        endif()

        set(CUDA_FOUND FALSE)
        option(__GEM_USE_CUDA__   "Use CUDA"   OFF)
        option(__GEM_USE_THRUST__ "Use Thrust" OFF)
    endif()
endif()

#-----------------------------------------------------------------------
# GPP4 1.3.1 & TCLAP & TIFF 3.9.6 & HDF5
check_library_exists("gpp4"     ccp4_cmap_open   ${3RDPARTY_LIB_DIR} HAVE_GPP4)
if(HAVE_GPP4)
    message(STATUS "Found GPP4")
    option(__GEM_USE_GPP4__ "Use GPP4" ON)
    if(GEM_LINK_STATIC)
        set(GPP4_LIBRARY gpp4.a)
    else(GEM_LINK_STATIC)
        set(GPP4_LIBRARY gpp4)
    endif(GEM_LINK_STATIC)
    set(EXTRA_LIBS ${EXTRA_LIBS} ${GPP4_LIBRARY})
else()
    message(FATAL_ERROR "ERROR: unable to find GPP4 library")
endif()

check_library_exists("tiff"     TIFFOpen         ${3RDPARTY_LIB_DIR} HAVE_TIFF)
if(HAVE_TIFF)
    message(STATUS "Found TIFF")
    if(GEM_LINK_STATIC)
        set(TIFF_LIBRARY tiff.a jpeg.a z.a)
    else(GEM_LINK_STATIC)
        set(TIFF_LIBRARY tiff)
    endif(GEM_LINK_STATIC)
    set(EXTRA_LIBS ${EXTRA_LIBS} ${TIFF_LIBRARY})
else()
    message(FATAL_ERROR "ERROR: unable to find TIFF library")
endif()

pkg_check_modules(TCLAP REQUIRED tclap)
if(TCLAP_FOUND)
    if(GEM_LINK_STATIC)
        set(EXTRA_LIBS ${EXTRA_LIBS} ${TCLAP_STATIC_LIBRARIES})
    else()
        set(EXTRA_LIBS ${EXTRA_LIBS} ${TCLAP_LIBRARIES})
    endif()
endif()

if(GEM_FIND_HDF5)
    if(GEM_LINK_STATIC)
        set(HDF5_USE_STATIC_LIBRARIES ON)
    endif()
    find_package(HDF5)
    if(HDF5_FOUND)
        message(STATUS "Found HDF5")
        option(__GEM_USE_HDF5__ "Use HDF5" ON)
        include_directories(${HDF5_INCLUDE_DIRS})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${HDF5_LIBRARIES})
    endif()
endif()

#-----------------------------------------------------------------------
# Qt4
if(GEM_FIND_QT4)
    find_package(Qt4 4.6)
    if(QT4_FOUND)
        message(STATUS "Found Qt4 " ${QT_VERSION})
        option(__GEM_USE_QT4__ "Use Qt4" ON)
        include(${QT_USE_FILE})
        include_directories(${QT_INCLUDES})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${QT_LIBRARIES})
    endif()
endif()

#-----------------------------------------------------------------------
# OpenCV and GSL
if(GEM_FIND_OpenCV)
    #set(OpenCV_DIR "/usr/share/OpenCV")
    set(OpenCV_DIR "${3RDPARTY_ROOT_DIR}/share/OpenCV")
    find_package(OpenCV)
    if(OpenCV_FOUND)
        message(STATUS "Found OpenCV " ${OpenCV_VERSION})
        option(__GEM_USE_OPENCV__ "Use OpenCV" ON)
        include_directories(${OpenCV_INCLUDE_DIRS})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${OpenCV_LIBS})
    endif()
endif()

if(GEM_FIND_GSL)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})
    find_package(GSL REQUIRED)
    if(GSL_FOUND)
        message(STATUS "Found GSL")
        option(__GEM_USE_GSL__ "Use GSL" ON)
        include_directories(${GSL_INCLUDE_DIRS})
        set(EXTRA_LIBS ${EXTRA_LIBS} ${GSL_LIBRARIES})
    endif()
endif()

#-----------------------------------------------------------------------
# Armadillo
if(GEM_FIND_Armadillo)
    set(ARMADILLO_INCLUDE_DIR ${3RDPARTY_INCLUDE_DIR})
    set(ARMADILLO_LIBRARY     ${3RDPARTY_LIB_DIR})

    if(GEM_LINK_STATIC)
        set(BLA_STATIC ON)
    endif()
    find_package(BLAS)
    if(BLAS_FOUND)
        message(STATUS "  Found BLAS")
        set(EXTRA_LIBS ${EXTRA_LIBS} ${BLAS_LIBRARIES} ${BLAS95_LIBRARIES})
    endif(BLAS_FOUND)

    find_package(LAPACK)
    if(LAPACK_FOUND)
        message(STATUS "  Found LAPACK")
        set(EXTRA_LIBS ${EXTRA_LIBS} ${LAPACK_LIBRARIES} ${LAPACK95_LIBRARIES})
    endif(LAPACK_FOUND)

    find_package(Armadillo)
    if(ARMADILLO_FOUND)
        message(STATUS "Found Armadillo " ${ARMADILLO_VERSION_STRING})
        option(__GEM_USE_ARMADILLO__ "Use Armadillo" ON)
        include_directories(${ARMADILLO_INCLUDE_DIRS})
        find_library(ARMADILLO_LIBS
                     NAMES armadillo
                     PATHS ${3RDPARTY_LIB_DIR})
        if(GEM_LINK_STATIC)
#            set(ARMADILLO_LIBS armadillo.a)
        else(GEM_LINK_STATIC)
#            set(ARMADILLO_LIBS armadillo   lapack  blas   arpack)
        endif(GEM_LINK_STATIC)
        set(EXTRA_LIBS ${EXTRA_LIBS} ${ARMADILLO_LIBS})
    endif(ARMADILLO_FOUND)
endif()

#-----------------------------------------------------------------------
# Used libraries
set(MSG_USE     "Use libraries: ")
set(MSG_NOTUSE  "Do not use libraries: ")

gem_use_notuse(__GEM_USE_STD11__        "C++11")
gem_use_notuse(__GEM_USE_HDF5__         "HDF5")
gem_use_notuse(__GEM_USE_QT4__          "Qt4")
gem_use_notuse(__GEM_USE_OPENCV__       "OpenCV")
gem_use_notuse(__GEM_USE_GSL__          "GSL")
gem_use_notuse(__GEM_USE_MPI__          "MPI")
gem_use_notuse(__GEM_USE_OPENMP__       "OpenMP")
gem_use_notuse(__GEM_USE_CUDA__         "CUDA")
gem_use_notuse(__GEM_USE_ARMADILLO__    "Armadillo")
gem_use_notuse(__GEM_USE_HEX__          "HEX")

message(STATUS ${MSG_USE})
message(STATUS ${MSG_NOTUSE})

#-----------------------------------------------------------------------
# Compiler options
# -ansi -pedantic -Wall -Wextra -Wconversion -Wshadow
# -Wcast-qual -Wwrite-strings -fPIC -O3 -Wno-long-long"
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi -Wall -Wextra -Wconversion")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wshadow -Wcast-qual -Wwrite-strings")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -O3 -Wno-long-long")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEIGEN_NO_DEBUG -DEIGEN_DEFAULT_TO_ROW_MAJOR")

if(CMAKE_COMPILER_IS_GNUCC)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                    OUTPUT_VARIABLE GCC_VERSION)
    string(REPLACE "\n" "" GCC_VERSION "${GCC_VERSION}")
endif()

message(STATUS "Compiler: ${GCC_VERSION}")
message(STATUS "  CMAKE_C_COMPILER:   " ${CMAKE_C_COMPILER})
message(STATUS "  CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})
message(STATUS "  CMAKE_CXX_FLAGS:    " ${CMAKE_CXX_FLAGS})

#-----------------------------------------------------------------------
# Compilation directives
include_directories(${PROJECT_BINARY_DIR})
configure_file("${PROJECT_SOURCE_DIR}/config.hpp.in"
               "${PROJECT_BINARY_DIR}/config.hpp")

#-----------------------------------------------------------------------
# Debug vs Release
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "Debug compilation")
else()
    message(STATUS "Release compilation")
endif()

#-----------------------------------------------------------------------
# GEM files
gem_find_dir(DirHeaderCommonLib "${PROJECT_SOURCE_DIR}/CommonLib/*.hpp")
include_directories(${DirHeaderCommonLib})
add_subdirectory(CommonLib)
set(EXTRA_LIBS ${EXTRA_LIBS} GemCommon)

if(CUDA_FOUND)
    gem_find_dir(DirHeaderCommonLibCuda "${PROJECT_SOURCE_DIR}/CommonLibCuda/*.cuh")
    include_directories(${DirHeaderCommonLibCuda})
    add_subdirectory(CommonLibCuda)
    set(EXTRA_LIBS ${EXTRA_LIBS} GemCuda)
endif()

if(GEM_FIND_VTK AND VTK_FOUND)
    include(${VTK_USE_FILE})
    set(EXTRA_LIBS ${EXTRA_LIBS} ${VTK_LIBRARIES})
endif()

gem_find_dir(DirHeaderCommonLibData "${PROJECT_SOURCE_DIR}/CommonLibData/*.hpp")
include_directories(${DirHeaderCommonLibData})
add_subdirectory(CommonLibData)
set(EXTRA_LIBS ${EXTRA_LIBS} GemData)

#-----------------------------------------------------------------------
# programs
add_subdirectory(programs)

#-----------------------------------------------------------------------
# Testing using dashboard
enable_testing()
include(CTest)

#-----------------------------------------------------------------------
# Personal testing scripts
add_subdirectory(test)
