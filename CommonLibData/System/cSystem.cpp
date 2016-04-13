/***********************************************************************
 *  File:       cSystem.cpp
 *
 *  Purpose:    Implementation of a system-info class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cSystem.hpp"

#include "fftw3.h"
#include "array.hpp"

#ifdef __GEM_USE_CUDA__
#include "array.cuh"
#include <cuda_runtime.h>
#endif

#ifdef __GEM_USE_OPENMP__
#include "omp.h"
#endif

#ifdef __unix__
#include <unistd.h>
#elif __win32__
#include <Windows.h>
#include <tchar.h>
#endif

namespace gem {

const std::string cSystem::getComputerName(void)
{
#ifdef __unix__
    char            res[255];
    std::string     strName;

    memset(res, 0, 255);
    if (gethostname(res, 255) == 0) {
        strName = res;
    }
    else {
        strName = "Unknown hostname";
    }

    return strName;
#elif __win32__
    char            res[255];
    std::string     strName;
    TCHAR           infoBuf[255];
    DWORD           bufCharCount = 255;

    memset(Name, 0, 255);
    if (GetComputerName(infoBuf, &bufCharCount)) {
        for(size_t i = 0; i < 255; i++) {
            res[i] = infoBuf[i];
        }

        strName = res;
    }
    else {
        strName = "Unknown hostname";
    }

#else
    ERROR("getComputerName", "this OS is not yest supported");
#endif
}

unsigned int cSystem::getNumProcessors(void)
{
#ifdef __unix__
    FILE     *fp;
    char     res[128] = {0};

    if (!(fp = popen("/bin/cat /proc/cpuinfo | grep -c '^processor'", "r"))) {
        std::cout << "getNumProcessors(): popen() reading error " << std::endl;
    }
    if (!(fgets(res, sizeof(res), fp))) {
        std::cout << "getNumProcessors(): fgets() reading error " << std::endl;
    }
    fclose(fp);

    return atoi(res);
#else
    ERROR("getNumProcessors", "this OS is not yest supported");
#endif
}

unsigned int cSystem::getNumCoresPerCPU(void)
{
#ifdef __unix__
    std::ifstream   cpuinfo("/proc/cpuinfo");
    std::string     line;
    size_t          start, end;

    while (!cpuinfo.eof()) {
        std::getline(cpuinfo,line);

        if (!line.size()) {
            continue;
        }
        if (line.find("cpu cores") != 0) {
            continue;
        }

        for (start = line.find(':'); line[start] < '0' || line[start] > '9'; start++);
        for (end = start; line[end] >= '0' && line[end] <= '9'; end++);
        line = line.substr(start, end-start);

        break;
    }

    return atoi(line.c_str());
#else
    ERROR("getNumCoresPerCPU", "this OS is not yest supported");
#endif
}

unsigned int cSystem::getNumGPUs(void)
{
    int     nGPU;

#ifdef __GEM_USE_CUDA__
    CUDA_SAFE_CALL(cudaGetDeviceCount(&nGPU));
#else
    nGPU = 0;
#endif

    return (unsigned int) nGPU;
}

bool cSystem::getEndianness(void)
{
    unsigned int    n = 1;

    // little endian if true
    if(*(char *)&n == 1) {
        return true;
    }
    else {
        return false;
    }
}

void cSystem::getCacheInfo(void)
{
    long    size, assoc, line;

    size  = sysconf(_SC_LEVEL1_ICACHE_SIZE);
    assoc = sysconf(_SC_LEVEL1_ICACHE_ASSOC);
    line  = sysconf(_SC_LEVEL1_ICACHE_LINESIZE);
    printf("level 1 icache size = %ldK, assoc = %ld, line size = %ld\n",
           size, assoc, line);

    size  = sysconf(_SC_LEVEL1_DCACHE_SIZE);
    assoc = sysconf(_SC_LEVEL1_DCACHE_ASSOC);
    line  = sysconf(_SC_LEVEL1_DCACHE_LINESIZE);
    printf("level 1 dcache size = %ldK, assoc = %ld, line size = %ld\n",
           size, assoc, line);

    size  = sysconf(_SC_LEVEL2_CACHE_SIZE);
    assoc = sysconf(_SC_LEVEL2_CACHE_ASSOC);
    line  = sysconf(_SC_LEVEL2_CACHE_LINESIZE);
    printf("level 2 cache size = %ldK, assoc = %ld, line size = %ld\n",
           size, assoc, line);

    size  = sysconf(_SC_LEVEL3_CACHE_SIZE);
    assoc = sysconf(_SC_LEVEL3_CACHE_ASSOC);
    line  = sysconf(_SC_LEVEL3_CACHE_LINESIZE);
    printf("level 3 cache size = %ldK, assoc = %ld, line size = %ld\n",
           size, assoc, line);
}

void cSystem::startOpenMP(void)
{
#ifdef __GEM_USE_OPENMP__
    if (omp_get_num_procs() == 1) {
        WARNING("startOpenMP", "cannot start OpenMP since there is only 1 online processor");
    }
    else {
        omp_set_num_threads((int) std::ceil((double) omp_get_num_procs()/2));
        omp_set_nested(false);

        std::cout.precision(0);
        std::cout << "OpenMP is turned on: " << std::ceil(omp_get_num_procs()/2) << " threads\n\n";
    }
#else
    std::cout << "OpenMP is not supported\n\n";
#endif
}

void cSystem::stopOpenMP(void)
{
#ifdef __GEM_USE_OPENMP__
    omp_set_num_threads(1);
    std::cout << "\nOpenMP is turned off\n\n";
#else
    std::cout << "\nOpenMP is not supported\n\n";
#endif
}

void cSystem::startNthreadsFFTW(void)
{
    require( fftw_init_threads() != 0, "void cSystem::startNthreadsFFTW(void)");
    require(fftwf_init_threads() != 0, "void cSystem::startNthreadsFFTW(void)");

     fftw_plan_with_nthreads(getNumProcessors());
    fftwf_plan_with_nthreads(getNumProcessors());

    std::cout << "FFTW multithreading is turned on: " << getNumProcessors() << " threads\n\n";
}

void cSystem::stopNthreadsFFTW(void)
{
     fftw_cleanup_threads();
    fftwf_cleanup_threads();

    std::cout << "FFTW multithreading is turned off\n\n";
}

void cSystem::checkComputingResource(unsigned int nGPU, unsigned int nCPU)
{
    std::string     message;
    std::string     hostname = getComputerName();
    unsigned int    nGPUmax = 0, nCPUmax = 0;

    nGPUmax = getNumGPUs();
    message = "(" + hostname + ") ";
    message += "has only " + num2str(nGPUmax) + " GPUs!";
    require(nGPU <= nGPUmax, message.c_str());

    nCPUmax = getNumProcessors();
    message = "(" + hostname + ") ";
    message += "has only " + num2str(nCPUmax) + " CPU cores!";
    require(nCPU <= nCPUmax, message.c_str());
}

void cSystem::checkMemoryLeak(void)
{
    array_memcheck();

#ifdef __GEM_USE_CUDA__

    cuda_array_memcheck();

#endif
}

} // namespace gem
