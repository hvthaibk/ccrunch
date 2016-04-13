/***********************************************************************
 *  File:       cuSystem.hpp
 *
 *  Purpose:    Header file for a CUDA-info class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUSYSTEM_HPP__
#define __GEM_CUSYSTEM_HPP__

#include "config.hpp"
#include "macro.hpp"

#include <iostream>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>

namespace gem {

class cuSystem
{
private:
    int                 _numDevice;
    int                 _verDriver;
    int                 _verRunTime;
    cudaDeviceProp*     _deviceProp;

public:
     cuSystem(void);
    ~cuSystem();

    void resetDevice(void);

    void printDevInfo(void);
    void printMemInfo(int devID);
    void printThroughput(size_t numElements, float time);

    int  getNumDevice     (void) { return _numDevice;  }
    int  getDriverVersion (void) { return _verDriver;  }
    int  getRuntimeVersion(void) { return _verRunTime; }
};

} // namespace gem

#endif
