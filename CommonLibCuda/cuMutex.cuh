/***********************************************************************
 *  File:       cMutex.cuh
 *
 *  Purpose:    Header file for mutex object
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CMUTEXCUDA_CUH__
#define __GEM_CMUTEXCUDA_CUH__

#include "macro.hpp"

#include <cuda.h>

namespace gem {

struct cMutexCUDA
{
    int     *mutex;

    cMutexCUDA(void)
    {
        int state = 0;

        cuda_arrayDev_new(mutex, 1);
        cuda_memcpy_h2d(mutex, &state, 1));
    }

    ~cMutexCUDA(void)
    {
        cuda_arrayDev_delete(mutex);
    }

    __device__ void lock(void)
    {
        while (atomicCAS(mutex, 0, 1) != 0);
    }

    __device__ void unlock(void)
    {
        atomicExch(mutex, 0);
    }
}

} // namespace gem

#endif
