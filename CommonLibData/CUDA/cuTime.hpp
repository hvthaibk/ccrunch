/***********************************************************************
 *  File:       cuTime.hpp
 *
 *  Purpose:    Header file for a high-resolution timer
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUTIME_HPP__
#define __GEM_CUTIME_HPP__

#include "config.hpp"
#include "macro.hpp"

#include <iostream>
#include <cstdio>
#include <cuda.h>
#include <cuda_runtime.h>

/* USAGE
 *
 * avoid using this timer object in production code since
 * it might cause segmentation error
 *      cuTime      timerCUDA;
 *
 * use pointer instead
 *      cuTime      *timerCUDA;

#ifdef __GEM_USE_CUDA__
    cuTime      *timerCUDA;
#endif
    timerCUDA = new cuTime;
    delete timerCUDA;
*/

namespace gem {

class cuTime
{
private:
    cudaEvent_t  timeStart, timeStop;
    float        timeElapsed;

public:
     cuTime(void);
    ~cuTime();

    void tic(cudaStream_t s = 0);
    void toc(std::string str = "", cudaStream_t s = 0);

    float getElapsedTime();
};

} // namespace gem

#endif
