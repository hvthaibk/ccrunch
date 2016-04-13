/***********************************************************************
 *  File:       cuTime.cu
 *
 *  Purpose:    Implementation of a high-resolution timer
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuTime.hpp"

namespace gem {

cuTime::cuTime(void)
{
    CUDA_SAFE_CALL(cudaEventCreate(&timeStart));
    CUDA_SAFE_CALL(cudaEventCreate(&timeStop));
}

cuTime::~cuTime()
{
    CUDA_SAFE_CALL(cudaEventDestroy(timeStart));
    CUDA_SAFE_CALL(cudaEventDestroy(timeStop));
}

void cuTime::tic(cudaStream_t s)
{
    CUDA_SAFE_CALL(cudaEventRecord(timeStart, s));
}

void cuTime::toc(std::string str, cudaStream_t s)
{
    CUDA_SAFE_CALL(cudaEventRecord(timeStop, s));

    CUDA_SAFE_CALL(cudaEventSynchronize(timeStop));
    CUDA_SAFE_CALL(cudaEventElapsedTime(&timeElapsed, timeStart, timeStop));

    // convert from milliseconds to seconds
    timeElapsed /= 1000;

    if (str.empty()) {
        std::cout << "     elapsed time: "
                  << timeElapsed << std::endl;
    }
    else if (str.compare("noshow") == 0) {
    }
    else {
        std::cout << str << timeElapsed << std::endl;
    }
}

float cuTime::getElapsedTime()
{
    return timeElapsed;
}

} // namespace gem
