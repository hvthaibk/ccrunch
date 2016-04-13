/***********************************************************************
 *  File:       fftw_lock.cpp
 *
 *  Purpose:    Implementation of fftw-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "fftw_lock.hpp"

#ifdef __GEM_USE_STD11__
    #include <mutex>

namespace gem {

    std::mutex      mutexFFTW;
    std::mutex      mutexMergeResult;

} // namespace gem

#else
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wconversion"
    #include "boost/thread.hpp"
    #pragma GCC diagnostic pop

namespace gem {

    boost::mutex    mutexFFTW;
    boost::mutex    mutexMergeResult;

} // namespace gem

#endif
