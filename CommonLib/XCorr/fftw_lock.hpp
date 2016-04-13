/***********************************************************************
 *  File:       fftw_lock.hpp
 *
 *  Purpose:    Header file for fftw-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include <fftw3.h>

#ifndef __GEM_LIB_FFTW_LOCK_HPP__
#define __GEM_LIB_FFTW_LOCK_HPP__

#ifdef __GEM_USE_STD11__
    #include <mutex>

namespace gem {

    extern  std::mutex      mutexFFTW;
    extern  std::mutex      mutexMergeResult;

    #define FFTW_LOCK(command) \
    { \
        std::lock_guard<std::mutex>     lockHandle(mutexFFTW); \
        command; \
    }

} // namespace gem

#else
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wconversion"
    #include "boost/thread.hpp"
    #pragma GCC diagnostic pop

namespace gem {

    extern  boost::mutex    mutexFFTW;
    extern  boost::mutex    mutexMergeResult;

    #define FFTW_LOCK(command) \
    { \
        boost::mutex::scoped_lock       lockHandle(mutexFFTW); \
        command; \
    }

} // namespace gem

#endif

#endif
