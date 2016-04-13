/***********************************************************************
 *  File:       cuFFTFilter.cuh
 *
 *  Purpose:    Header file for an CUFFT-filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUFFT_FILTER_HPP__
#define __GEM_CUFFT_FILTER_HPP__

#include "array.cuh"
#include <cufft.h>

namespace gem {

class cuFFTFilter
{
private:
    // single-precision
    float           *_kernelSingle,
                    *_dataSingle;
    cuFloatComplex  *_dataSingleFFT;

    // double-precision
    double          *_kernelDouble,
                    *_dataDouble;
    cuDoubleComplex *_dataDoubleFFT;

    // plan handle
    cufftHandle     _fftPlanFwd,
                    _fftPlanInv;

    // size information
    size_t          _nrow, _ncol, _nsec;
    size_t          _fftSize, _resSize;

public:
     cuFFTFilter(void);
    ~cuFFTFilter()  {};

    // 2D
    void    prepare(const float*  const kernel, size_t nrow, size_t ncol);
    void    prepare(const double* const kernel, size_t nrow, size_t ncol);

    // perform FFT and IFFT
    void    perform(float*  const data);
    void    perform(double* const data);

    // garbage collection
    void    clean  (const float*  const kernel);
    void    clean  (const double* const kernel);
};

} // namespace gem

#endif
