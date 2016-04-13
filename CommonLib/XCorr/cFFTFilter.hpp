/***********************************************************************
 *  File:       cFFTFilter.hpp
 *
 *  Purpose:    Header file for an FFTW3-filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CFFTW3_FILTER_HPP__
#define __GEM_CFFTW3_FILTER_HPP__

#include "array.hpp"

#include "fftw_lock.hpp"

namespace gem {

class cFFTFilter
{
private:
    // single-precision
    float           *_kernelSingle,
                    *_dataSingle;
    fftwf_complex   *_dataSingleFFT;

    fftwf_plan      _fftPlanFwdSingle,
                    _fftPlanInvSingle;

    // double-precision
    double          *_kernelDouble,
                    *_dataDouble;
    fftw_complex    *_dataDoubleFFT;

    fftw_plan       _fftPlanFwdDouble,
                    _fftPlanInvDouble;

    // size information
    size_t          _nrow, _ncol, _nsec;
    size_t          _fftSize, _resSize;

public:
     cFFTFilter(void);
    ~cFFTFilter()  {};

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
