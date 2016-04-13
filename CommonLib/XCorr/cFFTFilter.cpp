/***********************************************************************
 *  File:       cFFTFilter.cpp
 *
 *  Purpose:    Implementation of an FFTW3-filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFFTFilter.hpp"

namespace gem {

cFFTFilter::cFFTFilter(void)
{
    _kernelSingle = NULL;
    _kernelDouble = NULL;

    _dataSingle = NULL;
    _dataDouble = NULL;

    _dataSingleFFT = NULL;
    _dataDoubleFFT = NULL;
}

void cFFTFilter::prepare(const float*  const kernel, size_t nrow, size_t ncol)
{
    _nrow    = nrow;
    _ncol    = ncol;
    _fftSize = _nrow * (_ncol / 2 + 1);
    _resSize = _nrow * _ncol;

    // memory
    array_new(_kernelSingle, _fftSize);
    array_memcpy(_kernelSingle, kernel, _fftSize);

    array_new(_dataSingle, _resSize);
    _dataSingleFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * _fftSize);

    // plans
    mutexFFTW.lock();
    _fftPlanFwdSingle = fftwf_plan_dft_r2c_2d((int) _nrow, (int) _ncol,
                                _dataSingle, _dataSingleFFT, FFTW_ESTIMATE);
    _fftPlanInvSingle = fftwf_plan_dft_c2r_2d((int) _nrow, (int) _ncol,
                                _dataSingleFFT, _dataSingle, FFTW_ESTIMATE);
    mutexFFTW.unlock();
}

void cFFTFilter::prepare(const double* const kernel, size_t nrow, size_t ncol)
{
    _nrow    = nrow;
    _ncol    = ncol;
    _fftSize = _nrow * (_ncol / 2 + 1);
    _resSize = _nrow * _ncol;

    // memory
    array_new(_kernelDouble, _fftSize);
    array_memcpy(_kernelDouble, kernel, _fftSize);

    array_new(_dataDouble, _resSize);
    _dataDoubleFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * _fftSize);

    // plans
    mutexFFTW.lock();
    _fftPlanFwdDouble = fftw_plan_dft_r2c_2d((int) _nrow, (int) _ncol,
                                _dataDouble, _dataDoubleFFT, FFTW_ESTIMATE);
    _fftPlanInvDouble = fftw_plan_dft_c2r_2d((int) _nrow, (int) _ncol,
                                _dataDoubleFFT, _dataDouble, FFTW_ESTIMATE);
    mutexFFTW.unlock();
}

void cFFTFilter::perform(float*  const data)
{
    array_memcpy(_dataSingle, data, _resSize);

    fftwf_execute(_fftPlanFwdSingle);

    #pragma omp parallel for
    for (size_t i = 0; i < _fftSize; i++) {
        _dataSingleFFT[i][0] = _dataSingleFFT[i][0] * _kernelSingle[i] / (float) _resSize;
        _dataSingleFFT[i][1] = _dataSingleFFT[i][1] * _kernelSingle[i] / (float) _resSize;
    }

    fftwf_execute(_fftPlanInvSingle);

    array_memcpy(data, _dataSingle, _resSize);
}

void cFFTFilter::perform(double* const data)
{
    array_memcpy(_dataDouble, data, _resSize);

    fftw_execute(_fftPlanFwdDouble);

    #pragma omp parallel for
    for (size_t i = 0; i < _fftSize; i++) {
        _dataDoubleFFT[i][0] = _dataDoubleFFT[i][0] * _kernelDouble[i] / (double) _resSize;
        _dataDoubleFFT[i][1] = _dataDoubleFFT[i][1] * _kernelDouble[i] / (double) _resSize;
    }

    fftw_execute(_fftPlanInvDouble);

    array_memcpy(data, _dataDouble, _resSize);
}

void cFFTFilter::clean(const float* const kernel)
{
    require(kernel != 0);

    array_delete(_kernelSingle);
    array_delete(_dataSingle);

    mutexFFTW.lock();
    fftwf_free(_dataSingleFFT);
    fftwf_destroy_plan(_fftPlanFwdSingle);
    fftwf_destroy_plan(_fftPlanInvSingle);
    mutexFFTW.unlock();
}

void cFFTFilter::clean(const double* const kernel)
{
    require(kernel != 0);

    array_delete(_kernelDouble);
    array_delete(_dataDouble);

    mutexFFTW.lock();
    fftw_free(_dataDoubleFFT);
    fftw_destroy_plan(_fftPlanFwdDouble);
    fftw_destroy_plan(_fftPlanInvDouble);
    mutexFFTW.unlock();
}

} // namespace gem
