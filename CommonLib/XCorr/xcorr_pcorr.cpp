/***********************************************************************
 *  File:       xcorr_pcorr.cpp
 *
 *  Purpose:    Implementation of xcorr-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "fftw_lock.hpp"
#include "xcorr.hpp"

namespace gem {

// size(tplData) = size(refData) = size(resData) = refNrow
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);

    size_t    fftSize = refNrow / 2 + 1;

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftwf_complex *refFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *tplFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *resFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_1d((int) refNrow,
                                const_cast<float*>(refData), refFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_1d((int) refNrow,
                                const_cast<float*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_1d((int) refNrow,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanRef);
    fftwf_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (float) refNrow);

    // obtain the phase correlation array
    fftwf_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(refFFT);
    fftwf_free(tplFFT);
    fftwf_free(resFFT);

    mutexFFTW.unlock();
}

// size(tplData) = size(refData) = size(resData) = refNrow
void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);

    size_t    fftSize = refNrow / 2 + 1;

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftw_complex *refFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *tplFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *resFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_1d((int) refNrow,
                                const_cast<double*>(refData), refFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_1d((int) refNrow,
                                const_cast<double*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_1d((int) refNrow,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanRef);
    fftw_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (double) refNrow);

    // obtain the phase correlation array
    fftw_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(refFFT);
    fftw_free(tplFFT);
    fftw_free(resFFT);

    mutexFFTW.unlock();
}

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol]
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow, size_t refNcol)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);
    assert(refNcol > 0);

    size_t    fftSize = refNrow * (refNcol / 2 + 1);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftwf_complex *refFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *tplFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *resFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                const_cast<float*>(refData), refFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                const_cast<float*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanRef);
    fftwf_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (float) (refNrow*refNcol));

    // obtain the phase correlation array
    fftwf_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(refFFT);
    fftwf_free(tplFFT);
    fftwf_free(resFFT);

    mutexFFTW.unlock();
}

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol]
void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow, size_t refNcol)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);
    assert(refNcol > 0);

    size_t    fftSize = refNrow * (refNcol / 2 + 1);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftw_complex *refFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *tplFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *resFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                const_cast<double*>(refData), refFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                const_cast<double*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanRef);
    fftw_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (double) (refNrow*refNcol));

    // obtain the phase correlation array
    fftw_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(refFFT);
    fftw_free(tplFFT);
    fftw_free(resFFT);

    mutexFFTW.unlock();
}

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol refNsec]
void pcorr(const float*  const tplData,
           const float*  const refData,
                 float*  const resData,
           size_t refNrow, size_t refNcol, size_t refNsec)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);
    assert(refNcol > 0);
    assert(refNsec > 0);

    size_t    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftwf_complex *refFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *tplFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    fftwf_complex *resFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                const_cast<float*>(refData), refFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                const_cast<float*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanRef);
    fftwf_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (float) (refNrow*refNcol*refNsec));

    // obtain the phase correlation array
    fftwf_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(refFFT);
    fftwf_free(tplFFT);
    fftwf_free(resFFT);

    mutexFFTW.unlock();
}

// size(tplData) = size(refData) = size(resData) = [refNrow refNcol refNsec]
void pcorr(const double* const tplData,
           const double* const refData,
                 double* const resData,
           size_t refNrow, size_t refNcol, size_t refNsec)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(refNrow > 0);
    assert(refNcol > 0);
    assert(refNsec > 0);

    size_t    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    fftw_complex *refFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *tplFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    fftw_complex *resFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                const_cast<double*>(refData), refFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                const_cast<double*>(tplData), tplFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resFFT, resData, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanRef);
    fftw_execute(fftPlanTpl);

    // obtain the cross power spectrum and normalize the data
    pcorrModulateAndNormalize(tplFFT, refFFT, resFFT,
                              fftSize,
                              (double) (refNrow*refNcol*refNsec));

    // obtain the phase correlation array
    fftw_execute(fftPlanRes);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(refFFT);
    fftw_free(tplFFT);
    fftw_free(resFFT);

    mutexFFTW.unlock();
}

} // namespace gem
