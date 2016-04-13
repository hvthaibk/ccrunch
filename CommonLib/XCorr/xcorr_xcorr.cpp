/***********************************************************************
 *  File:       xcorr_xcorr.cpp
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

/***************************************************
 * Cross-correlation
 **************************************************/

// 1D float
void xcorr(const float* const tplData, size_t tplNrow,
           const float* const refData, size_t refNrow,
                 float* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t          fftSize, resNrow, resSize;
    float           *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftwf_complex   *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow / 2 + 1;
    resSize = resNrow;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_1d((int) resNrow,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_1d((int) resNrow,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_1d((int) resNrow,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow,
                      tplDataPadded, resNrow,
                                     0);
            array_pad(refData,       refNrow,
                      refDataPadded, resNrow,
                                     tplNrow-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow,
                      tplDataPadded, resNrow,
                                     0);
            array_pad(refData,       refNrow,
                      refDataPadded, resNrow,
                                     0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanTpl);
    fftwf_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) resSize);

    // obtain the correlation array
    fftwf_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow,
                       resData,       resNrow,
                                      0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,
                       resData,       resNrow-tplNrow+1,
                                      0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 1D double
void xcorr(const double* const tplData, size_t tplNrow,
           const double* const refData, size_t refNrow,
                 double* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t          fftSize, resNrow, resSize;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow / 2 + 1;
    resSize = resNrow;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_1d((int) resNrow,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_1d((int) resNrow,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_1d((int) resNrow,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow,
                      tplDataPadded, resNrow,
                                     0);
            array_pad(refData,       refNrow,
                      refDataPadded, resNrow,
                                     tplNrow-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow,
                      tplDataPadded, resNrow,
                                     0);
            array_pad(refData,       refNrow,
                      refDataPadded, resNrow,
                                     0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanTpl);
    fftw_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) resSize);

    // obtain the correlation array
    fftw_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow,
                       resData,       resNrow,
                                      0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,
                       resData,       resNrow-tplNrow+1,
                                      0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(tplDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 2D float
void xcorr(const float* const tplData, size_t tplNrow, size_t tplNcol,
           const float* const refData, size_t refNrow, size_t refNcol,
                 float* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t          fftSize, resNrow, resNcol, resSize;
    float           *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftwf_complex   *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            resNcol = refNcol;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * (resNcol / 2 + 1);
    resSize = resNrow * resNcol;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_2d((int) resNrow, (int) resNcol,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_2d((int) resNrow, (int) resNcol,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_2d((int) resNrow, (int) resNcol,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow, tplNcol,
                      tplDataPadded, resNrow, resNcol,
                                     0,       0);
            array_pad(refData,       refNrow,   refNcol,
                      refDataPadded, resNrow,   resNcol,
                                     tplNrow-1, tplNcol-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow, tplNcol,
                      tplDataPadded, resNrow, resNcol,
                                     0,       0);
            array_pad(refData,       refNrow, refNcol,
                      refDataPadded, resNrow, resNcol,
                                     0,       0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanTpl);
    fftwf_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) resSize);

    // obtain the correlation array
    fftwf_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow, resNcol,
                       resData,       resNrow, resNcol,
                                      0,       0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,           resNcol,
                       resData,       resNrow-tplNrow+1, resNcol-tplNcol+1,
                                      0,                 0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 2D double
void xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol,
           const double* const refData, size_t refNrow, size_t refNcol,
                 double* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t          fftSize, resNrow, resNcol, resSize;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            resNcol = refNcol;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * (resNcol / 2 + 1);
    resSize = resNrow * resNcol;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_2d((int) resNrow, (int) resNcol,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_2d((int) resNrow, (int) resNcol,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_2d((int) resNrow, (int) resNcol,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow, tplNcol,
                      tplDataPadded, resNrow, resNcol,
                                     0,       0);
            array_pad(refData,       refNrow,   refNcol,
                      refDataPadded, resNrow,   resNcol,
                                     tplNrow-1, tplNcol-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow, tplNcol,
                      tplDataPadded, resNrow, resNcol,
                                     0,       0);
            array_pad(refData,       refNrow, refNcol,
                      refDataPadded, resNrow, resNcol,
                                     0,       0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanTpl);
    fftw_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) resSize);

    // obtain the correlation array
    fftw_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow, resNcol,
                       resData,       resNrow, resNcol,
                                      0,       0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,           resNcol,
                       resData,       resNrow-tplNrow+1, resNcol-tplNcol+1,
                                      0,                 0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(tplDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 3D float
void xcorr(const float* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           const float* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                 float* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t          fftSize, resNrow, resNcol, resNsec, resSize;
    float           *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftwf_complex   *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            resNsec = refNsec + tplNsec - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            resNcol = refNcol;
            resNsec = refNsec;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * resNcol * (resNsec / 2 + 1);
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // setup FFTW plans
    fftwf_plan fftPlanTpl = fftwf_plan_dft_r2c_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef = fftwf_plan_dft_r2c_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRes = fftwf_plan_dft_c2r_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                      tplDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            array_pad(refData,       refNrow,   refNcol,   refNsec,
                      refDataPadded, resNrow,   resNcol,   resNsec,
                                     tplNrow-1, tplNcol-1, tplNsec-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                      tplDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            array_pad(refData,       refNrow, refNcol, refNsec,
                      refDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftwf_execute(fftPlanTpl);
    fftwf_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) resSize);

    // obtain the correlation array
    fftwf_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow, resNcol, resNsec,
                       resData,       resNrow, resNcol, resNsec,
                                      0,       0,       0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,           resNcol,           resNsec,
                       resData,       resNrow-tplNrow+1, resNcol-tplNcol+1, resNsec-tplNsec+1,
                                      0,                 0,                 0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRes);
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 3D double
void xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           const double* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                 double* const resData,
           eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t          fftSize, resNrow, resNcol, resNsec, resSize;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            resNsec = refNsec + tplNsec - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            resNcol = refNcol;
            resNsec = refNsec;
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * resNcol * (resNsec / 2 + 1);
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays
    array_new(tplDataPadded, resSize);
    array_new(refDataPadded, resSize);
    array_new(resDataPadded, resSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // setup FFTW plans
    fftw_plan fftPlanTpl = fftw_plan_dft_r2c_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef = fftw_plan_dft_r2c_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRes = fftw_plan_dft_c2r_3d((int) resNrow, (int) resNcol, (int) resNsec,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // load data to FFTW inputs
    switch (shape) {
        case XCORR_RES_FULL:
            array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                      tplDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            array_pad(refData,       refNrow,   refNcol,   refNsec,
                      refDataPadded, resNrow,   resNcol,   resNsec,
                                     tplNrow-1, tplNcol-1, tplNsec-1);
            break;
        case XCORR_RES_VALID:
            array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                      tplDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            array_pad(refData,       refNrow, refNcol, refNsec,
                      refDataPadded, resNrow, resNcol, resNsec,
                                     0,       0,       0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // obtain the FFT of ref and tpl
    fftw_execute(fftPlanTpl);
    fftw_execute(fftPlanRef);

    // obtain the cross power spectrum and normalize the data
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) resSize);

    // obtain the correlation array
    fftw_execute(fftPlanRes);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            array_crop(resDataPadded, resNrow, resNcol, resNsec,
                       resData,       resNrow, resNcol, resNsec,
                                      0,       0,       0);
            break;
        case XCORR_RES_VALID:
            array_crop(resDataPadded, resNrow,           resNcol,           resNsec,
                       resData,       resNrow-tplNrow+1, resNcol-tplNcol+1, resNsec-tplNsec+1,
                                      0,                 0,                 0);
            break;
        default:
            ERROR("xcorr", "unsupported xcorr mode");
    }

    // deallocate arrays
    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRes);
    fftw_free(tplDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

} // namespace gem
