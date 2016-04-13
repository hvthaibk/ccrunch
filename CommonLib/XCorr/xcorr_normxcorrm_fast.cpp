/***********************************************************************
 *  File:       xcorr_normxcorrm_fast.cpp
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
 * Combination for fast computation
 **************************************************/

// 1D float
void normxcorrm_fast(float *tplData, size_t tplNrow,
                     float *refData, size_t refNrow,
                     float *resData,
                     float *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t           resNrow, tplSize, refSize, fftSize;
    float            *tplDataTmp = NULL;
    float            *tplDataPadded = NULL,
                     *refDataPadded = NULL,
                     *mskDataPadded = NULL,
                     *resDataPadded = NULL;
    fftwf_complex    *tplDataPaddedFFT = NULL,
                     *refDataPaddedFFT = NULL,
                     *refDataPaddedFFT2 = NULL,
                     *mskDataPaddedFFT = NULL,
                     *resDataPaddedFFT = NULL;
    float            *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t           mskArea;

    resNrow = refNrow - tplNrow + 1;
    tplSize = tplNrow;
    refSize = refNrow;
    fftSize = refNrow / 2 + 1;

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT2 = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    mskDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // FFT plans (FLOAT)
    fftwf_plan fftPlanTpl  = fftwf_plan_dft_r2c_1d((int) refNrow,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef  = fftwf_plan_dft_r2c_1d((int) refNrow,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef2 = fftwf_plan_dft_r2c_1d((int) refNrow,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftwf_plan fftPlanMsk  = fftwf_plan_dft_r2c_1d((int) refNrow,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv1 = fftwf_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv2 = fftwf_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv3 = fftwf_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(float));
    memcpy(refDataPadded, refData, refSize*sizeof(float));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftwf_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftwf_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow,
              tplDataPadded, refNrow,
                             0);
    fftwf_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow,
              mskDataPadded, refNrow,
                             0);
    fftwf_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow,
               resData,       resNrow,
                              0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRef2);
    fftwf_destroy_plan(fftPlanMsk);
    fftwf_destroy_plan(fftPlanInv1);
    fftwf_destroy_plan(fftPlanInv2);
    fftwf_destroy_plan(fftPlanInv3);

    fftwf_free(tplDataPaddedFFT);
    fftwf_free(mskDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(refDataPaddedFFT2);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 1D double
void normxcorrm_fast(double *tplData, size_t tplNrow,
                     double *refData, size_t refNrow,
                     double *resData,
                     double *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t          resNrow, tplSize, refSize, fftSize;
    double          *tplDataTmp = NULL;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *mskDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *refDataPaddedFFT2 = NULL,
                    *mskDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;
    double          *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t          mskArea;

    resNrow = refNrow - tplNrow + 1;
    tplSize = tplNrow;
    refSize = refNrow;
    fftSize = refNrow / 2 + 1;

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    mskDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // FFT plans (DOUBLE)
    fftw_plan fftPlanTpl  = fftw_plan_dft_r2c_1d((int) refNrow,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef  = fftw_plan_dft_r2c_1d((int) refNrow,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef2 = fftw_plan_dft_r2c_1d((int) refNrow,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftw_plan fftPlanMsk  = fftw_plan_dft_r2c_1d((int) refNrow,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanInv1 = fftw_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftw_plan fftPlanInv2 = fftw_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftw_plan fftPlanInv3 = fftw_plan_dft_c2r_1d((int) refNrow,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(double));
    memcpy(refDataPadded, refData, refSize*sizeof(double));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftw_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftw_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow,
              tplDataPadded, refNrow,
                             0);
    fftw_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow,
              mskDataPadded, refNrow,
                             0);
    fftw_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow,
               resData,       resNrow,
                              0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRef2);
    fftw_destroy_plan(fftPlanMsk);
    fftw_destroy_plan(fftPlanInv1);
    fftw_destroy_plan(fftPlanInv2);
    fftw_destroy_plan(fftPlanInv3);

    fftw_free(tplDataPaddedFFT);
    fftw_free(mskDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(refDataPaddedFFT2);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 2D float
void normxcorrm_fast(float *tplData, size_t tplNrow, size_t tplNcol,
                     float *refData, size_t refNrow, size_t refNcol,
                     float *resData,
                     float *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t           resNrow, resNcol, tplSize, refSize, fftSize;
    float            *tplDataTmp = NULL;
    float            *tplDataPadded = NULL,
                     *refDataPadded = NULL,
                     *mskDataPadded = NULL,
                     *resDataPadded = NULL;
    fftwf_complex    *tplDataPaddedFFT = NULL,
                     *refDataPaddedFFT = NULL,
                     *refDataPaddedFFT2 = NULL,
                     *mskDataPaddedFFT = NULL,
                     *resDataPaddedFFT = NULL;
    float            *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t           mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    fftSize = refNrow * (refNcol / 2 + 1);

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT2 = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    mskDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // FFT plans (FLOAT)
    fftwf_plan fftPlanTpl  = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef  = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef2 = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftwf_plan fftPlanMsk  = fftwf_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv1 = fftwf_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv2 = fftwf_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv3 = fftwf_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(float));
    memcpy(refDataPadded, refData, refSize*sizeof(float));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftwf_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftwf_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow, tplNcol,
              tplDataPadded, refNrow, refNcol,
                             0,       0);
    fftwf_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow, tplNcol,
              mskDataPadded, refNrow, refNcol,
                             0,       0);
    fftwf_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow, refNcol,
               resData,       resNrow, resNcol,
                              0,       0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRef2);
    fftwf_destroy_plan(fftPlanMsk);
    fftwf_destroy_plan(fftPlanInv1);
    fftwf_destroy_plan(fftPlanInv2);
    fftwf_destroy_plan(fftPlanInv3);

    fftwf_free(tplDataPaddedFFT);
    fftwf_free(mskDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(refDataPaddedFFT2);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 2D double
void normxcorrm_fast(double *tplData, size_t tplNrow, size_t tplNcol,
                     double *refData, size_t refNrow, size_t refNcol,
                     double *resData,
                     double *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t          resNrow, resNcol, tplSize, refSize, fftSize;
    double          *tplDataTmp = NULL;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *mskDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *refDataPaddedFFT2 = NULL,
                    *mskDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;
    double          *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t          mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    fftSize = refNrow * (refNcol / 2 + 1);

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    mskDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // FFT plans (DOUBLE)
    fftw_plan fftPlanTpl  = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef  = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef2 = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftw_plan fftPlanMsk  = fftw_plan_dft_r2c_2d((int) refNrow, (int) refNcol,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanInv1 = fftw_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftw_plan fftPlanInv2 = fftw_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftw_plan fftPlanInv3 = fftw_plan_dft_c2r_2d((int) refNrow, (int) refNcol,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(double));
    memcpy(refDataPadded, refData, refSize*sizeof(double));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftw_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftw_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow, tplNcol,
              tplDataPadded, refNrow, refNcol,
                             0,       0);
    fftw_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow, tplNcol,
              mskDataPadded, refNrow, refNcol,
                             0,       0);
    fftw_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow, refNcol,
               resData,       resNrow, resNcol,
                              0,       0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRef2);
    fftw_destroy_plan(fftPlanMsk);
    fftw_destroy_plan(fftPlanInv1);
    fftw_destroy_plan(fftPlanInv2);
    fftw_destroy_plan(fftPlanInv3);

    fftw_free(tplDataPaddedFFT);
    fftw_free(mskDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(refDataPaddedFFT2);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 3D float
void normxcorrm_fast(float *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     float *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     float *resData,
                     float *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t           resNrow, resNcol, resNsec, tplSize, refSize, fftSize;
    float            *tplDataTmp = NULL;
    float            *tplDataPadded = NULL,
                     *refDataPadded = NULL,
                     *mskDataPadded = NULL,
                     *resDataPadded = NULL;
    fftwf_complex    *tplDataPaddedFFT = NULL,
                     *refDataPaddedFFT = NULL,
                     *refDataPaddedFFT2 = NULL,
                     *mskDataPaddedFFT = NULL,
                     *resDataPaddedFFT = NULL;
    float            *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t           mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    resNsec = refNsec - tplNsec + 1;
    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT2 = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    mskDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * fftSize);

    // FFT plans (FLOAT)
    fftwf_plan fftPlanTpl  = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef  = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanRef2 = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftwf_plan fftPlanMsk  = fftwf_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv1 = fftwf_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv2 = fftwf_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftwf_plan fftPlanInv3 = fftwf_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(float));
    memcpy(refDataPadded, refData, refSize*sizeof(float));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftwf_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftwf_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow, tplNcol, tplNsec,
              tplDataPadded, refNrow, refNcol, refNsec,
                             0,       0,       0);
    fftwf_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow, tplNcol, tplNsec,
              mskDataPadded, refNrow, refNcol, refNsec,
                             0,       0,       0);
    fftwf_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow, refNcol, refNsec,
               resData,       resNrow, resNcol, resNsec,
                              0,       0,       0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRef2);
    fftwf_destroy_plan(fftPlanMsk);
    fftwf_destroy_plan(fftPlanInv1);
    fftwf_destroy_plan(fftPlanInv2);
    fftwf_destroy_plan(fftPlanInv3);

    fftwf_free(tplDataPaddedFFT);
    fftwf_free(mskDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(refDataPaddedFFT2);
    fftwf_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

// 3D double
void normxcorrm_fast(double *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                     double *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                     double *resData,
                     double *mskData)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t          resNrow, resNcol, resNsec, tplSize, refSize, fftSize;
    double          *tplDataTmp = NULL;
    double          *tplDataPadded = NULL,
                    *refDataPadded = NULL,
                    *mskDataPadded = NULL,
                    *resDataPadded = NULL;
    fftw_complex    *tplDataPaddedFFT = NULL,
                    *refDataPaddedFFT = NULL,
                    *refDataPaddedFFT2 = NULL,
                    *mskDataPaddedFFT = NULL,
                    *resDataPaddedFFT = NULL;
    double          *fftTmp1 = NULL, *fftTmp2 = NULL;
    size_t          mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    resNsec = refNsec - tplNsec + 1;
    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    // allocate arrays
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    mutexFFTW.lock();

    // allocate FFTW input and output arrays
    tplDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    mskDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fftSize);

    // FFT plans (DOUBLE)
    fftw_plan fftPlanTpl  = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                tplDataPadded, tplDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef  = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                refDataPadded, refDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanRef2 = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                refDataPadded, refDataPaddedFFT2, FFTW_ESTIMATE);
    fftw_plan fftPlanMsk  = fftw_plan_dft_r2c_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                mskDataPadded, mskDataPaddedFFT, FFTW_ESTIMATE);
    fftw_plan fftPlanInv1 = fftw_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, resDataPadded, FFTW_ESTIMATE);
    fftw_plan fftPlanInv2 = fftw_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, fftTmp1, FFTW_ESTIMATE);
    fftw_plan fftPlanInv3 = fftw_plan_dft_c2r_3d((int) refNrow, (int) refNcol, (int) refNsec,
                                resDataPaddedFFT, fftTmp2, FFTW_ESTIMATE);

    mutexFFTW.unlock();

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    memcpy(tplDataTmp,    tplData, tplSize*sizeof(double));
    memcpy(refDataPadded, refData, refSize*sizeof(double));

    // mask size and template normalization
    mskArea = xcorrNormTemplateNCC(tplDataTmp, tplSize, mskData);

    // reference
    fftw_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftw_execute(fftPlanRef2);

    // template
    array_pad(tplDataTmp,    tplNrow, tplNcol, tplNsec,
              tplDataPadded, refNrow, refNcol, refNsec,
                             0,       0,       0);
    fftw_execute(fftPlanTpl);

    // mask
    array_pad(mskData,       tplNrow, tplNcol, tplNsec,
              mskDataPadded, refNrow, refNcol, refNsec,
                             0,       0,       0);
    fftw_execute(fftPlanMsk);

    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv1);

    // xcorr 2
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv2);

    // xcorr 3
    xcorrModulateAndNormalize(mskDataPaddedFFT,
                              refDataPaddedFFT2,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv3);

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2, mskArea, refSize);

    // crop
    array_crop(resDataPadded, refNrow, refNcol, refNsec,
               resData,       resNrow, resNcol, resNsec,
                              0,       0,       0);

    // deallocate arrays
    array_delete(tplDataTmp);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    mutexFFTW.lock();

    // deallocate FFTW arrays and plans
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRef2);
    fftw_destroy_plan(fftPlanMsk);
    fftw_destroy_plan(fftPlanInv1);
    fftw_destroy_plan(fftPlanInv2);
    fftw_destroy_plan(fftPlanInv3);

    fftw_free(tplDataPaddedFFT);
    fftw_free(mskDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(refDataPaddedFFT2);
    fftw_free(resDataPaddedFFT);

    mutexFFTW.unlock();
}

} // namespace gem
