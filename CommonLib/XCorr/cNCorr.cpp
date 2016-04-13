/***********************************************************************
 *  File:       cNCorr.cpp
 *
 *  Purpose:    Implementation of NCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cNCorr.hpp"

namespace gem {

cNCorrSingle::cNCorrSingle(size_t refNrowInput,
                           size_t refNcolInput,
                           size_t refNsecInput)
{
    assert(refNrowInput > 0 && refNcolInput > 0 && refNsecInput > 0);

    refNrow = refNrowInput;
    refNcol = refNcolInput;
    refNsec = refNsecInput;
    refSize = refNrow * refNcol * refNsec;

    if (refNsec == 1) {       // 2D
        fftSize = refNrow * (refNcol / 2 + 1);
    }
    else {                    // 3D
        fftSize = refNrow * refNcol * (refNsec / 2 + 1);
    }

    mskArea      = 0;
    bMskSingle   = false;
    bMskComputed = false;

    tplData = NULL;
    mskData = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    mskDataPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT  = NULL;
    refDataPaddedFFT  = NULL;
    refDataPaddedFFT2 = NULL;
    mskDataPaddedFFT  = NULL;
    resDataPaddedFFT  = NULL;

    fftTmp1 = NULL;
    fftTmp2 = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
    mskDataRotInit = NULL;
}

cNCorrSingle::~cNCorrSingle()
{
}

void cNCorrSingle::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    array_new_zero(resDataAbsMaxPadded, refSize);
    array_new     (resDataMaxIndPadded, refSize);

    mutexFFTW.lock();

    // complex data
    tplDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT2 = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    mskDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);

    // FFTW: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        fftPlanTpl  = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef  = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef2 = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT2,
                        FFTW_ESTIMATE);
        fftPlanMsk  = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        mskDataPadded, mskDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv1 = fftwf_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
        fftPlanInv2 = fftwf_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, fftTmp1,
                        FFTW_ESTIMATE);
        fftPlanInv3 = fftwf_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, fftTmp2,
                        FFTW_ESTIMATE);
    }
    else {                    // 3D
        // plans
        fftPlanTpl  = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef  = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef2 = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT2,
                        FFTW_ESTIMATE);
        fftPlanMsk  = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        mskDataPadded, mskDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv1 = fftwf_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
        fftPlanInv2 = fftwf_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, fftTmp1,
                        FFTW_ESTIMATE);
        fftPlanInv3 = fftwf_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, fftTmp2,
                        FFTW_ESTIMATE);
    }

    mutexFFTW.unlock();
}

void cNCorrSingle::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);
    if (mskData != NULL)  array_delete(mskData);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(mskDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(refDataPaddedFFT2);
    fftwf_free(resDataPaddedFFT);

    // FFTW
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanRef2);
    fftwf_destroy_plan(fftPlanMsk);
    fftwf_destroy_plan(fftPlanInv1);
    fftwf_destroy_plan(fftPlanInv2);
    fftwf_destroy_plan(fftPlanInv3);

    mutexFFTW.unlock();

    // 3D
    if (tplDataRotInit != NULL) array_delete(tplDataRotInit);
    if (mskDataRotInit != NULL) array_delete(mskDataRotInit);
}

void cNCorrSingle::prepareRef(float *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftwf_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftwf_execute(fftPlanRef2);
}

void cNCorrSingle::setSizeTpl(size_t tplNrowInput,
                              size_t tplNcolInput,
                              size_t tplNsecInput)
{
    assert(tplNrowInput > 0 && tplNcolInput > 0 && tplNsecInput > 0);
    assert(!(refNsec == 1 && tplNsecInput != 1));
    assert(!(refNsec != 1 && tplNsecInput == 1));

    tplNrow = tplNrowInput;
    tplNcol = tplNcolInput;
    tplNsec = tplNsecInput;
    tplSize = tplNrow * tplNcol * tplNsec;

    if (tplData == NULL && mskData == NULL) {
        array_new(tplData, tplSize);
        array_new(mskData, tplSize);
    }
}

void cNCorrSingle::copyTplMsk(float *tplDataInput, float *mskDataInput)
{
    array_memcpy(tplData, tplDataInput, tplSize);
    array_memcpy(mskData, mskDataInput, tplSize);
}

void cNCorrSingle::normalizeTpl(size_t mskAreaInput)
{
    if (mskAreaInput > 0) {
        mskArea = mskAreaInput;
    }
    else {
        mskArea = xcorrNormTemplateNCC(tplData, tplSize, mskData);
    }
}

void cNCorrSingle::prepareTplMsk(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftwf_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(mskData,       tplNrow, tplNcol,
                      mskDataPadded, refNrow, refNcol,
                                     0,       0);
            fftwf_execute(fftPlanMsk);
        }
    }
    else {                    // 3D
        // template
        array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                  tplDataPadded, refNrow, refNcol, refNsec,
                                 0,       0,       0);
        fftwf_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(mskData,       tplNrow, tplNcol, tplNsec,
                      mskDataPadded, refNrow, refNcol, refNsec,
                                     0,       0,       0);
            fftwf_execute(fftPlanMsk);
        }
    }
}

void cNCorrSingle::computeCorr(void)
{
    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv1);

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
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

        // lock the next iteration
        bMskComputed = true;
    }

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2,
                          mskArea, refSize);
}

float cNCorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

float* cNCorrSingle::getResultRef(void)
{
    return resDataPadded;
}

void cNCorrSingle::getResult(float *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cNCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cNCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                                     size_t *resDataMaxIndPaddedGlobal,
                                     eXCorrMerge bAbs)
{
    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cNCorrSingle::copyTplMskRot3D(float *tplDataInput, float *mskDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplDataRotInit == NULL);
    assert(mskDataRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_new(mskDataRotInit, tplSize);
    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
    array_memcpy(mskDataRotInit, mskDataInput, tplSize);
}

void cNCorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cNCorrSingle::mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(mskDataRotInit != NULL && mskData != NULL);

    transform_rotate(mskDataRotInit, mskData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
    array_threshold(mskData, tplSize, 0.5f);
}

cNCorrDouble::cNCorrDouble(size_t refNrowInput,
                           size_t refNcolInput,
                           size_t refNsecInput)
{
    assert(refNrowInput > 0 && refNcolInput > 0 && refNsecInput > 0);

    refNrow = refNrowInput;
    refNcol = refNcolInput;
    refNsec = refNsecInput;
    refSize = refNrow * refNcol * refNsec;

    if (refNsec == 1) {       // 2D
        fftSize = refNrow * (refNcol / 2 + 1);
    }
    else {                    // 3D
        fftSize = refNrow * refNcol * (refNsec / 2 + 1);
    }

    mskArea      = 0;
    bMskSingle   = false;
    bMskComputed = false;

    tplData = NULL;
    mskData = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    mskDataPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT  = NULL;
    refDataPaddedFFT  = NULL;
    refDataPaddedFFT2 = NULL;
    mskDataPaddedFFT  = NULL;
    resDataPaddedFFT  = NULL;

    fftTmp1 = NULL;
    fftTmp2 = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
    mskDataRotInit = NULL;
}

cNCorrDouble::~cNCorrDouble()
{
}

void cNCorrDouble::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(mskDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new(fftTmp1, refSize);
    array_new(fftTmp2, refSize);

    array_new_zero(resDataAbsMaxPadded, refSize);
    array_new     (resDataMaxIndPadded, refSize);

    mutexFFTW.lock();

    // complex data
    tplDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT2 = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    mskDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);

    // FFTW: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        fftPlanTpl  = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef  = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef2 = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT2,
                        FFTW_ESTIMATE);
        fftPlanMsk  = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        mskDataPadded, mskDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv1 = fftw_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
        fftPlanInv2 = fftw_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, fftTmp1,
                        FFTW_ESTIMATE);
        fftPlanInv3 = fftw_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, fftTmp2,
                        FFTW_ESTIMATE);
    }
    else {                    // 3D
        // plans
        fftPlanTpl  = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef  = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef2 = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT2,
                        FFTW_ESTIMATE);
        fftPlanMsk  = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        mskDataPadded, mskDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv1 = fftw_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
        fftPlanInv2 = fftw_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, fftTmp1,
                        FFTW_ESTIMATE);
        fftPlanInv3 = fftw_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, fftTmp2,
                        FFTW_ESTIMATE);
    }

    mutexFFTW.unlock();
}

void cNCorrDouble::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);
    if (mskData != NULL)  array_delete(mskData);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(mskDataPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftw_free(tplDataPaddedFFT);
    fftw_free(mskDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(refDataPaddedFFT2);
    fftw_free(resDataPaddedFFT);

    // FFTW
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanRef2);
    fftw_destroy_plan(fftPlanMsk);
    fftw_destroy_plan(fftPlanInv1);
    fftw_destroy_plan(fftPlanInv2);
    fftw_destroy_plan(fftPlanInv3);

    mutexFFTW.unlock();

    // 3D
    if (tplDataRotInit != NULL) array_delete(tplDataRotInit);
    if (mskDataRotInit != NULL) array_delete(mskDataRotInit);
}

void cNCorrDouble::prepareRef(double *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftw_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftw_execute(fftPlanRef2);
}

void cNCorrDouble::setSizeTpl(size_t tplNrowInput,
                              size_t tplNcolInput,
                              size_t tplNsecInput)
{
    assert(tplNrowInput > 0 && tplNcolInput > 0 && tplNsecInput > 0);
    assert(!(refNsec == 1 && tplNsecInput != 1));
    assert(!(refNsec != 1 && tplNsecInput == 1));

    tplNrow = tplNrowInput;
    tplNcol = tplNcolInput;
    tplNsec = tplNsecInput;
    tplSize = tplNrow * tplNcol * tplNsec;

    if (tplData == NULL && mskData == NULL) {
        array_new(tplData, tplSize);
        array_new(mskData, tplSize);
    }
}

void cNCorrDouble::copyTplMsk(double *tplDataInput, double *mskDataInput)
{
    array_memcpy(tplData, tplDataInput, tplSize);
    array_memcpy(mskData, mskDataInput, tplSize);
}

void cNCorrDouble::normalizeTpl(size_t mskAreaInput)
{
    if (mskAreaInput > 0) {
        mskArea = mskAreaInput;
    }
    else {
        mskArea = xcorrNormTemplateNCC(tplData, tplSize, mskData);
    }
}

void cNCorrDouble::prepareTplMsk(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftw_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(mskData,       tplNrow, tplNcol,
                      mskDataPadded, refNrow, refNcol,
                                     0,       0);
            fftw_execute(fftPlanMsk);
        }
    }
    else {                    // 3D
        // template
        array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                  tplDataPadded, refNrow, refNcol, refNsec,
                                 0,       0,       0);
        fftw_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(mskData,       tplNrow, tplNcol, tplNsec,
                      mskDataPadded, refNrow, refNcol, refNsec,
                                     0,       0,       0);
            fftw_execute(fftPlanMsk);
        }
    }
}

void cNCorrDouble::computeCorr(void)
{
    // xcorr 1
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv1);

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
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

        // lock the next iteration
        bMskComputed = true;
    }

    // NCC
    xcorrCombineResultNCC(resDataPadded, fftTmp1, fftTmp2,
                          mskArea, refSize);
}

double cNCorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

double* cNCorrDouble::getResultRef(void)
{
    return resDataPadded;
}

void cNCorrDouble::getResult(double *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cNCorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cNCorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
                                     size_t *resDataMaxIndPaddedGlobal,
                                     eXCorrMerge bAbs)
{
    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cNCorrDouble::copyTplMskRot3D(double *tplDataInput, double *mskDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplDataRotInit == NULL);
    assert(mskDataRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_new(mskDataRotInit, tplSize);
    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
    array_memcpy(mskDataRotInit, mskDataInput, tplSize);
}

void cNCorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cNCorrDouble::mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(mskDataRotInit != NULL && mskData != NULL);

    transform_rotate(mskDataRotInit, mskData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
    array_threshold(mskData, tplSize, 0.5);
}

} // namespace gem
