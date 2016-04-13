/***********************************************************************
 *  File:       cWCorr.cpp
 *
 *  Purpose:    Implementation of WCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cWCorr.hpp"

namespace gem {

cWCorrSingle::cWCorrSingle(size_t refNrowInput,
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

    bMskSingle   = false;
    bMskComputed = false;

    tplData = NULL;
    mskData = NULL;
    tplWght = NULL;
    refWght = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    refWghtPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT  = NULL;
    refDataPaddedFFT  = NULL;
    refDataPaddedFFT2 = NULL;
    refWghtPaddedFFT  = NULL;
    resDataPaddedFFT  = NULL;

    fftTmp1 = NULL;
    fftTmp2 = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
    mskDataRotInit = NULL;
    tplWghtRotInit = NULL;
    refWghtRotInit = NULL;
}

cWCorrSingle::~cWCorrSingle()
{
}

void cWCorrSingle::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(refWghtPadded, refSize);
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
    refWghtPaddedFFT  = (fftwf_complex*) fftwf_malloc(
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
                        refWghtPadded, refWghtPaddedFFT,
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
                        refWghtPadded, refWghtPaddedFFT,
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

void cWCorrSingle::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);
    if (mskData != NULL)  array_delete(mskData);
    if (tplWght != NULL)  array_delete(tplWght);
    if (refWght != NULL)  array_delete(refWght);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(refWghtPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(refWghtPaddedFFT);
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
    if (tplWghtRotInit != NULL) array_delete(tplWghtRotInit);
    if (refWghtRotInit != NULL) array_delete(refWghtRotInit);
}

void cWCorrSingle::prepareRef(float *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftwf_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftwf_execute(fftPlanRef2);
}

void cWCorrSingle::setSizeTpl(size_t tplNrowInput,
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

    if (tplData == NULL && mskData == NULL && tplWght == NULL && refWght == NULL) {
        array_new(tplData, tplSize);
        array_new(mskData, tplSize);
        array_new(tplWght, tplSize);
        array_new(refWght, tplSize);
    }
}

void cWCorrSingle::copyTplMskWght(float *tplDataInput,
                                  float *mskDataInput,
                                  float *tplWghtInput,
                                  float *refWghtInput)
{
    // template
    array_memcpy(tplData, tplDataInput, tplSize);

    // mask and weights
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        array_memcpy(mskData, mskDataInput, tplSize);
        array_memcpy(tplWght, tplWghtInput, tplSize);
        array_memcpy(refWght, refWghtInput, tplSize);
    }
}

void cWCorrSingle::normalizeTpl()
{
    xcorrNormTemplateWCC2(tplData, tplSize,
                          mskData,
                          tplWght,
                          refWght);
}

void cWCorrSingle::prepareTplRefWght(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftwf_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(refWght,       tplNrow, tplNcol,
                      refWghtPadded, refNrow, refNcol,
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
            array_pad(refWght,       tplNrow, tplNcol, tplNsec,
                      refWghtPadded, refNrow, refNcol, refNsec,
                                     0,       0,       0);
            fftwf_execute(fftPlanMsk);
        }
    }
}

void cWCorrSingle::computeCorr(void)
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
        xcorrModulateAndNormalize(refWghtPaddedFFT,
                                  refDataPaddedFFT,
                                  resDataPaddedFFT,
                                  fftSize,
                                  (float) refSize);
        fftwf_execute(fftPlanInv2);

        // xcorr 3
        xcorrModulateAndNormalize(refWghtPaddedFFT,
                                  refDataPaddedFFT2,
                                  resDataPaddedFFT,
                                  fftSize,
                                  (float) refSize);
        fftwf_execute(fftPlanInv3);

        // lock the next iteration
        bMskComputed = true;
    }

    // WCC
    xcorrCombineResultWCC(resDataPadded,
                          fftTmp1,
                          fftTmp2,
                          array_reduce_sum(tplData, tplSize),
                          refSize);
}

float cWCorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

float* cWCorrSingle::getResultRef(void)
{
    return resDataPadded;
}

void cWCorrSingle::getResult(float *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cWCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cWCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
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

void cWCorrSingle::copyTplMskWghtRot3D(float *tplDataInput,
                                       float *mskDataInput,
                                       float *tplWghtInput,
                                       float *refWghtInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplWghtInput != NULL);
    assert(refWghtInput != NULL);
    assert(tplDataRotInit == NULL);
    assert(mskDataRotInit == NULL);
    assert(tplWghtRotInit == NULL);
    assert(refWghtRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_new(mskDataRotInit, tplSize);
    array_new(tplWghtRotInit, tplSize);
    array_new(refWghtRotInit, tplSize);

    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
    array_memcpy(mskDataRotInit, mskDataInput, tplSize);
    array_memcpy(tplWghtRotInit, tplWghtInput, tplSize);
    array_memcpy(refWghtRotInit, refWghtInput, tplSize);
}

void cWCorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cWCorrSingle::mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(mskDataRotInit != NULL && mskData != NULL);

    transform_rotate(mskDataRotInit, mskData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
    array_threshold(mskData, tplSize, 0.5f);
}

void cWCorrSingle::wghtTplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplWghtRotInit != NULL && tplWght != NULL);

    transform_rotate(tplWghtRotInit, tplWght,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cWCorrSingle::wghtRefCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(refWghtRotInit != NULL && refWght != NULL);

    transform_rotate(refWghtRotInit, refWght,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

cWCorrDouble::cWCorrDouble(size_t refNrowInput,
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

    bMskSingle   = false;
    bMskComputed = false;

    tplData = NULL;
    mskData = NULL;
    tplWght = NULL;
    refWght = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    refWghtPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT  = NULL;
    refDataPaddedFFT  = NULL;
    refDataPaddedFFT2 = NULL;
    refWghtPaddedFFT  = NULL;
    resDataPaddedFFT  = NULL;

    fftTmp1 = NULL;
    fftTmp2 = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
    mskDataRotInit = NULL;
    tplWghtRotInit = NULL;
    refWghtRotInit = NULL;
}

cWCorrDouble::~cWCorrDouble()
{
}

void cWCorrDouble::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(refWghtPadded, refSize);
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
    refWghtPaddedFFT  = (fftw_complex*) fftw_malloc(
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
                        refWghtPadded, refWghtPaddedFFT,
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
                        refWghtPadded, refWghtPaddedFFT,
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

void cWCorrDouble::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);
    if (mskData != NULL)  array_delete(mskData);
    if (tplWght != NULL)  array_delete(tplWght);
    if (refWght != NULL)  array_delete(refWght);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(refWghtPadded);
    array_delete(resDataPadded);

    array_delete(fftTmp1);
    array_delete(fftTmp2);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftw_free(tplDataPaddedFFT);
    fftw_free(refWghtPaddedFFT);
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
    if (tplWghtRotInit != NULL) array_delete(tplWghtRotInit);
    if (refWghtRotInit != NULL) array_delete(refWghtRotInit);
}

void cWCorrDouble::prepareRef(double *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftw_execute(fftPlanRef);
    array_math_sqr(refDataPadded, refSize);
    fftw_execute(fftPlanRef2);
}

void cWCorrDouble::setSizeTpl(size_t tplNrowInput,
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

    if (tplData == NULL && mskData == NULL && tplWght == NULL && refWght == NULL) {
        array_new(tplData, tplSize);
        array_new(mskData, tplSize);
        array_new(tplWght, tplSize);
        array_new(refWght, tplSize);
    }
}

void cWCorrDouble::copyTplMskWght(double *tplDataInput,
                                  double *mskDataInput,
                                  double *tplWghtInput,
                                  double *refWghtInput)
{
    // template
    array_memcpy(tplData, tplDataInput, tplSize);

    // mask and weights
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        array_memcpy(mskData, mskDataInput, tplSize);
        array_memcpy(tplWght, tplWghtInput, tplSize);
        array_memcpy(refWght, refWghtInput, tplSize);
    }
}

void cWCorrDouble::normalizeTpl()
{
    xcorrNormTemplateWCC2(tplData, tplSize,
                          mskData,
                          tplWght,
                          refWght);
}

void cWCorrDouble::prepareTplRefWght(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftw_execute(fftPlanTpl);

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            array_pad(refWght,       tplNrow, tplNcol,
                      refWghtPadded, refNrow, refNcol,
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
            array_pad(refWght,       tplNrow, tplNcol, tplNsec,
                      refWghtPadded, refNrow, refNcol, refNsec,
                                     0,       0,       0);
            fftw_execute(fftPlanMsk);
        }
    }
}

void cWCorrDouble::computeCorr(void)
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
        xcorrModulateAndNormalize(refWghtPaddedFFT,
                                  refDataPaddedFFT,
                                  resDataPaddedFFT,
                                  fftSize,
                                  (double) refSize);
        fftw_execute(fftPlanInv2);

        // xcorr 3
        xcorrModulateAndNormalize(refWghtPaddedFFT,
                                  refDataPaddedFFT2,
                                  resDataPaddedFFT,
                                  fftSize,
                                  (double) refSize);
        fftw_execute(fftPlanInv3);

        // lock the next iteration
        bMskComputed = true;
    }

    // WCC
    xcorrCombineResultWCC(resDataPadded,
                          fftTmp1,
                          fftTmp2,
                          array_reduce_sum(tplData, tplSize),
                          refSize);
}

double cWCorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

double* cWCorrDouble::getResultRef(void)
{
    return resDataPadded;
}

void cWCorrDouble::getResult(double *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cWCorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cWCorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
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

void cWCorrDouble::copyTplMskWghtRot3D(double *tplDataInput,
                                       double *mskDataInput,
                                       double *tplWghtInput,
                                       double *refWghtInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplWghtInput != NULL);
    assert(refWghtInput != NULL);
    assert(tplDataRotInit == NULL);
    assert(mskDataRotInit == NULL);
    assert(tplWghtRotInit == NULL);
    assert(refWghtRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_new(mskDataRotInit, tplSize);
    array_new(tplWghtRotInit, tplSize);
    array_new(refWghtRotInit, tplSize);

    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
    array_memcpy(mskDataRotInit, mskDataInput, tplSize);
    array_memcpy(tplWghtRotInit, tplWghtInput, tplSize);
    array_memcpy(refWghtRotInit, refWghtInput, tplSize);
}

void cWCorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cWCorrDouble::mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(mskDataRotInit != NULL && mskData != NULL);

    transform_rotate(mskDataRotInit, mskData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
    array_threshold(mskData, tplSize, 0.5);
}

void cWCorrDouble::wghtTplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplWghtRotInit != NULL && tplWght != NULL);

    transform_rotate(tplWghtRotInit, tplWght,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

void cWCorrDouble::wghtRefCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(refWghtRotInit != NULL && refWght != NULL);

    transform_rotate(refWghtRotInit, refWght,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

} // namespace gem
