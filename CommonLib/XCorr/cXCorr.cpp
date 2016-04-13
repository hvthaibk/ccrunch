/***********************************************************************
 *  File:       cXCorr.cpp
 *
 *  Purpose:    Implementation of XCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cXCorr.hpp"

namespace gem {

cXCorrSingle::cXCorrSingle(size_t refNrowInput,
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

    tplData = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT = NULL;
    refDataPaddedFFT = NULL;
    resDataPaddedFFT = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
}

cXCorrSingle::~cXCorrSingle()
{
}

void cXCorrSingle::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new_zero(resDataAbsMaxPadded, refSize);
    array_new     (resDataMaxIndPadded, refSize);

    mutexFFTW.lock();

    // complex data
    tplDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    refDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);
    resDataPaddedFFT  = (fftwf_complex*) fftwf_malloc(
                                sizeof(fftwf_complex) * fftSize);

    // FFTW: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        fftPlanTpl = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef = fftwf_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv = fftwf_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
    }
    else {                    // 3D
        // plans
        fftPlanTpl = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef = fftwf_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv = fftwf_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
    }

    mutexFFTW.unlock();
}

void cXCorrSingle::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftwf_free(tplDataPaddedFFT);
    fftwf_free(refDataPaddedFFT);
    fftwf_free(resDataPaddedFFT);

    // FFTW
    fftwf_destroy_plan(fftPlanTpl);
    fftwf_destroy_plan(fftPlanRef);
    fftwf_destroy_plan(fftPlanInv);

    mutexFFTW.unlock();

    // 3D
    if (tplDataRotInit != NULL) array_delete(tplDataRotInit);
}

void cXCorrSingle::prepareRef(float *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftwf_execute(fftPlanRef);
}

void cXCorrSingle::setSizeTpl(size_t tplNrowInput,
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

    if (tplData == NULL) {
        array_new(tplData, tplSize);
    }
}

void cXCorrSingle::copyTpl(float *tplDataInput)
{
    array_memcpy(tplData, tplDataInput, tplSize);
}

void cXCorrSingle::prepareTpl(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftwf_execute(fftPlanTpl);
    }
    else {                    // 3D
        // template
        array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                  tplDataPadded, refNrow, refNcol, refNsec,
                                 0,       0,       0);
        fftwf_execute(fftPlanTpl);
    }
}

void cXCorrSingle::computeCorr(void)
{
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (float) refSize);
    fftwf_execute(fftPlanInv);
}

float cXCorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

float* cXCorrSingle::getResultRef(void)
{
    return resDataPadded;
}

void cXCorrSingle::getResult(float *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cXCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cXCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
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

void cXCorrSingle::copyTplRot3D(float *tplDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(tplDataRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
}

void cXCorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

cXCorrDouble::cXCorrDouble(size_t refNrowInput,
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

    tplData = NULL;

    tplDataPadded = NULL;
    refDataPadded = NULL;
    resDataPadded = NULL;

    tplDataPaddedFFT = NULL;
    refDataPaddedFFT = NULL;
    resDataPaddedFFT = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    tplDataRotInit = NULL;
}

cXCorrDouble::~cXCorrDouble()
{
}

void cXCorrDouble::memAlloc(void)
{
    // memory allocation
    array_new(tplDataPadded, refSize);
    array_new(refDataPadded, refSize);
    array_new(resDataPadded, refSize);

    array_new_zero(resDataAbsMaxPadded, refSize);
    array_new     (resDataMaxIndPadded, refSize);

    mutexFFTW.lock();

    // complex data
    tplDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    refDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);
    resDataPaddedFFT  = (fftw_complex*) fftw_malloc(
                                sizeof(fftw_complex) * fftSize);

    // FFTW: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        fftPlanTpl = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef = fftw_plan_dft_r2c_2d(
                        (int) refNrow, (int) refNcol,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv = fftw_plan_dft_c2r_2d(
                        (int) refNrow, (int) refNcol,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
    }
    else {                    // 3D
        // plans
        fftPlanTpl = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        tplDataPadded, tplDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanRef = fftw_plan_dft_r2c_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        refDataPadded, refDataPaddedFFT,
                        FFTW_ESTIMATE);
        fftPlanInv = fftw_plan_dft_c2r_3d(
                        (int) refNrow, (int) refNcol, (int) refNsec,
                        resDataPaddedFFT, resDataPadded,
                        FFTW_ESTIMATE);
    }

    mutexFFTW.unlock();
}

void cXCorrDouble::memFree(void)
{
    // memory deallocation
    if (tplData != NULL)  array_delete(tplData);

    array_delete(tplDataPadded);
    array_delete(refDataPadded);
    array_delete(resDataPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    mutexFFTW.lock();

    // complex data
    fftw_free(tplDataPaddedFFT);
    fftw_free(refDataPaddedFFT);
    fftw_free(resDataPaddedFFT);

    // FFTW
    fftw_destroy_plan(fftPlanTpl);
    fftw_destroy_plan(fftPlanRef);
    fftw_destroy_plan(fftPlanInv);

    mutexFFTW.unlock();

    // 3D
    if (tplDataRotInit != NULL) array_delete(tplDataRotInit);
}

void cXCorrDouble::prepareRef(double *refDataInput)
{
    array_memcpy(refDataPadded, refDataInput, refSize);

    fftw_execute(fftPlanRef);
}

void cXCorrDouble::setSizeTpl(size_t tplNrowInput,
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

    if (tplData == NULL) {
        array_new(tplData, tplSize);
    }
}

void cXCorrDouble::copyTpl(double *tplDataInput)
{
    array_memcpy(tplData, tplDataInput, tplSize);
}

void cXCorrDouble::prepareTpl(void)
{
    if (refNsec == 1) {       // 2D
        // template
        array_pad(tplData,       tplNrow, tplNcol,
                  tplDataPadded, refNrow, refNcol,
                                 0,       0);
        fftw_execute(fftPlanTpl);
    }
    else {                    // 3D
        // template
        array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                  tplDataPadded, refNrow, refNcol, refNsec,
                                 0,       0,       0);
        fftw_execute(fftPlanTpl);
    }
}

void cXCorrDouble::computeCorr(void)
{
    xcorrModulateAndNormalize(tplDataPaddedFFT,
                              refDataPaddedFFT,
                              resDataPaddedFFT,
                              fftSize,
                              (double) refSize);
    fftw_execute(fftPlanInv);
}

double cXCorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return array_reduce_min(resDataPadded, refSize);
    }
    else {                    // 3D
        return array_reduce_max(resDataPadded, refSize);
    }
}

double* cXCorrDouble::getResultRef(void)
{
    return resDataPadded;
}

void cXCorrDouble::getResult(double *resDataPaddedOut)
{
    array_memcpy(resDataPaddedOut, resDataPadded, refSize);
}

void cXCorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    xcorrMergeResult(resDataPadded,
                     resDataAbsMaxPadded,
                     resDataMaxIndPadded,
                     refSize,
                     indx,
                     bAbs);
}

void cXCorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
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

void cXCorrDouble::copyTplRot3D(double *tplDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(tplDataRotInit == NULL);

    array_new(tplDataRotInit, tplSize);
    array_memcpy(tplDataRotInit, tplDataInput, tplSize);
}

void cXCorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(tplDataRotInit != NULL && tplData != NULL);

    transform_rotate(tplDataRotInit, tplData,
                     tplNrow, tplNcol, tplNsec,
                     alpha, beta, gamma,
                     inter);
}

} // namespace gem
