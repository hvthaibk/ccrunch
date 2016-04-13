/***********************************************************************
 *  File:       cuWCorr.cu
 *
 *  Purpose:    Implementation of WCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuWCorr.cuh"

namespace gem {

cuWCorrSingle::cuWCorrSingle(size_t refNrowInput,
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

    dev_tplData = NULL;
    dev_mskData = NULL;
    dev_tplWght = NULL;
    dev_refWght = NULL;

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_refWghtPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT2 = NULL;
    dev_refWghtPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_fftTmp1 = NULL;
    dev_fftTmp2 = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    bRotTexture = false;
    dev_tplDataRotInit = NULL;
    dev_mskDataRotInit = NULL;
    dev_tplWghtRotInit = NULL;
    dev_refWghtRotInit = NULL;
    dev_tplDataRotInit_tex3D = NULL;
    dev_mskDataRotInit_tex3D = NULL;
    dev_tplWghtRotInit_tex3D = NULL;
    dev_refWghtRotInit_tex3D = NULL;
}

cuWCorrSingle::~cuWCorrSingle()
{
}

void cuWCorrSingle::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_refWghtPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new(dev_fftTmp1, refSize);
    cuda_arrayDev_new(dev_fftTmp2, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_refWghtPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT , fftSize);

    // cuFFT: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        CUFFT_SAFE_CALL(cufftPlan2d(&fftPlanFwd,
                                    (int) refNrow,
                                    (int) refNcol,
                                    CUFFT_R2C));
        CUFFT_SAFE_CALL(cufftPlan2d(&fftPlanInv,
                                    (int) refNrow,
                                    (int) refNcol,
                                    CUFFT_C2R));
    }
    else {                    // 3D
        // plans
        CUFFT_SAFE_CALL(cufftPlan3d(&fftPlanFwd,
                                    (int) refNrow,
                                    (int) refNcol,
                                    (int) refNsec,
                                    CUFFT_R2C));
        CUFFT_SAFE_CALL(cufftPlan3d(&fftPlanInv,
                                    (int) refNrow,
                                    (int) refNcol,
                                    (int) refNsec,
                                    CUFFT_C2R));
    }
}

void cuWCorrSingle::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);
    if (dev_mskData != NULL)  cuda_arrayDev_delete(dev_mskData);
    if (dev_tplWght != NULL)  cuda_arrayDev_delete(dev_tplWght);
    if (dev_refWght != NULL)  cuda_arrayDev_delete(dev_refWght);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_refWghtPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_refWghtPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);

    // cuFFT
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanInv));

    // 3D
    if (bRotTexture) {
        if (dev_tplDataRotInit_tex3D != NULL)  cuda_arrayArr_delete(dev_tplDataRotInit_tex3D);
        if (dev_mskDataRotInit_tex3D != NULL)  cuda_arrayArr_delete(dev_mskDataRotInit_tex3D);
        if (dev_tplWghtRotInit_tex3D != NULL)  cuda_arrayArr_delete(dev_tplWghtRotInit_tex3D);
        if (dev_refWghtRotInit_tex3D != NULL)  cuda_arrayArr_delete(dev_refWghtRotInit_tex3D);
    }
    else {
        if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
        if (dev_mskDataRotInit != NULL)  cuda_arrayDev_delete(dev_mskDataRotInit);
        if (dev_tplWghtRotInit != NULL)  cuda_arrayDev_delete(dev_tplWghtRotInit);
        if (dev_refWghtRotInit != NULL)  cuda_arrayDev_delete(dev_refWghtRotInit);
    }
}

void cuWCorrSingle::prepareRef(float *refDataInput)
{
    cuda_array_memcpy_h2d(dev_refDataPadded, refDataInput, refSize);

    CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                 (cufftReal *)    dev_refDataPadded,
                                 (cufftComplex *) dev_refDataPaddedFFT));

    cuda_array_math_sqr(dev_refDataPadded, refSize);

    CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                 (cufftReal *)    dev_refDataPadded,
                                 (cufftComplex *) dev_refDataPaddedFFT2));
}

void cuWCorrSingle::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL && dev_mskData == NULL &&
        dev_tplWght == NULL && dev_refWght == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
        cuda_arrayDev_new(dev_mskData, tplSize);
        cuda_arrayDev_new(dev_tplWght, tplSize);
        cuda_arrayDev_new(dev_refWght, tplSize);
    }
}

void cuWCorrSingle::copyTplMskWght(float *tplDataInput,
                                   float *mskDataInput,
                                   float *tplWghtInput,
                                   float *refWghtInput)
{
    // template
    cuda_array_memcpy_h2d(dev_tplData, tplDataInput, tplSize);

    // mask and weights
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        cuda_array_memcpy_h2d(dev_mskData, mskDataInput, tplSize);
        cuda_array_memcpy_h2d(dev_tplWght, tplWghtInput, tplSize);
        cuda_array_memcpy_h2d(dev_refWght, refWghtInput, tplSize);
    }
}

void cuWCorrSingle::normalizeTpl()
{
    cuda_xcorrNormTemplateWCC2(dev_tplData, tplSize,
                               dev_mskData,
                               dev_tplWght,
                               dev_refWght);
}

void cuWCorrSingle::prepareTplRefWght()
{
    if (refNsec == 1) {       // 2D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol,
                       dev_tplDataPadded, refNrow, refNcol,
                                          0,       0);
        CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                     (cufftReal *)    dev_tplDataPadded,
                                     (cufftComplex *) dev_tplDataPaddedFFT));

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            cuda_array_pad(dev_refWght,       tplNrow, tplNcol,
                           dev_refWghtPadded, refNrow, refNcol,
                                              0,       0);
            CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                         (cufftReal *)    dev_refWghtPadded,
                                         (cufftComplex *) dev_refWghtPaddedFFT));
        }
    }
    else {                    // 3D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol, tplNsec,
                       dev_tplDataPadded, refNrow, refNcol, refNsec,
                                          0,       0,       0);
        CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                     (cufftReal *)    dev_tplDataPadded,
                                     (cufftComplex *) dev_tplDataPaddedFFT));

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            cuda_array_pad(dev_refWght,       tplNrow, tplNcol, tplNsec,
                           dev_refWghtPadded, refNrow, refNcol, refNsec,
                                              0,       0,       0);
            CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                         (cufftReal *)    dev_refWghtPadded,
                                         (cufftComplex *) dev_refWghtPaddedFFT));
        }
    }
}

void cuWCorrSingle::computeCorr()
{
    // xcorr 1
    cuda_xcorrModulateAndNormalize
        (dev_tplDataPaddedFFT,
         dev_refDataPaddedFFT,
         dev_resDataPaddedFFT,
         fftSize,
         (float) refSize);
    CUFFT_SAFE_CALL(cufftExecC2R(fftPlanInv,
                                 (cufftComplex *) dev_resDataPaddedFFT,
                                 (cufftReal *)    dev_resDataPadded));

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        // xcorr 2
        cuda_xcorrModulateAndNormalize
            (dev_refWghtPaddedFFT,
             dev_refDataPaddedFFT,
             dev_resDataPaddedFFT,
             fftSize,
             (float) refSize);
        CUFFT_SAFE_CALL(cufftExecC2R(fftPlanInv,
                                     (cufftComplex *) dev_resDataPaddedFFT,
                                     (cufftReal *)    dev_fftTmp1));

        // xcorr 3
        cuda_xcorrModulateAndNormalize
            (dev_refWghtPaddedFFT,
             dev_refDataPaddedFFT2,
             dev_resDataPaddedFFT,
             fftSize,
             (float) refSize);
        CUFFT_SAFE_CALL(cufftExecC2R(fftPlanInv,
                                     (cufftComplex *) dev_resDataPaddedFFT,
                                     (cufftReal *)    dev_fftTmp2));

        // lock the next iteration
        bMskComputed = true;
    }

    // WCC
    cuda_xcorrCombineResultWCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               cuda_array_reduce_sum(dev_tplData, tplSize),
                               refSize);
}

float cuWCorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

float* cuWCorrSingle::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuWCorrSingle::getResult(float *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuWCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuWCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
                                      size_t *resDataMaxIndPaddedGlobal,
                                      eXCorrMerge bAbs)
{
    cuda_array_memcpy_d2h(resDataAbsMaxPadded, dev_resDataAbsMaxPadded, refSize);
    cuda_array_memcpy_d2h(resDataMaxIndPadded, dev_resDataMaxIndPadded, refSize);

    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cuWCorrSingle::copyTplMskWghtRot3D(float *tplDataInput,
                                        float *mskDataInput,
                                        float *tplWghtInput,
                                        float *refWghtInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplWghtInput != NULL);
    assert(refWghtInput != NULL);
    assert(dev_tplDataRotInit == NULL && dev_tplDataRotInit_tex3D == NULL);
    assert(dev_mskDataRotInit == NULL && dev_mskDataRotInit_tex3D == NULL);
    assert(dev_tplWghtRotInit == NULL && dev_tplWghtRotInit_tex3D == NULL);
    assert(dev_refWghtRotInit == NULL && dev_refWghtRotInit_tex3D == NULL);

    if (bRotTexture) {
        cuda_arrayArr_new(dev_tplDataRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_arrayArr_new(dev_mskDataRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_arrayArr_new(dev_tplWghtRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_arrayArr_new(dev_refWghtRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_tplDataRotInit_tex3D, tplDataInput, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_mskDataRotInit_tex3D, mskDataInput, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_tplWghtRotInit_tex3D, tplWghtInput, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_refWghtRotInit_tex3D, refWghtInput, tplNrow, tplNcol, tplNsec);
    }
    else {
        cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
        cuda_arrayDev_new(dev_mskDataRotInit, tplSize);
        cuda_arrayDev_new(dev_tplWghtRotInit, tplSize);
        cuda_arrayDev_new(dev_refWghtRotInit, tplSize);
        cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
        cuda_array_memcpy_h2d(dev_mskDataRotInit, mskDataInput, tplSize);
        cuda_array_memcpy_h2d(dev_tplWghtRotInit, tplWghtInput, tplSize);
        cuda_array_memcpy_h2d(dev_refWghtRotInit, refWghtInput, tplSize);
    }
}

void cuWCorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplDataRotInit != NULL && dev_tplData != NULL);

    if (bRotTexture) {
        cuda_transform_rotate(dev_tplDataRotInit_tex3D, dev_tplData,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
    else {
        cuda_transform_rotate(dev_tplDataRotInit, dev_tplData,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
}

void cuWCorrSingle::mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_mskDataRotInit != NULL && dev_mskData != NULL);

    if (bRotTexture) {
        cuda_transform_rotate(dev_mskDataRotInit_tex3D, dev_mskData,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
    else {
        cuda_transform_rotate(dev_mskDataRotInit, dev_mskData,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
    cuda_array_threshold(dev_mskData, tplSize, 0.5f);
}

void cuWCorrSingle::wghtTplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplWghtRotInit != NULL && dev_tplWght != NULL);

    if (bRotTexture) {
        cuda_transform_rotate(dev_tplWghtRotInit_tex3D, dev_tplWght,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
    else {
        cuda_transform_rotate(dev_tplWghtRotInit, dev_tplWght,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
}

void cuWCorrSingle::wghtRefCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_refWghtRotInit != NULL && dev_refWght != NULL);

    if (bRotTexture) {
        cuda_transform_rotate(dev_refWghtRotInit_tex3D, dev_refWght,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
    else {
        cuda_transform_rotate(dev_refWghtRotInit, dev_refWght,
                              tplNrow, tplNcol, tplNsec,
                              alpha, beta, gamma,
                              inter);
    }
}

cuWCorrDouble::cuWCorrDouble(size_t refNrowInput,
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

    dev_tplData = NULL;
    dev_mskData = NULL;
    dev_tplWght = NULL;
    dev_refWght = NULL;

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_refWghtPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT2 = NULL;
    dev_refWghtPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_fftTmp1 = NULL;
    dev_fftTmp2 = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    dev_tplDataRotInit = NULL;
    dev_mskDataRotInit = NULL;
    dev_tplWghtRotInit = NULL;
    dev_refWghtRotInit = NULL;
}

cuWCorrDouble::~cuWCorrDouble()
{
}

void cuWCorrDouble::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_refWghtPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new(dev_fftTmp1, refSize);
    cuda_arrayDev_new(dev_fftTmp2, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_refWghtPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT , fftSize);

    // cuFFT: single-precision
    if (refNsec == 1) {       // 2D
        // plans
        CUFFT_SAFE_CALL(cufftPlan2d(&fftPlanFwd,
                                    (int) refNrow,
                                    (int) refNcol,
                                    CUFFT_D2Z));
        CUFFT_SAFE_CALL(cufftPlan2d(&fftPlanInv,
                                    (int) refNrow,
                                    (int) refNcol,
                                    CUFFT_Z2D));
    }
    else {                    // 3D
        // plans
        CUFFT_SAFE_CALL(cufftPlan3d(&fftPlanFwd,
                                    (int) refNrow,
                                    (int) refNcol,
                                    (int) refNsec,
                                    CUFFT_D2Z));
        CUFFT_SAFE_CALL(cufftPlan3d(&fftPlanInv,
                                    (int) refNrow,
                                    (int) refNcol,
                                    (int) refNsec,
                                    CUFFT_Z2D));
    }
}

void cuWCorrDouble::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);
    if (dev_mskData != NULL)  cuda_arrayDev_delete(dev_mskData);
    if (dev_tplWght != NULL)  cuda_arrayDev_delete(dev_tplWght);
    if (dev_refWght != NULL)  cuda_arrayDev_delete(dev_refWght);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_refWghtPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_refWghtPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);

    // cuFFT
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanInv));

    // 3D
    if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
    if (dev_mskDataRotInit != NULL)  cuda_arrayDev_delete(dev_mskDataRotInit);
    if (dev_tplWghtRotInit != NULL)  cuda_arrayDev_delete(dev_tplWghtRotInit);
    if (dev_refWghtRotInit != NULL)  cuda_arrayDev_delete(dev_refWghtRotInit);
}

void cuWCorrDouble::prepareRef(double *refDataInput)
{
    cuda_array_memcpy_h2d(dev_refDataPadded, refDataInput, refSize);

    CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                 (cufftDoubleReal *)    dev_refDataPadded,
                                 (cufftDoubleComplex *) dev_refDataPaddedFFT));

    cuda_array_math_sqr(dev_refDataPadded, refSize);

    CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                 (cufftDoubleReal *)    dev_refDataPadded,
                                 (cufftDoubleComplex *) dev_refDataPaddedFFT2));
}

void cuWCorrDouble::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL && dev_mskData == NULL &&
        dev_tplWght == NULL && dev_refWght == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
        cuda_arrayDev_new(dev_mskData, tplSize);
        cuda_arrayDev_new(dev_tplWght, tplSize);
        cuda_arrayDev_new(dev_refWght, tplSize);
    }
}

void cuWCorrDouble::copyTplMskWght(double *tplDataInput,
                                   double *mskDataInput,
                                   double *tplWghtInput,
                                   double *refWghtInput)
{
    // template
    cuda_array_memcpy_h2d(dev_tplData, tplDataInput, tplSize);

    // mask and weights
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        cuda_array_memcpy_h2d(dev_mskData, mskDataInput, tplSize);
        cuda_array_memcpy_h2d(dev_tplWght, tplWghtInput, tplSize);
        cuda_array_memcpy_h2d(dev_refWght, refWghtInput, tplSize);
    }
}

void cuWCorrDouble::normalizeTpl()
{
    cuda_xcorrNormTemplateWCC2(dev_tplData, tplSize,
                               dev_mskData,
                               dev_tplWght,
                               dev_refWght);
}

void cuWCorrDouble::prepareTplRefWght()
{
    if (refNsec == 1) {       // 2D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol,
                       dev_tplDataPadded, refNrow, refNcol,
                                          0,       0);
        CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                     (cufftDoubleReal *)    dev_tplDataPadded,
                                     (cufftDoubleComplex *) dev_tplDataPaddedFFT));

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            cuda_array_pad(dev_refWght,       tplNrow, tplNcol,
                           dev_refWghtPadded, refNrow, refNcol,
                                              0,       0);
            CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                         (cufftDoubleReal *)    dev_refWghtPadded,
                                         (cufftDoubleComplex *) dev_refWghtPaddedFFT));
        }
    }
    else {                    // 3D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol, tplNsec,
                       dev_tplDataPadded, refNrow, refNcol, refNsec,
                                          0,       0,       0);
        CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                     (cufftDoubleReal *)    dev_tplDataPadded,
                                     (cufftDoubleComplex *) dev_tplDataPaddedFFT));

        // mask
        if ((bMskSingle && !bMskComputed) || !bMskSingle) {
            cuda_array_pad(dev_refWght,       tplNrow, tplNcol, tplNsec,
                           dev_refWghtPadded, refNrow, refNcol, refNsec,
                                              0,       0,       0);
            CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                         (cufftDoubleReal *)    dev_refWghtPadded,
                                         (cufftDoubleComplex *) dev_refWghtPaddedFFT));
        }
    }
}

void cuWCorrDouble::computeCorr()
{
    // xcorr 1
    cuda_xcorrModulateAndNormalize
        (dev_tplDataPaddedFFT,
         dev_refDataPaddedFFT,
         dev_resDataPaddedFFT,
         fftSize,
         (double) refSize);
    CUFFT_SAFE_CALL(cufftExecZ2D(fftPlanInv,
                                 (cufftDoubleComplex *) dev_resDataPaddedFFT,
                                 (cufftDoubleReal *)    dev_resDataPadded));

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
        // xcorr 2
        cuda_xcorrModulateAndNormalize
            (dev_refWghtPaddedFFT,
             dev_refDataPaddedFFT,
             dev_resDataPaddedFFT,
             fftSize,
             (double) refSize);
        CUFFT_SAFE_CALL(cufftExecZ2D(fftPlanInv,
                                     (cufftDoubleComplex *) dev_resDataPaddedFFT,
                                     (cufftDoubleReal *)    dev_fftTmp1));

        // xcorr 3
        cuda_xcorrModulateAndNormalize
            (dev_refWghtPaddedFFT,
             dev_refDataPaddedFFT2,
             dev_resDataPaddedFFT,
             fftSize,
             (double) refSize);
        CUFFT_SAFE_CALL(cufftExecZ2D(fftPlanInv,
                                     (cufftDoubleComplex *) dev_resDataPaddedFFT,
                                     (cufftDoubleReal *)    dev_fftTmp2));

        // lock the next iteration
        bMskComputed = true;
    }

    // WCC
    cuda_xcorrCombineResultWCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               cuda_array_reduce_sum(dev_tplData, tplSize),
                               refSize);
}

double cuWCorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

double* cuWCorrDouble::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuWCorrDouble::getResult(double *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuWCorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuWCorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
                                      size_t *resDataMaxIndPaddedGlobal,
                                      eXCorrMerge bAbs)
{
    cuda_array_memcpy_d2h(resDataAbsMaxPadded, dev_resDataAbsMaxPadded, refSize);
    cuda_array_memcpy_d2h(resDataMaxIndPadded, dev_resDataMaxIndPadded, refSize);

    xcorrMergeResultGlobal(resDataAbsMaxPadded,
                           resDataMaxIndPadded,
                           resDataAbsMaxPaddedGlobal,
                           resDataMaxIndPaddedGlobal,
                           refSize,
                           bAbs);
}

void cuWCorrDouble::copyTplMskWghtRot3D(double *tplDataInput,
                                        double *mskDataInput,
                                        double *tplWghtInput,
                                        double *refWghtInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(tplWghtInput != NULL);
    assert(refWghtInput != NULL);
    assert(dev_tplDataRotInit == NULL);
    assert(dev_mskDataRotInit == NULL);
    assert(dev_tplWghtRotInit == NULL);
    assert(dev_refWghtRotInit == NULL);

    cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
    cuda_arrayDev_new(dev_mskDataRotInit, tplSize);
    cuda_arrayDev_new(dev_tplWghtRotInit, tplSize);
    cuda_arrayDev_new(dev_refWghtRotInit, tplSize);
    cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
    cuda_array_memcpy_h2d(dev_mskDataRotInit, mskDataInput, tplSize);
    cuda_array_memcpy_h2d(dev_tplWghtRotInit, tplWghtInput, tplSize);
    cuda_array_memcpy_h2d(dev_refWghtRotInit, refWghtInput, tplSize);
}

void cuWCorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplDataRotInit != NULL && dev_tplData != NULL);

    cuda_transform_rotate(dev_tplDataRotInit, dev_tplData,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
}

void cuWCorrDouble::mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_mskDataRotInit != NULL && dev_mskData != NULL);

    cuda_transform_rotate(dev_mskDataRotInit, dev_mskData,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
    cuda_array_threshold(dev_mskData, tplSize, 0.5);
}

void cuWCorrDouble::wghtTplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplWghtRotInit != NULL && dev_tplWght != NULL);

    cuda_transform_rotate(dev_tplWghtRotInit, dev_tplWght,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
}

void cuWCorrDouble::wghtRefCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_refWghtRotInit != NULL && dev_refWght != NULL);

    cuda_transform_rotate(dev_refWghtRotInit, dev_refWght,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
}

} // namespace gem
