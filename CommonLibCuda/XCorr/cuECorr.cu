/***********************************************************************
 *  File:       cuECorr.cu
 *
 *  Purpose:    Implementation of ECC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuECorr.cuh"

namespace gem {

cuECorrSingle::cuECorrSingle(size_t refNrowInput,
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

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_mskDataPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT2 = NULL;
    dev_mskDataPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_fftTmp1 = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    bRotTexture = false;
    dev_tplDataRotInit = NULL;
    dev_mskDataRotInit = NULL;
    dev_tplDataRotInit_tex3D = NULL;
    dev_mskDataRotInit_tex3D = NULL;
}

cuECorrSingle::~cuECorrSingle()
{
}

void cuECorrSingle::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_mskDataPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new(dev_fftTmp1, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT , fftSize);
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

void cuECorrSingle::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);
    if (dev_mskData != NULL)  cuda_arrayDev_delete(dev_mskData);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_fftTmp1);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
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
    }
    else {
        if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
        if (dev_mskDataRotInit != NULL)  cuda_arrayDev_delete(dev_mskDataRotInit);
    }
}

void cuECorrSingle::prepareRef(float *refDataInput)
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

void cuECorrSingle::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL && dev_mskData == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
        cuda_arrayDev_new(dev_mskData, tplSize);
    }
}

void cuECorrSingle::copyTplMsk(float *tplDataInput, float *mskDataInput)
{
    // template
#ifdef __GEM_USE_CUDA_PINNEDMEM__
    CUDA_SAFE_CALL( cudaMemcpyAsync(dev_tplData,
                                    tplDataInput,
                                    tplSize*sizeof(float),
                                    cudaMemcpyHostToDevice) );
#else
    cuda_array_memcpy_h2d(dev_tplData, tplDataInput, tplSize);
#endif

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
#ifdef __GEM_USE_CUDA_PINNEDMEM__
        CUDA_SAFE_CALL( cudaMemcpyAsync(dev_mskData,
                                        mskDataInput,
                                        tplSize*sizeof(float),
                                        cudaMemcpyHostToDevice) );
#else
        cuda_array_memcpy_h2d(dev_mskData, mskDataInput, tplSize);
#endif
    }
}

void cuECorrSingle::normalizeTpl(void)
{
    //cuda_array_math_mul(dev_tplData, dev_mskData, tplSize);
    //cuda_array_math_div(dev_tplData, std::sqrt(cuda_array_reduce_sum2(dev_tplData, tplSize)), tplSize);
    cuda_xcorrNormTemplateECC(dev_tplData, tplSize, dev_mskData);
}

void cuECorrSingle::prepareTplMsk(void)
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
            cuda_array_pad(dev_mskData,       tplNrow, tplNcol,
                           dev_mskDataPadded, refNrow, refNcol,
                                              0,       0);
            CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                         (cufftReal *)    dev_mskDataPadded,
                                         (cufftComplex *) dev_mskDataPaddedFFT));
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
            cuda_array_pad(dev_mskData,       tplNrow, tplNcol, tplNsec,
                           dev_mskDataPadded, refNrow, refNcol, refNsec,
                                              0,       0,       0);
            CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                         (cufftReal *)    dev_mskDataPadded,
                                         (cufftComplex *) dev_mskDataPaddedFFT));
        }
    }
}

void cuECorrSingle::computeCorr(void)
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
            (dev_mskDataPaddedFFT,
             dev_refDataPaddedFFT2,
             dev_resDataPaddedFFT,
             fftSize,
             (float) refSize);
        CUFFT_SAFE_CALL(cufftExecC2R(fftPlanInv,
                                     (cufftComplex *) dev_resDataPaddedFFT,
                                     (cufftReal *)    dev_fftTmp1));

        // lock the next iteration
        bMskComputed = true;
    }

    // ECC
    //cuda_array_math_sqrt(dev_fftTmp1, refSize);
    //cuda_array_math_div(dev_resDataPadded, dev_fftTmp1, refSize);
    cuda_xcorrCombineResultECC(dev_resDataPadded, dev_fftTmp1, refSize);
}

float cuECorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

float* cuECorrSingle::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuECorrSingle::getResult(float *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuECorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuECorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
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

void cuECorrSingle::copyTplMskRot3D(float *tplDataInput, float *mskDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(dev_tplDataRotInit == NULL && dev_tplDataRotInit_tex3D == NULL);
    assert(dev_mskDataRotInit == NULL && dev_mskDataRotInit_tex3D == NULL);

    if (bRotTexture) {
        cuda_arrayArr_new(dev_tplDataRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_arrayArr_new(dev_mskDataRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_tplDataRotInit_tex3D, tplDataInput, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_mskDataRotInit_tex3D, mskDataInput, tplNrow, tplNcol, tplNsec);
    }
    else {
        cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
        cuda_arrayDev_new(dev_mskDataRotInit, tplSize);
        cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
        cuda_array_memcpy_h2d(dev_mskDataRotInit, mskDataInput, tplSize);
    }
}

void cuECorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
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

void cuECorrSingle::mskCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
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

cuECorrDouble::cuECorrDouble(size_t refNrowInput,
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

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_mskDataPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT2 = NULL;
    dev_mskDataPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_fftTmp1 = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    dev_tplDataRotInit = NULL;
    dev_mskDataRotInit = NULL;
}

cuECorrDouble::~cuECorrDouble()
{
}

void cuECorrDouble::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_mskDataPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new(dev_fftTmp1, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT , fftSize);
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

void cuECorrDouble::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);
    if (dev_mskData != NULL)  cuda_arrayDev_delete(dev_mskData);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_fftTmp1);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);

    // cuFFT
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanInv));

    // 3D
    if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
    if (dev_mskDataRotInit != NULL)  cuda_arrayDev_delete(dev_mskDataRotInit);
}

void cuECorrDouble::prepareRef(double *refDataInput)
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

void cuECorrDouble::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL && dev_mskData == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
        cuda_arrayDev_new(dev_mskData, tplSize);
    }
}

void cuECorrDouble::copyTplMsk(double *tplDataInput, double *mskDataInput)
{
    // template
#ifdef __GEM_USE_CUDA_PINNEDMEM__
    CUDA_SAFE_CALL( cudaMemcpyAsync(dev_tplData,
                                    tplDataInput,
                                    tplSize*sizeof(double),
                                    cudaMemcpyHostToDevice) );
#else
    cuda_array_memcpy_h2d(dev_tplData, tplDataInput, tplSize);
#endif

    // mask
    if ((bMskSingle && !bMskComputed) || !bMskSingle) {
#ifdef __GEM_USE_CUDA_PINNEDMEM__
        CUDA_SAFE_CALL( cudaMemcpyAsync(dev_mskData,
                                        mskDataInput,
                                        tplSize*sizeof(double),
                                        cudaMemcpyHostToDevice) );
#else
        cuda_array_memcpy_h2d(dev_mskData, mskDataInput, tplSize);
#endif
    }
}

void cuECorrDouble::normalizeTpl(void)
{
    //cuda_array_math_mul(dev_tplData, dev_mskData, tplSize);
    //cuda_array_math_div(dev_tplData, std::sqrt(cuda_array_reduce_sum2(dev_tplData, tplSize)), tplSize);
    cuda_xcorrNormTemplateECC(dev_tplData, tplSize, dev_mskData);
}

void cuECorrDouble::prepareTplMsk(void)
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
            cuda_array_pad(dev_mskData,       tplNrow, tplNcol,
                           dev_mskDataPadded, refNrow, refNcol,
                                              0,       0);
            CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                         (cufftDoubleReal *)    dev_mskDataPadded,
                                         (cufftDoubleComplex *) dev_mskDataPaddedFFT));
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
            cuda_array_pad(dev_mskData,       tplNrow, tplNcol, tplNsec,
                           dev_mskDataPadded, refNrow, refNcol, refNsec,
                                              0,       0,       0);
            CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                         (cufftDoubleReal *)    dev_mskDataPadded,
                                         (cufftDoubleComplex *) dev_mskDataPaddedFFT));
        }
    }
}

void cuECorrDouble::computeCorr(void)
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
            (dev_mskDataPaddedFFT,
             dev_refDataPaddedFFT2,
             dev_resDataPaddedFFT,
             fftSize,
             (double) refSize);
        CUFFT_SAFE_CALL(cufftExecZ2D(fftPlanInv,
                                     (cufftDoubleComplex *) dev_resDataPaddedFFT,
                                     (cufftDoubleReal *)    dev_fftTmp1));

        // lock the next iteration
        bMskComputed = true;
    }

    // ECC
    //cuda_array_math_sqrt(dev_fftTmp1, refSize);
    //cuda_array_math_div(dev_resDataPadded, dev_fftTmp1, refSize);
    cuda_xcorrCombineResultECC(dev_resDataPadded, dev_fftTmp1, refSize);
}

double cuECorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

double* cuECorrDouble::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuECorrDouble::getResult(double *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuECorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuECorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
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

void cuECorrDouble::copyTplMskRot3D(double *tplDataInput, double *mskDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(mskDataInput != NULL);
    assert(dev_tplDataRotInit == NULL);
    assert(dev_mskDataRotInit == NULL);

    cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
    cuda_arrayDev_new(dev_mskDataRotInit, tplSize);
    cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
    cuda_array_memcpy_h2d(dev_mskDataRotInit, mskDataInput, tplSize);
}

void cuECorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplDataRotInit != NULL && dev_tplData != NULL);

    cuda_transform_rotate(dev_tplDataRotInit, dev_tplData,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
}

void cuECorrDouble::mskCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_mskDataRotInit != NULL && dev_mskData != NULL);

    cuda_transform_rotate(dev_mskDataRotInit, dev_mskData,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
    cuda_array_threshold(dev_mskData, tplSize, 0.5);
}

} // namespace gem
