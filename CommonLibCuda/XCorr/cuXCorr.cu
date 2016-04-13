/***********************************************************************
 *  File:       cuXCorr.cu
 *
 *  Purpose:    Implementation of XCC computing classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuXCorr.cuh"

namespace gem {

cuXCorrSingle::cuXCorrSingle(size_t refNrowInput,
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

    dev_tplData = NULL;

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    bRotTexture = false;
    dev_tplDataRotInit = NULL;
    dev_tplDataRotInit_tex3D = NULL;
}

cuXCorrSingle::~cuXCorrSingle()
{
}

void cuXCorrSingle::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
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

void cuXCorrSingle::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);

    // cuFFT
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanInv));

    // 3D
    if (bRotTexture) {
        if (dev_tplDataRotInit_tex3D != NULL)  cuda_arrayArr_delete(dev_tplDataRotInit_tex3D);
    }
    else {
        if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
    }
}

void cuXCorrSingle::prepareRef(float *refDataInput)
{
    cuda_array_memcpy_h2d(dev_refDataPadded, refDataInput, refSize);

    CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                 (cufftReal *)    dev_refDataPadded,
                                 (cufftComplex *) dev_refDataPaddedFFT));
}

void cuXCorrSingle::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
    }
}

void cuXCorrSingle::copyTpl(float *tplDataInput)
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
}

void cuXCorrSingle::prepareTpl(void)
{
    if (refNsec == 1) {       // 2D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol,
                       dev_tplDataPadded, refNrow, refNcol,
                                          0,       0);
        CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                     (cufftReal *)    dev_tplDataPadded,
                                     (cufftComplex *) dev_tplDataPaddedFFT));
    }
    else {                    // 3D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol, tplNsec,
                       dev_tplDataPadded, refNrow, refNcol, refNsec,
                                          0,       0,       0);
        CUFFT_SAFE_CALL(cufftExecR2C(fftPlanFwd,
                                     (cufftReal *)    dev_tplDataPadded,
                                     (cufftComplex *) dev_tplDataPaddedFFT));
    }
}

void cuXCorrSingle::computeCorr(void)
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
}

float cuXCorrSingle::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

float* cuXCorrSingle::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuXCorrSingle::getResult(float *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuXCorrSingle::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuXCorrSingle::mergeResultGlobal(float  *resDataAbsMaxPaddedGlobal,
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

void cuXCorrSingle::copyTplRot3D(float *tplDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(dev_tplDataRotInit == NULL && dev_tplDataRotInit_tex3D == NULL);

    if (bRotTexture) {
        cuda_arrayArr_new(dev_tplDataRotInit_tex3D, tplNrow, tplNcol, tplNsec);
        cuda_array_memcpy_h2d(dev_tplDataRotInit_tex3D, tplDataInput, tplNrow, tplNcol, tplNsec);
    }
    else {
        cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
        cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
    }
}

void cuXCorrSingle::tplCreateByRot3D(float alpha, float beta, float gamma, eInter inter)
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

cuXCorrDouble::cuXCorrDouble(size_t refNrowInput,
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

    dev_tplData = NULL;

    dev_tplDataPadded = NULL;
    dev_refDataPadded = NULL;
    dev_resDataPadded = NULL;

    dev_tplDataPaddedFFT  = NULL;
    dev_refDataPaddedFFT  = NULL;
    dev_resDataPaddedFFT  = NULL;

    dev_resDataAbsMaxPadded = NULL;
    dev_resDataMaxIndPadded = NULL;

    resDataAbsMaxPadded = NULL;
    resDataMaxIndPadded = NULL;

    // 3D
    dev_tplDataRotInit = NULL;
}

cuXCorrDouble::~cuXCorrDouble()
{
}

void cuXCorrDouble::memAlloc(void)
{
    // memory allocation
    cuda_arrayDev_new(dev_tplDataPadded, refSize);
    cuda_arrayDev_new(dev_refDataPadded, refSize);
    cuda_arrayDev_new(dev_resDataPadded, refSize);

    cuda_arrayDev_new_zero(dev_resDataAbsMaxPadded, refSize);
    cuda_arrayDev_new     (dev_resDataMaxIndPadded, refSize);

    array_new(resDataAbsMaxPadded, refSize);
    array_new(resDataMaxIndPadded, refSize);

    // complex data
    cuda_arrayDev_new(dev_tplDataPaddedFFT , fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT , fftSize);
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

void cuXCorrDouble::memFree(void)
{
    // memory deallocation
    if (dev_tplData != NULL)  cuda_arrayDev_delete(dev_tplData);

    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);

    cuda_arrayDev_delete(dev_resDataAbsMaxPadded);
    cuda_arrayDev_delete(dev_resDataMaxIndPadded);

    array_delete(resDataAbsMaxPadded);
    array_delete(resDataMaxIndPadded);

    // complex data
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);

    // cuFFT
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(fftPlanInv));

    // 3D
    if (dev_tplDataRotInit != NULL)  cuda_arrayDev_delete(dev_tplDataRotInit);
}

void cuXCorrDouble::prepareRef(double *refDataInput)
{
    cuda_array_memcpy_h2d(dev_refDataPadded, refDataInput, refSize);

    CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                 (cufftDoubleReal *)    dev_refDataPadded,
                                 (cufftDoubleComplex *) dev_refDataPaddedFFT));
}

void cuXCorrDouble::setSizeTpl(size_t tplNrowInput,
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

    if (dev_tplData == NULL) {
        cuda_arrayDev_new(dev_tplData, tplSize);
    }
}

void cuXCorrDouble::copyTpl(double *tplDataInput)
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
}

void cuXCorrDouble::prepareTpl(void)
{
    if (refNsec == 1) {       // 2D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol,
                       dev_tplDataPadded, refNrow, refNcol,
                                          0,       0);
        CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                     (cufftDoubleReal *)    dev_tplDataPadded,
                                     (cufftDoubleComplex *) dev_tplDataPaddedFFT));
    }
    else {                    // 3D
        // template
        cuda_array_pad(dev_tplData,       tplNrow, tplNcol, tplNsec,
                       dev_tplDataPadded, refNrow, refNcol, refNsec,
                                          0,       0,       0);
        CUFFT_SAFE_CALL(cufftExecD2Z(fftPlanFwd,
                                     (cufftDoubleReal *)    dev_tplDataPadded,
                                     (cufftDoubleComplex *) dev_tplDataPaddedFFT));
    }
}

void cuXCorrDouble::computeCorr(void)
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
}

double cuXCorrDouble::getResultMax(void)
{
    if (refNsec == 1) {       // 2D
        return cuda_array_reduce_min(dev_resDataPadded, refSize);
    }
    else {                    // 3D
        return cuda_array_reduce_max(dev_resDataPadded, refSize);
    }
}

double* cuXCorrDouble::getResultRef(void)
{
    return dev_resDataPadded;
}

void cuXCorrDouble::getResult(double *resDataPaddedOut)
{
    cuda_array_memcpy_d2h(resDataPaddedOut, dev_resDataPadded, refSize);
}

void cuXCorrDouble::mergeResult(size_t indx, eXCorrMerge bAbs)
{
    cuda_xcorrMergeResult(dev_resDataPadded,
                          dev_resDataAbsMaxPadded,
                          dev_resDataMaxIndPadded,
                          refSize,
                          indx,
                          bAbs);
}

void cuXCorrDouble::mergeResultGlobal(double *resDataAbsMaxPaddedGlobal,
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

void cuXCorrDouble::copyTplRot3D(double *tplDataInput)
{
    assert(tplSize > 0);
    assert(tplDataInput != NULL);
    assert(dev_tplDataRotInit == NULL);

    cuda_arrayDev_new(dev_tplDataRotInit, tplSize);
    cuda_array_memcpy_h2d(dev_tplDataRotInit, tplDataInput, tplSize);
}

void cuXCorrDouble::tplCreateByRot3D(double alpha, double beta, double gamma, eInter inter)
{
    assert(tplSize > 0);
    assert(dev_tplDataRotInit != NULL && dev_tplData != NULL);

    cuda_transform_rotate(dev_tplDataRotInit, dev_tplData,
                          tplNrow, tplNcol, tplNsec,
                          alpha, beta, gamma,
                          inter);
}

} // namespace gem
