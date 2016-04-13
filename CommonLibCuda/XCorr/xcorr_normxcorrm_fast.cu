/***********************************************************************
 *  File:       xcorr_normxcorrm_fast.cu
 *
 *  Purpose:    Implementation of xcorr-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "xcorr.cuh"

namespace gem {

/***********************************************************************
 * Combination for fast computation
 **********************************************************************/

// 1D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow,
                          float *dev_refData, size_t refNrow,
                          float *dev_resData,
                          float *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t              resNrow, tplSize, refSize, fftSize;
    float               *dev_tplDataTmp = NULL;
    float               *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuFloatComplex      *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    float               *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    tplSize = tplNrow;
    refSize = refNrow;
    fftSize = refNrow / 2 + 1;

    // FFT plans (FLOAT)
    cufftPlan1d(&fftPlanFwd, (int) refNrow, CUFFT_R2C, 1);
    cufftPlan1d(&fftPlanInv, (int) refNrow, CUFFT_C2R, 1);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow,
                   dev_tplDataPadded, refNrow,
                   0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_tplDataPadded,
                             (cufftComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow,
                   dev_mskDataPadded, refNrow,
                   0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_mskDataPadded,
                             (cufftComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow,
                    dev_resData,       resNrow,
                                       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// 1D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow,
                          double *dev_refData, size_t refNrow,
                          double *dev_resData,
                          double *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t              resNrow, tplSize, refSize, fftSize;
    double              *dev_tplDataTmp = NULL;
    double              *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuDoubleComplex     *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    double              *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    tplSize = tplNrow;
    refSize = refNrow;
    fftSize = refNrow / 2 + 1;

    // FFT plans (DOUBLE)
    cufftPlan1d(&fftPlanFwd, (int) refNrow, CUFFT_D2Z, 1);
    cufftPlan1d(&fftPlanInv, (int) refNrow, CUFFT_Z2D, 1);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow,
                   dev_tplDataPadded, refNrow,
                                      0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_tplDataPadded,
                             (cufftDoubleComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow,
                   dev_mskDataPadded, refNrow,
                                      0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_mskDataPadded,
                             (cufftDoubleComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow,
                    dev_resData,       resNrow,
                                       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// 2D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow, size_t tplNcol,
                          float *dev_refData, size_t refNrow, size_t refNcol,
                          float *dev_resData,
                          float *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t              resNrow, resNcol, tplSize, refSize, fftSize;
    float               *dev_tplDataTmp = NULL;
    float               *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuFloatComplex      *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    float               *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    fftSize = refNrow * (refNcol / 2 + 1);

    // FFT plans (FLOAT)
    cufftPlan2d(&fftPlanFwd, (int) refNrow, (int) refNcol, CUFFT_R2C);
    cufftPlan2d(&fftPlanInv, (int) refNrow, (int) refNcol, CUFFT_C2R);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow, tplNcol,
                   dev_tplDataPadded, refNrow, refNcol,
                                      0,       0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_tplDataPadded,
                             (cufftComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow, tplNcol,
                   dev_mskDataPadded, refNrow, refNcol,
                                      0,       0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_mskDataPadded,
                             (cufftComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow, refNcol,
                    dev_resData,       resNrow, resNcol,
                                       0,       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// 2D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow, size_t tplNcol,
                          double *dev_refData, size_t refNrow, size_t refNcol,
                          double *dev_resData,
                          double *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t              resNrow, resNcol, tplSize, refSize, fftSize;
    double              *dev_tplDataTmp = NULL;
    double              *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuDoubleComplex     *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    double              *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    fftSize = refNrow * (refNcol / 2 + 1);

    // FFT plans (DOUBLE)
    cufftPlan2d(&fftPlanFwd, (int) refNrow, (int) refNcol, CUFFT_D2Z);
    cufftPlan2d(&fftPlanInv, (int) refNrow, (int) refNcol, CUFFT_Z2D);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow, tplNcol,
                   dev_tplDataPadded, refNrow, refNcol,
                                      0,       0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_tplDataPadded,
                             (cufftDoubleComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow, tplNcol,
                   dev_mskDataPadded, refNrow, refNcol,
                                      0,       0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_mskDataPadded,
                             (cufftDoubleComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow, refNcol,
                    dev_resData,       resNrow, resNcol,
                                       0,       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// 3D float
void cuda_normxcorrm_fast(float *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                          float *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                          float *dev_resData,
                          float *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t              resNrow, resNcol, resNsec, tplSize, refSize, fftSize;
    float               *dev_tplDataTmp = NULL;
    float               *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuFloatComplex      *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    float               *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    resNsec = refNsec - tplNsec + 1;
    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    // FFT plans (FLOAT)
    cufftPlan3d(&fftPlanFwd, (int) refNrow, (int) refNcol, (int) refNsec, CUFFT_R2C);
    cufftPlan3d(&fftPlanInv, (int) refNrow, (int) refNcol, (int) refNsec, CUFFT_C2R);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_refDataPadded,
                             (cufftComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow, tplNcol, tplNsec,
                   dev_tplDataPadded, refNrow, refNcol, refNsec,
                                      0,       0,       0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_tplDataPadded,
                             (cufftComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow, tplNcol, tplNsec,
                   dev_mskDataPadded, refNrow, refNcol, refNsec,
                                      0,       0,       0);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    dev_mskDataPadded,
                             (cufftComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (float) refSize);
    cufftExecC2R(fftPlanInv, (cufftComplex *) dev_resDataPaddedFFT,
                             (cufftReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow, refNcol, refNsec,
                    dev_resData,       resNrow, resNcol, resNsec,
                                       0,       0,       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

// 3D double
void cuda_normxcorrm_fast(double *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                          double *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                          double *dev_resData,
                          double *dev_mskData)
{
    assert(dev_tplData != NULL);
    assert(dev_refData != NULL);
    assert(dev_resData != NULL);
    assert(dev_mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t              resNrow, resNcol, resNsec, tplSize, refSize, fftSize;
    double              *dev_tplDataTmp = NULL;
    double              *dev_tplDataPadded = NULL,
                        *dev_refDataPadded = NULL,
                        *dev_mskDataPadded = NULL,
                        *dev_resDataPadded = NULL;
    cuDoubleComplex     *dev_tplDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT = NULL,
                        *dev_refDataPaddedFFT2 = NULL,
                        *dev_mskDataPaddedFFT = NULL,
                        *dev_resDataPaddedFFT = NULL;
    cufftHandle         fftPlanFwd, fftPlanInv;
    double              *dev_fftTmp1 = NULL, *dev_fftTmp2 = NULL;
    size_t              mskArea;

    resNrow = refNrow - tplNrow + 1;
    resNcol = refNcol - tplNcol + 1;
    resNsec = refNsec - tplNsec + 1;
    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    fftSize = refNrow * refNcol * (refNsec / 2 + 1);

    // FFT plans (DOUBLE)
    cufftPlan3d(&fftPlanFwd, (int) refNrow, (int) refNcol, (int) refNsec, CUFFT_D2Z);
    cufftPlan3d(&fftPlanInv, (int) refNrow, (int) refNcol, (int) refNsec, CUFFT_Z2D);

    cuda_arrayDev_new(dev_tplDataPadded,     refSize);
    cuda_arrayDev_new(dev_refDataPadded,     refSize);
    cuda_arrayDev_new(dev_mskDataPadded,     refSize);
    cuda_arrayDev_new(dev_resDataPadded,     refSize);
    cuda_arrayDev_new(dev_tplDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_refDataPaddedFFT2, fftSize);
    cuda_arrayDev_new(dev_mskDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_resDataPaddedFFT,  fftSize);
    cuda_arrayDev_new(dev_fftTmp1,           refSize);
    cuda_arrayDev_new(dev_fftTmp2,           refSize);

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp,    dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataPadded, dev_refData, refSize);

    // mask size and template normalization
    mskArea = cuda_xcorrNormTemplateNCC(dev_tplDataTmp, tplSize, dev_mskData);

    // reference
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT);
    cuda_array_math_sqr(dev_refDataPadded, refSize);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_refDataPadded,
                             (cufftDoubleComplex *) dev_refDataPaddedFFT2);

    // template
    cuda_array_pad(dev_tplDataTmp,    tplNrow, tplNcol, tplNsec,
                   dev_tplDataPadded, refNrow, refNcol, refNsec,
                                      0,       0,       0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_tplDataPadded,
                             (cufftDoubleComplex *) dev_tplDataPaddedFFT);

    // mask
    cuda_array_pad(dev_mskData,       tplNrow, tplNcol, tplNsec,
                   dev_mskDataPadded, refNrow, refNcol, refNsec,
                                      0,       0,       0);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    dev_mskDataPadded,
                             (cufftDoubleComplex *) dev_mskDataPaddedFFT);

    // xcorr 1
    cuda_xcorrModulateAndNormalize(dev_tplDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_resDataPadded);

    // xcorr 2
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp1);

    // xcorr 3
    cuda_xcorrModulateAndNormalize(dev_mskDataPaddedFFT,
                                   dev_refDataPaddedFFT2,
                                   dev_resDataPaddedFFT,
                                   fftSize,
                                   (double) refSize);
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) dev_resDataPaddedFFT,
                             (cufftDoubleReal *)    dev_fftTmp2);

    // NCC
    cuda_xcorrCombineResultNCC(dev_resDataPadded,
                               dev_fftTmp1,
                               dev_fftTmp2,
                               mskArea,
                               refSize);

    // crop
    cuda_array_crop(dev_resDataPadded, refNrow, refNcol, refNsec,
                    dev_resData,       resNrow, resNcol, resNsec,
                                       0,       0,       0);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);

    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_tplDataPadded);
    cuda_arrayDev_delete(dev_refDataPadded);
    cuda_arrayDev_delete(dev_mskDataPadded);
    cuda_arrayDev_delete(dev_resDataPadded);
    cuda_arrayDev_delete(dev_tplDataPaddedFFT);
    cuda_arrayDev_delete(dev_mskDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT);
    cuda_arrayDev_delete(dev_refDataPaddedFFT2);
    cuda_arrayDev_delete(dev_resDataPaddedFFT);
    cuda_arrayDev_delete(dev_fftTmp1);
    cuda_arrayDev_delete(dev_fftTmp2);
}

} // namespace gem
