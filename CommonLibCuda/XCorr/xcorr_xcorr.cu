/***********************************************************************
 *  File:       xcorr_xcorr.cu
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
 * Cross-correlation
 **********************************************************************/

// 1D float
void cuda_xcorr(const float* const tplData, size_t tplNrow,
                const float* const refData, size_t refNrow,
                      float* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t            fftSize, resNrow, resSize;
    float             *tplDataPadded = NULL,
                      *refDataPadded = NULL,
                      *resDataPadded = NULL;
    cuFloatComplex    *tplDataPaddedFFT = NULL,
                      *refDataPaddedFFT = NULL,
                      *resDataPaddedFFT = NULL;
    cufftHandle       fftPlanFwd, fftPlanInv;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow / 2 + 1;
    resSize = resNrow;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow,
                           tplDataPadded, resNrow,
                                          0);
            cuda_array_pad(refData,       refNrow,
                           refDataPadded, resNrow,
                                          tplNrow-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow,
                           tplDataPadded, resNrow,
                                          0);
            cuda_array_pad(refData,       refNrow,
                           refDataPadded, resNrow,
                                          0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (FLOAT)
    cufftPlan1d(&fftPlanFwd, (int) resNrow, CUFFT_R2C, 1);
    cufftPlan1d(&fftPlanInv, (int) resNrow, CUFFT_C2R, 1);

    // FFT forward (FLOAT)
    cufftExecR2C(fftPlanFwd, (cufftReal *)    tplDataPadded,
                             (cufftComplex *) tplDataPaddedFFT);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    refDataPadded,
                             (cufftComplex *) refDataPaddedFFT);

    // power spectrum (FLOAT)
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (float) resSize);

    // FFT inverse (FLOAT)
    cufftExecC2R(fftPlanInv, (cufftComplex *) resDataPaddedFFT,
                             (cufftReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow,
                            resData,       resNrow,
                                           0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,
                            resData,       resNrow-tplNrow+1,
                                           0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

// 1D double
void cuda_xcorr(const double* const tplData, size_t tplNrow,
                const double* const refData, size_t refNrow,
                      double* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    size_t             fftSize, resNrow, resSize;
    double             *tplDataPadded = NULL,
                       *refDataPadded = NULL,
                       *resDataPadded = NULL;
    cuDoubleComplex    *tplDataPaddedFFT = NULL,
                       *refDataPaddedFFT = NULL,
                       *resDataPaddedFFT = NULL;
    cufftHandle        fftPlanFwd, fftPlanInv;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow;
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow / 2 + 1;
    resSize = resNrow;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow,
                           tplDataPadded, resNrow,
                                          0);
            cuda_array_pad(refData,       refNrow,
                           refDataPadded, resNrow,
                                          tplNrow-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow,
                           tplDataPadded, resNrow,
                                          0);
            cuda_array_pad(refData,       refNrow,
                           refDataPadded, resNrow,
                                          0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (DOUBLE)
    cufftPlan1d(&fftPlanFwd, (int) resNrow, CUFFT_D2Z, 1);
    cufftPlan1d(&fftPlanInv, (int) resNrow, CUFFT_Z2D, 1);

    // FFT forward (DOUBLE)
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    tplDataPadded,
                             (cufftDoubleComplex *) tplDataPaddedFFT);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    refDataPadded,
                             (cufftDoubleComplex *) refDataPaddedFFT);

    // power spectrum (DOUBLE)
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (double) resSize);

    // FFT inverse (DOUBLE)
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) resDataPaddedFFT,
                             (cufftDoubleReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow,
                            resData,       resNrow,
                                           0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,
                            resData,       resNrow-tplNrow+1,
                                           0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

// 2D float
void cuda_xcorr(const float* const tplData, size_t tplNrow, size_t tplNcol,
                const float* const refData, size_t refNrow, size_t refNcol,
                      float* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t            fftSize, resNrow, resNcol, resSize;
    float             *tplDataPadded = NULL,
                      *refDataPadded = NULL,
                      *resDataPadded = NULL;
    cuFloatComplex    *tplDataPaddedFFT = NULL,
                      *refDataPaddedFFT = NULL,
                      *resDataPaddedFFT = NULL;
    cufftHandle       fftPlanFwd, fftPlanInv;

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
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * (resNcol / 2 + 1);
    resSize = resNrow * resNcol;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow, tplNcol,
                           tplDataPadded, resNrow, resNcol,
                                          0,       0);
            cuda_array_pad(refData,       refNrow,   refNcol,
                           refDataPadded, resNrow,   resNcol,
                                          tplNrow-1, tplNcol-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow, tplNcol,
                           tplDataPadded, resNrow, resNcol,
                                          0,       0);
            cuda_array_pad(refData,       refNrow, refNcol,
                           refDataPadded, resNrow, resNcol,
                                          0,       0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (FLOAT)
    cufftPlan2d(&fftPlanFwd, (int) resNrow, (int) resNcol, CUFFT_R2C);
    cufftPlan2d(&fftPlanInv, (int) resNrow, (int) resNcol, CUFFT_C2R);

    // FFT forward (FLOAT)
    cufftExecR2C(fftPlanFwd, (cufftReal *)    tplDataPadded,
                             (cufftComplex *) tplDataPaddedFFT);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    refDataPadded,
                             (cufftComplex *) refDataPaddedFFT);

    // power spectrum
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (float) resSize);

    // FFT inverse (FLOAT)
    cufftExecC2R(fftPlanInv, (cufftComplex *) resDataPaddedFFT,
                             (cufftReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow, resNcol,
                            resData,       resNrow, resNcol,
                                           0,       0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,           resNcol,
                            resData,       resNrow-tplNrow+1, resNcol-tplNcol+1,
                                           0,                 0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

// 2D double
void cuda_xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol,
                const double* const refData, size_t refNrow, size_t refNcol,
                      double* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    size_t             fftSize, resNrow, resNcol, resSize;
    double             *tplDataPadded = NULL,
                       *refDataPadded = NULL,
                       *resDataPadded = NULL;
    cuDoubleComplex    *tplDataPaddedFFT = NULL,
                       *refDataPaddedFFT = NULL,
                       *resDataPaddedFFT = NULL;
    cufftHandle        fftPlanFwd, fftPlanInv;

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
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * (resNcol / 2 + 1);
    resSize = resNrow * resNcol;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow, tplNcol,
                           tplDataPadded, resNrow, resNcol,
                                          0,       0);
            cuda_array_pad(refData,       refNrow,   refNcol,
                           refDataPadded, resNrow,   resNcol,
                                          tplNrow-1, tplNcol-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow, tplNcol,
                           tplDataPadded, resNrow, resNcol,
                                          0,       0);
            cuda_array_pad(refData,       refNrow, refNcol,
                           refDataPadded, resNrow, resNcol,
                                          0,       0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (DOUBLE)
    cufftPlan2d(&fftPlanFwd, (int) resNrow, (int) resNcol, CUFFT_D2Z);
    cufftPlan2d(&fftPlanInv, (int) resNrow, (int) resNcol, CUFFT_Z2D);

    // FFT forward (DOUBLE)
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    tplDataPadded,
                             (cufftDoubleComplex *) tplDataPaddedFFT);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    refDataPadded,
                             (cufftDoubleComplex *) refDataPaddedFFT);

    // power spectrum (DOUBLE)
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (double) resSize);

    // FFT inverse (DOUBLE)
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) resDataPaddedFFT,
                             (cufftDoubleReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow, resNcol,
                            resData,       resNrow, resNcol,
                                           0,       0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,           resNcol,
                            resData,       resNrow-tplNrow+1, resNcol-tplNcol+1,
                                           0,                 0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

// 3D float
void cuda_xcorr(const float* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                const float* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                      float* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t            fftSize, resNrow, resNcol, resNsec, resSize;
    float             *tplDataPadded = NULL,
                      *refDataPadded = NULL,
                      *resDataPadded = NULL;
    cuFloatComplex    *tplDataPaddedFFT = NULL,
                      *refDataPaddedFFT = NULL,
                      *resDataPaddedFFT = NULL;
    cufftHandle       fftPlanFwd, fftPlanInv;

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
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * resNcol * (resNsec / 2 + 1);
    resSize = resNrow * resNcol * resNsec;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                           tplDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            cuda_array_pad(refData,       refNrow,   refNcol,   refNsec,
                           refDataPadded, resNrow,   resNcol,   resNsec,
                                          tplNrow-1, tplNcol-1, tplNsec-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                           tplDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            cuda_array_pad(refData,       refNrow, refNcol, refNsec,
                           refDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (FLOAT)
    cufftPlan3d(&fftPlanFwd, (int) resNrow, (int) resNcol, (int) resNsec, CUFFT_R2C);
    cufftPlan3d(&fftPlanInv, (int) resNrow, (int) resNcol, (int) resNsec, CUFFT_C2R);

    // FFT forward (FLOAT)
    cufftExecR2C(fftPlanFwd, (cufftReal *)    tplDataPadded,
                             (cufftComplex *) tplDataPaddedFFT);
    cufftExecR2C(fftPlanFwd, (cufftReal *)    refDataPadded,
                             (cufftComplex *) refDataPaddedFFT);

    // power spectrum
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (float) resSize);

    // FFT inverse (FLOAT)
    cufftExecC2R(fftPlanInv, (cufftComplex *) resDataPaddedFFT,
                             (cufftReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow, resNcol, resNsec,
                            resData,       resNrow, resNcol, resNsec,
                                           0,       0,       0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,           resNcol,           resNsec,
                            resData,       resNrow-tplNrow+1, resNcol-tplNcol+1, resNsec-tplNsec+1,
                                           0,                 0,                 0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

// 3D double
void cuda_xcorr(const double* const tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                const double* const refData, size_t refNrow, size_t refNcol, size_t refNsec,
                      double* const resData,
                eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    size_t             fftSize, resNrow, resNcol, resNsec, resSize;
    double             *tplDataPadded = NULL,
                       *refDataPadded = NULL,
                       *resDataPadded = NULL;
    cuDoubleComplex    *tplDataPaddedFFT = NULL,
                       *refDataPaddedFFT = NULL,
                       *resDataPaddedFFT = NULL;
    cufftHandle        fftPlanFwd, fftPlanInv;

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
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    fftSize = resNrow * resNcol * (resNsec / 2 + 1);
    resSize = resNrow * resNcol * resNsec;

    // allocate GPU memory
    cuda_arrayDev_new(tplDataPadded,    resSize);
    cuda_arrayDev_new(refDataPadded,    resSize);
    cuda_arrayDev_new(resDataPadded,    resSize);
    cuda_arrayDev_new(tplDataPaddedFFT, fftSize);
    cuda_arrayDev_new(refDataPaddedFFT, fftSize);
    cuda_arrayDev_new(resDataPaddedFFT, fftSize);

    // padding arrays in GPU
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                           tplDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            cuda_array_pad(refData,       refNrow,   refNcol,   refNsec,
                           refDataPadded, resNrow,   resNcol,   resNsec,
                                          tplNrow-1, tplNcol-1, tplNsec-1);
            break;
        case XCORR_RES_VALID:
            cuda_array_pad(tplData,       tplNrow, tplNcol, tplNsec,
                           tplDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            cuda_array_pad(refData,       refNrow, refNcol, refNsec,
                           refDataPadded, resNrow, resNcol, resNsec,
                                          0,       0,       0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // FFT plans (DOUBLE)
    cufftPlan3d(&fftPlanFwd, (int) resNrow, (int) resNcol, (int) resNsec, CUFFT_D2Z);
    cufftPlan3d(&fftPlanInv, (int) resNrow, (int) resNcol, (int) resNsec, CUFFT_Z2D);

    // FFT forward (DOUBLE)
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    tplDataPadded,
                             (cufftDoubleComplex *) tplDataPaddedFFT);
    cufftExecD2Z(fftPlanFwd, (cufftDoubleReal *)    refDataPadded,
                             (cufftDoubleComplex *) refDataPaddedFFT);

    // power spectrum (DOUBLE)
    cuda_xcorrModulateAndNormalize(tplDataPaddedFFT,
                                   refDataPaddedFFT,
                                   resDataPaddedFFT,
                                   fftSize,
                                   (double) resSize);

    // FFT inverse (DOUBLE)
    cufftExecZ2D(fftPlanInv, (cufftDoubleComplex *) resDataPaddedFFT,
                             (cufftDoubleReal *)    resDataPadded);

    // extract result from FFTW output array
    switch (shape) {
        case XCORR_RES_FULL:
            cuda_array_crop(resDataPadded, resNrow, resNcol, resNsec,
                            resData,       resNrow, resNcol, resNsec,
                                           0,       0,       0);
            break;
        case XCORR_RES_VALID:
            cuda_array_crop(resDataPadded, resNrow,           resNcol,           resNsec,
                            resData,       resNrow-tplNrow+1, resNcol-tplNcol+1, resNsec-tplNsec+1,
                                           0,                 0,                 0);
            break;
        default:
            ERROR("cuda_xcorr", "unsupported xcorr mode");
    }

    // free GPU memory
    cuda_arrayDev_delete(tplDataPadded);
    cuda_arrayDev_delete(refDataPadded);
    cuda_arrayDev_delete(resDataPadded);
    cuda_arrayDev_delete(tplDataPaddedFFT);
    cuda_arrayDev_delete(refDataPaddedFFT);
    cuda_arrayDev_delete(resDataPaddedFFT);

    cufftDestroy(fftPlanFwd);
    cufftDestroy(fftPlanInv);
}

} // namespace gem
