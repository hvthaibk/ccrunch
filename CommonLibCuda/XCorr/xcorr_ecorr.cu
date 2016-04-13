/***********************************************************************
 *  File:       xcorr_ecorr.cu
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
 * Normalized cross-correlation
 **********************************************************************/

// 1D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow,
                T *dev_refData, size_t refNrow,
                T *dev_resData,
                eXCorrRes shape)
{
    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_mskData = NULL;
    size_t    resNrow = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            break;
        default:
            ERROR("cuda_ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // template normalization
    cuda_arrayDev_new(dev_mskData, tplSize);
    cuda_array_setval(dev_mskData, (T) 1, tplSize);
    cuda_xcorrNormTemplateECC(dev_tplDataTmp, tplSize);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow,
               dev_refDataTmp, refNrow,
               dev_resData,
               shape);

    // denominator
    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow,
               dev_refDataTmp, refNrow,
               dev_fftTmp1,
               shape);

    // ECC
    cuda_xcorrCombineResultECC(dev_resData,
                               dev_fftTmp1,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_mskData);
    cuda_arrayDev_delete(dev_fftTmp1);
}

// instantiation
template
void cuda_ecorr<float >(float  *dev_tplData, size_t tplNrow,
                        float  *dev_refData, size_t refNrow,
                        float  *dev_resData,
                        eXCorrRes shape);
template
void cuda_ecorr<double>(double *dev_tplData, size_t tplNrow,
                        double *dev_refData, size_t refNrow,
                        double *dev_resData,
                        eXCorrRes shape);

// 2D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow, size_t tplNcol,
                T *dev_refData, size_t refNrow, size_t refNcol,
                T *dev_resData,
                eXCorrRes shape)
{
    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_mskData = NULL;
    size_t    resNrow = 0, resNcol = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            resNcol = refNcol - tplNcol + 1;
            break;
        default:
            ERROR("cuda_ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // template normalization
    cuda_arrayDev_new(dev_mskData, tplSize);
    cuda_array_setval(dev_mskData, (T) 1, tplSize);
    cuda_xcorrNormTemplateECC(dev_tplDataTmp, tplSize);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_resData,
               shape);

    // denominator
    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol,
               dev_refDataTmp, refNrow, refNcol,
               dev_fftTmp1,
               shape);

    // ECC
    cuda_xcorrCombineResultECC(dev_resData,
                               dev_fftTmp1,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_mskData);
    cuda_arrayDev_delete(dev_fftTmp1);
}

// instantiation
template
void cuda_ecorr<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol,
                        float  *dev_refData, size_t refNrow, size_t refNcol,
                        float  *dev_resData,
                        eXCorrRes shape);
template
void cuda_ecorr<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol,
                        double *dev_refData, size_t refNrow, size_t refNcol,
                        double *dev_resData,
                        eXCorrRes shape);

// 3D
template <typename T>
void cuda_ecorr(T *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                T *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                T *dev_resData,
                eXCorrRes shape)
{
    T         *dev_tplDataTmp = NULL,
              *dev_refDataTmp = NULL,
              *dev_fftTmp1 = NULL,
              *dev_mskData = NULL;
    size_t    resNrow = 0, resNcol = 0, resNsec = 0, tplSize, refSize, resSize;

    if (shape != XCORR_RES_FULL && shape != XCORR_RES_VALID) {
        shape = XCORR_RES_FULL;
    }
    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            resNcol = refNcol + tplNcol - 1;
            resNsec = refNsec + tplNsec - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            resNcol = refNcol - tplNcol + 1;
            resNsec = refNsec - tplNsec + 1;
            break;
        default:
            ERROR("cuda_ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays as copies of input arrays
    cuda_arrayDev_new(dev_tplDataTmp, tplSize);
    cuda_arrayDev_new(dev_refDataTmp, refSize);
    cuda_array_memcpy_d2d(dev_tplDataTmp, dev_tplData, tplSize);
    cuda_array_memcpy_d2d(dev_refDataTmp, dev_refData, refSize);

    // template normalization
    cuda_arrayDev_new(dev_mskData, tplSize);
    cuda_array_setval(dev_mskData, (T) 1, tplSize);
    cuda_xcorrNormTemplateECC(dev_tplDataTmp, tplSize);

    // allocate temporary arrays
    cuda_arrayDev_new(dev_fftTmp1, resSize);

    // numerator
    cuda_xcorr(dev_tplDataTmp, tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_resData,
               shape);

    // denominator
    cuda_array_math_sqr(dev_refDataTmp, refSize);
    cuda_xcorr(dev_mskData,    tplNrow, tplNcol, tplNsec,
               dev_refDataTmp, refNrow, refNcol, refNsec,
               dev_fftTmp1,
               shape);

    // ECC
    cuda_xcorrCombineResultECC(dev_resData,
                               dev_fftTmp1,
                               resSize);

    // deallocate
    cuda_arrayDev_delete(dev_tplDataTmp);
    cuda_arrayDev_delete(dev_refDataTmp);
    cuda_arrayDev_delete(dev_mskData);
    cuda_arrayDev_delete(dev_fftTmp1);
}

// instantiation
template
void cuda_ecorr<float >(float  *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                        float  *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                        float  *dev_resData,
                        eXCorrRes shape);
template
void cuda_ecorr<double>(double *dev_tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                        double *dev_refData, size_t refNrow, size_t refNcol, size_t refNsec,
                        double *dev_resData,
                        eXCorrRes shape);

} // namespace gem
