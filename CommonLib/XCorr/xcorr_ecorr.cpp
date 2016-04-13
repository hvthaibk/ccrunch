/***********************************************************************
 *  File:       xcorr_ecorr.cpp
 *
 *  Purpose:    Implementation of xcorr-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "xcorr.hpp"

namespace gem {

/***************************************************
 * Normalized cross-correlation
 **************************************************/

// 1D
template <typename T>
void ecorr(T *tplData, size_t tplNrow,
           T *refData, size_t refNrow,
           T *resData,
           eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL,
              *mskData = NULL;
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
            ERROR("ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow;
    refSize = refNrow;
    resSize = resNrow;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    // template normalization
    array_new(mskData, tplSize);
    array_setval(mskData, (T) 1, tplSize);
    xcorrNormTemplateECC(tplDataTmp, tplSize);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow,
          refDataTmp, refNrow,
          resData,
          shape);

    // denominator
    array_math_sqr(refDataTmp, refSize);
    xcorr(mskData,    tplNrow,
          refDataTmp, refNrow,
          fftTmp1,
          shape);

    // ECC
    xcorrCombineResultECC(resData, fftTmp1, resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(mskData);
    array_delete(fftTmp1);
}

// instantiation
template
void ecorr<float >(float  *tplData, size_t tplNrow,
                   float  *refData, size_t refNrow,
                   float  *resData,
                   eXCorrRes shape);
template
void ecorr<double>(double *tplData, size_t tplNrow,
                   double *refData, size_t refNrow,
                   double *resData,
                   eXCorrRes shape);

// 2D
template <typename T>
void ecorr(T *tplData, size_t tplNrow, size_t tplNcol,
           T *refData, size_t refNrow, size_t refNcol,
           T *resData,
           eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL,
              *mskData = NULL;
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
            ERROR("ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol;
    refSize = refNrow * refNcol;
    resSize = resNrow * resNcol;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    // template normalization
    array_new(mskData, tplSize);
    array_setval(mskData, (T) 1, tplSize);
    xcorrNormTemplateECC(tplDataTmp, tplSize);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow, tplNcol,
          refDataTmp, refNrow, refNcol,
          resData,
          shape);

    // denominator
    array_math_sqr(refDataTmp, refSize);
    xcorr(mskData,    tplNrow, tplNcol,
          refDataTmp, refNrow, refNcol,
          fftTmp1,
          shape);

    // ECC
    xcorrCombineResultECC(resData, fftTmp1, resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(mskData);
    array_delete(fftTmp1);
}

// instantiation
template
void ecorr<float >(float  *tplData, size_t tplNrow, size_t tplNcol,
                   float  *refData, size_t refNrow, size_t refNcol,
                   float  *resData,
                   eXCorrRes shape);
template
void ecorr<double>(double *tplData, size_t tplNrow, size_t tplNcol,
                   double *refData, size_t refNrow, size_t refNcol,
                   double *resData,
                   eXCorrRes shape);

// 3D
template <typename T>
void ecorr(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
           T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
           T *resData,
           eXCorrRes shape)
{
    assert(tplData != NULL && refData != NULL && resData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL,
              *mskData = NULL;
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
            ERROR("ecorr", "unsupported xcorr mode");
    }

    tplSize = tplNrow * tplNcol * tplNsec;
    refSize = refNrow * refNcol * refNsec;
    resSize = resNrow * resNcol * resNsec;

    // allocate arrays as copies of input arrays
    array_new(tplDataTmp, tplSize);
    array_new(refDataTmp, refSize);
    memcpy(tplDataTmp, tplData, tplSize*sizeof(T));
    memcpy(refDataTmp, refData, refSize*sizeof(T));

    // template normalization
    array_new(mskData, tplSize);
    array_setval(mskData, (T) 1, tplSize);
    xcorrNormTemplateECC(tplDataTmp, tplSize);

    // allocate temporary arrays
    array_new(fftTmp1, resSize);

    // numerator
    xcorr(tplDataTmp, tplNrow, tplNcol, tplNsec,
          refDataTmp, refNrow, refNcol, refNsec,
          resData,
          shape);

    // denominator
    array_math_sqr(refDataTmp, refSize);
    xcorr(mskData,    tplNrow, tplNcol, tplNsec,
          refDataTmp, refNrow, refNcol, refNsec,
          fftTmp1,
          shape);

    // ECC
    xcorrCombineResultECC(resData, fftTmp1, resSize);

    // deallocate
    array_delete(tplDataTmp);
    array_delete(refDataTmp);
    array_delete(mskData);
    array_delete(fftTmp1);
}

// instantiation
template
void ecorr<float >(float  *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                   float  *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                   float  *resData,
                   eXCorrRes shape);
template
void ecorr<double>(double *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                   double *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                   double *resData,
                   eXCorrRes shape);

} // namespace gem
