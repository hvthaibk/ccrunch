/***********************************************************************
 *  File:       xcorr_ecorrm.cpp
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
 * Normalized cross-correlation with mask
 **************************************************/

// 1D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow,
            T *refData, size_t refNrow,
            T *resData,
            T *mskData,
            eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0);
    assert(refNrow > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL;
    size_t    resNrow, tplSize, refSize, resSize;

    switch (shape) {
        case XCORR_RES_FULL:
            resNrow = refNrow + tplNrow - 1;
            break;
        case XCORR_RES_VALID:
            resNrow = refNrow - tplNrow + 1;
            break;
        default:
            ERROR("ecorrm", "unsupported xcorr mode");
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
    xcorrNormTemplateECC(tplDataTmp, tplSize, mskData);

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
    array_delete(fftTmp1);
}

// instantiation
template
void ecorrm<float >(float  *tplData, size_t tplNrow,
                    float  *refData, size_t refNrow,
                    float  *resData,
                    float  *mskData,
                    eXCorrRes shape);
template
void ecorrm<double>(double *tplData, size_t tplNrow,
                    double *refData, size_t refNrow,
                    double *resData,
                    double *mskData,
                    eXCorrRes shape);

// 2D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow, size_t tplNcol,
            T *refData, size_t refNrow, size_t refNcol,
            T *resData,
            T *mskData,
            eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0);
    assert(refNrow > 0 && refNcol > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL;
    size_t    resNrow, resNcol, tplSize, refSize, resSize;

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
            ERROR("ecorrm", "unsupported xcorr mode");
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
    xcorrNormTemplateECC(tplDataTmp, tplSize, mskData);

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
    array_delete(fftTmp1);
}

// instantiation
template
void ecorrm<float >(float  *tplData, size_t tplNrow, size_t tplNcol,
                    float  *refData, size_t refNrow, size_t refNcol,
                    float  *resData,
                    float  *mskData,
                    eXCorrRes shape);
template
void ecorrm<double>(double *tplData, size_t tplNrow, size_t tplNcol,
                    double *refData, size_t refNrow, size_t refNcol,
                    double *resData,
                    double *mskData,
                    eXCorrRes shape);

// 3D
template <typename T>
void ecorrm(T *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
            T *refData, size_t refNrow, size_t refNcol, size_t refNsec,
            T *resData,
            T *mskData,
            eXCorrRes shape)
{
    assert(tplData != NULL);
    assert(refData != NULL);
    assert(resData != NULL);
    assert(mskData != NULL);
    assert(tplNrow > 0 && tplNcol > 0 && tplNsec > 0);
    assert(refNrow > 0 && refNcol > 0 && refNsec > 0);

    T         *tplDataTmp = NULL,
              *refDataTmp = NULL,
              *fftTmp1 = NULL;
    size_t    resNrow, resNcol, resNsec, tplSize, refSize, resSize;

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
            ERROR("ecorrm", "unsupported xcorr mode");
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
    xcorrNormTemplateECC(tplDataTmp, tplSize, mskData);

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
    array_delete(fftTmp1);
}

// instantiation
template
void ecorrm<float >(float  *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                    float  *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                    float  *resData,
                    float  *mskData,
                    eXCorrRes shape);

template
void ecorrm<double>(double *tplData, size_t tplNrow, size_t tplNcol, size_t tplNsec,
                    double *refData, size_t refNrow, size_t refNcol, size_t refNsec,
                    double *resData,
                    double *mskData,
                    eXCorrRes shape);

} // namespace gem
