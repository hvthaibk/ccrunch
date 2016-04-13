/***********************************************************************
 *  File:       xcorr_helper.cpp
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
#include "fftw_lock.hpp"

namespace gem {

/**********************************************
 * Helper functions for computing correlation
 *********************************************/

// modulation
template <typename T1, typename T2>
void pcorrModulateAndNormalize(const T1* const tplFFT,
                               const T1* const refFFT,
                                     T1* const resFFT,
                               size_t fftSize,
                               T2     nElements)
{
    assert(tplFFT != NULL);
    assert(refFFT != NULL);
    assert(resFFT != NULL);
    assert(fftSize   > 0);
    assert(nElements > 0);

    // obtain the cross power spectrum and then normalize
    #pragma omp parallel for
    for (size_t i = 0; i < fftSize; i++) {
        resFFT[i][0] = tplFFT[i][0] * refFFT[i][0] +
                       tplFFT[i][1] * refFFT[i][1];
        resFFT[i][1] = tplFFT[i][1] * refFFT[i][0] -
                       tplFFT[i][0] * refFFT[i][1];

        T2 absFFT = std::sqrt(resFFT[i][0]*resFFT[i][0] +
                              resFFT[i][1]*resFFT[i][1]);

        resFFT[i][0] /= absFFT * nElements;
        resFFT[i][1] /= absFFT * nElements;
    }
}

// instantiation
template
void pcorrModulateAndNormalize<fftwf_complex,float >(const fftwf_complex* const tplFFT,
                                                     const fftwf_complex* const refFFT,
                                                           fftwf_complex* const resFFT,
                                                     size_t fftSize,
                                                     float  nElements);
template
void pcorrModulateAndNormalize<fftw_complex ,double>(const fftw_complex*  const tplFFT,
                                                     const fftw_complex*  const refFFT,
                                                           fftw_complex*  const resFFT,
                                                     size_t fftSize,
                                                     double nElements);

// modulation
template <typename T1, typename T2>
void xcorrModulateAndNormalize(const T1* const tplFFT,
                               const T1* const refFFT,
                                     T1* const resFFT,
                               size_t fftSize,
                               T2     nElements)
{
    assert(tplFFT != NULL);
    assert(refFFT != NULL);
    assert(resFFT != NULL);
    assert(fftSize   > 0);
    assert(nElements > 0);

    // obtain the cross power spectrum and then normalize
    #pragma omp parallel for
    for (size_t i = 0; i < fftSize; i++) {
        resFFT[i][0] = (refFFT[i][0] * tplFFT[i][0] +
                        refFFT[i][1] * tplFFT[i][1]) / nElements;
        resFFT[i][1] = (refFFT[i][1] * tplFFT[i][0] -
                        refFFT[i][0] * tplFFT[i][1]) / nElements;
    }
}

// instantiation
template
void xcorrModulateAndNormalize<fftwf_complex,float >(const fftwf_complex* const tplFFT,
                                                     const fftwf_complex* const refFFT,
                                                           fftwf_complex* const resFFT,
                                                     size_t fftSize,
                                                     float  nElements);
template
void xcorrModulateAndNormalize<fftw_complex ,double>(const fftw_complex*  const tplFFT,
                                                     const fftw_complex*  const refFFT,
                                                           fftw_complex*  const resFFT,
                                                     size_t fftSize,
                                                     double nElements);

// template normalization
template <typename T>
void xcorrNormTemplateECC(T* const tplData, size_t tplSize)
{
    assert(tplData != NULL);
    assert(tplSize > 0);

    array_math_div(tplData, std::sqrt(array_reduce_sum2(tplData, tplSize)), tplSize);
}

// instantiation
template
void xcorrNormTemplateECC<float >(float*  const tplData, size_t tplSize);
template
void xcorrNormTemplateECC<double>(double* const tplData, size_t tplSize);

// template normalization
template <typename T>
void xcorrNormTemplateECC(      T* const tplData, size_t tplSize,
                          const T* const mskData)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplSize > 0);

    array_math_mul(tplData, mskData, tplSize);
    array_math_div(tplData, std::sqrt(array_reduce_sum2(tplData, tplSize)), tplSize);
}

// instantiation
template
void xcorrNormTemplateECC<float >(      float*  const tplData, size_t tplSize,
                                  const float*  const mskData);
template
void xcorrNormTemplateECC<double>(      double* const tplData, size_t tplSize,
                                  const double* const mskData);

// template normalization
template <typename T>
size_t xcorrNormTemplateNCC(T* const tplData, size_t tplSize)
{
    assert(tplData != NULL);
    assert(tplSize > 0);

    array_math_norm(tplData, tplSize);

    return tplSize;
}

// instantiation
template
size_t xcorrNormTemplateNCC<float >(float*  const tplData, size_t tplSize);
template
size_t xcorrNormTemplateNCC<double>(double* const tplData, size_t tplSize);

// template normalization
template <typename T>
size_t xcorrNormTemplateNCC(      T* const tplData, size_t tplSize,
                            const T* const mskData)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplSize > 0);

/*#ifndef NDEBUG
    T   arraySet[2] = { 0.0, 1.0 }, *arrayMem = NULL;
    array_new(arrayMem, tplSize);
    array_ismember(mskData, tplSize, arraySet, 2, arrayMem);
    assert(array_reduce_sum(arrayMem, tplSize) == tplSize);
    array_delete(arrayMem);
#endif*/

    size_t       *indxArray = NULL;
    size_t       i, mskArea;
    T            sum, tmp, mean, std;

    // allocate memory
    array_new(indxArray, tplSize);

    // indices of non-zero elements
    mskArea = 0;
    sum     = 0;
    for (i = 0; i < tplSize; i++) {
        if (((int) mskData[i]) > 0) {
            indxArray[mskArea++] = i;
            sum += tplData[i];
        }
    }

    // normalize the array
    mean = sum / (T) mskArea;

    sum = 0;
    for (i = 0; i < mskArea; i++) {
        tmp = tplData[indxArray[i]] - mean;
        sum += tmp*tmp;
    }
    std = std::sqrt(sum / (T) mskArea);

    for (i = 0; i < tplSize; i++) {
        tplData[i] = mskData[i] * (tplData[i] - mean) / std;
    }

    // deallocate memory
    array_delete(indxArray);

    return mskArea;
}

// instantiation
template
size_t xcorrNormTemplateNCC<float >(      float*  const tplData, size_t tplSize,
                                    const float*  const mskData);
template
size_t xcorrNormTemplateNCC<double>(      double* const tplData, size_t tplSize,
                                    const double* const mskData);

// template normalization
template <typename T>
void xcorrNormTemplateWCC1(T* const tplData, size_t tplSize,
                           T* const tplWght)
{
    assert(tplData != NULL);
    assert(tplWght != NULL);
    assert(tplSize > 0);

    T    *dataTmp = NULL;

    array_new(dataTmp, tplSize);

    // preprocess the weight
    array_math_div(tplWght, array_reduce_sum(tplWght, tplSize), tplSize);

    // normalize the template
    array_math_mul(dataTmp, tplData, tplWght, tplSize);
    array_math_sub(tplData, array_reduce_sum(dataTmp, tplSize), tplSize);

    array_math_sqr(dataTmp, tplData, tplSize);
    array_math_mul(dataTmp, tplWght, tplSize);

    array_math_sqrt(tplWght, tplSize);
    array_math_mul(tplData, tplWght, tplSize);
    array_math_div(tplData, std::sqrt(array_reduce_sum(dataTmp, tplSize)), tplSize);
    //array_math_sqr(tplWght, tplSize);

    array_delete(dataTmp);
}

// instantiation
template
void xcorrNormTemplateWCC1<float >(float*  const tplData, size_t tplSize,
                                   float*  const tplWght);
template
void xcorrNormTemplateWCC1<double>(double* const tplData, size_t tplSize,
                                   double* const tplWght);

// template normalization
template <typename T>
void xcorrNormTemplateWCC2(T* const tplData, size_t tplSize,
                           T* const tplWght,
                           T* const refWght)
{
    assert(tplData != NULL);
    assert(tplWght != NULL);
    assert(refWght != NULL);
    assert(tplSize > 0);

    // preprocess the weights
    array_math_div(tplWght, array_reduce_sum(tplWght, tplSize), tplSize);
    array_math_div(refWght, array_reduce_sum(refWght, tplSize), tplSize);

    // normalize the template
    T    *dataTmp = NULL, *dataNom = NULL;

    array_new(dataTmp, tplSize);
    array_new(dataNom, tplSize);

    array_math_mul(dataTmp, tplWght, tplData, tplSize);
    array_math_sub(tplData, array_reduce_sum(dataTmp, tplSize), tplSize);

    array_math_mul(dataTmp, tplWght, refWght, tplSize);
    array_math_sqrt(dataTmp, tplSize);
    array_math_mul(dataNom, dataTmp, tplData, tplSize);

    array_math_sqr(tplData, tplSize);
    array_math_mul(dataTmp, tplWght, tplData, tplSize);

    array_math_div(tplData, dataNom,
                   std::sqrt(array_reduce_sum(dataTmp, tplSize)),
                   tplSize);

    array_delete(dataTmp);
    array_delete(dataNom);
}

// instantiation
template
void xcorrNormTemplateWCC2<float >(float*  const tplData, size_t tplSize,
                                   float*  const tplWght,
                                   float*  const refWght);
template
void xcorrNormTemplateWCC2<double>(double* const tplData, size_t tplSize,
                                   double* const tplWght,
                                   double* const refWght);

// template normalization
template <typename T>
void xcorrNormTemplateWCC1(      T* const tplData, size_t tplSize,
                           const T* const mskData,
                                 T* const tplWght)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplWght != NULL);
    assert(tplSize > 0);

    T    *dataTmp = NULL;

    array_new(dataTmp, tplSize);

    // preprocess the weight
    array_math_mul(tplWght, mskData, tplSize);
    array_math_div(tplWght, array_reduce_sum(tplWght, tplSize), tplSize);

    // normalize the template
    array_math_mul(dataTmp, tplData, tplWght, tplSize);
    array_math_sub(tplData, array_reduce_sum(dataTmp, tplSize), tplSize);
    array_math_mul(tplData, mskData, tplSize);

    array_math_sqr(dataTmp, tplData, tplSize);
    array_math_mul(dataTmp, tplWght, tplSize);

    array_math_sqrt(tplWght, tplSize);
    array_math_mul(tplData, tplWght, tplSize);
    array_math_div(tplData, std::sqrt(array_reduce_sum(dataTmp, tplSize)), tplSize);
    //array_math_sqr(tplWght, tplSize);

    array_delete(dataTmp);
}

// instantiation
template
void xcorrNormTemplateWCC1<float >(      float*  const tplData, size_t tplSize,
                                   const float*  const mskData,
                                         float*  const tplWght);
template
void xcorrNormTemplateWCC1<double>(      double* const tplData, size_t tplSize,
                                   const double* const mskData,
                                         double* const tplWght);

// template normalization
template <typename T>
void xcorrNormTemplateWCC2(      T* const tplData, size_t tplSize,
                           const T* const mskData,
                                 T* const tplWght,
                                 T* const refWght)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplWght != NULL);
    assert(refWght != NULL);
    assert(tplSize > 0);

    for (size_t i = 0; i < tplSize; i++) {
        if (((int) mskData[i]) == 0) {
            tplData[i] = tplWght[i] = refWght[i] = 0;
        }
    }

    // preprocess the weights
    array_math_div(tplWght, array_reduce_sum(tplWght, tplSize), tplSize);
    array_math_div(refWght, array_reduce_sum(refWght, tplSize), tplSize);

    // normalize the template
    T    *dataTmp = NULL, *dataNom = NULL;

    array_new(dataTmp, tplSize);
    array_new(dataNom, tplSize);

    array_math_mul(dataTmp, tplWght, tplData, tplSize);
    array_math_sub(tplData, array_reduce_sum(dataTmp, tplSize), tplSize);

    array_math_mul(dataTmp, tplWght, refWght, tplSize);
    array_math_sqrt(dataTmp, tplSize);
    array_math_mul(dataNom, dataTmp, tplData, tplSize);

    array_math_sqr(tplData, tplSize);
    array_math_mul(dataTmp, tplWght, tplData, tplSize);

    array_math_div(tplData, dataNom,
                   std::sqrt(array_reduce_sum(dataTmp, tplSize)),
                   tplSize);

    array_delete(dataTmp);
    array_delete(dataNom);
}

// instantiation
template
void xcorrNormTemplateWCC2<float >(      float*  const tplData, size_t tplSize,
                                   const float*  const mskData,
                                         float*  const tplWght,
                                         float*  const refWght);
template
void xcorrNormTemplateWCC2<double>(      double* const tplData, size_t tplSize,
                                   const double* const mskData,
                                         double* const tplWght,
                                         double* const refWght);

// combine correlations
template <typename T>
void xcorrCombineResultCCC(      T* const resData,
                           const T* const corrData1, T weight1,
                           const T* const corrData2, T weight2,
                           const T* const corrData3, T weight3,
                           size_t resSize)
{
    assert(resData   != NULL);
    assert(corrData1 != NULL);
    assert(corrData2 != NULL);
    assert(corrData3 != NULL);
    assert(resSize > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        resData[i] = weight1 * corrData1[i] +
                     weight2 * corrData2[i] -
                     weight3 * corrData3[i];
    }
}

// instantiation
template
void xcorrCombineResultCCC<float >(      float*  const resData,
                                   const float*  const corrData1, float  weight1,
                                   const float*  const corrData2, float  weight2,
                                   const float*  const corrData3, float  weight3,
                                   size_t resSize);
template
void xcorrCombineResultCCC<double>(      double* const resData,
                                   const double* const corrData1, double weight1,
                                   const double* const corrData2, double weight2,
                                   const double* const corrData3, double weight3,
                                   size_t resSize);

// combine correlations
template <typename T>
void xcorrCombineResultECC(      T* const resData,
                           const T* const fftTmp1,
                           size_t resSize)
{
    assert(resData != NULL);
    assert(fftTmp1 != NULL);
    assert(resSize > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        resData[i] /= std::sqrt(fftTmp1[i]);
    }
}

// instantiation
template
void xcorrCombineResultECC<float >(      float*  const resData,
                                   const float*  const fftTmp1,
                                   size_t resSize);
template
void xcorrCombineResultECC<double>(      double* const resData,
                                   const double* const fftTmp1,
                                   size_t resSize);

// combine correlations
template <typename T>
void xcorrCombineResultNCC(      T* const resData,
                           const T* const fftTmp1,
                           const T* const fftTmp2,
                           size_t mskArea,
                           size_t resSize)
{
    assert(resData != NULL);
    assert(fftTmp1 != NULL);
    assert(fftTmp2 != NULL);
    assert(mskArea > 0);
    assert(resSize > 0);

    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        resData[i] /= std::sqrt(fftTmp2[i] * (T) mskArea - fftTmp1[i]*fftTmp1[i]);
    }
}

// instantiation
template
void xcorrCombineResultNCC<float >(      float*  const resData,
                                   const float*  const fftTmp1,
                                   const float*  const fftTmp2,
                                   size_t mskArea,
                                   size_t resSize);
template
void xcorrCombineResultNCC<double>(      double* const resData,
                                   const double* const fftTmp1,
                                   const double* const fftTmp2,
                                   size_t mskArea,
                                   size_t resSize);

// combine correlations
template <typename T>
void xcorrCombineResultWCC(      T* const resData,
                           const T* const fftTmp1,
                           const T* const fftTmp2,
                           T      tplSum,
                           size_t resSize)
{
    for (size_t i = 0; i < resSize; i++) {
        resData[i] -= tplSum * fftTmp1[i];
        resData[i] /= std::sqrt(fftTmp2[i] - fftTmp1[i]*fftTmp1[i]);
    }
}

// instantiation
template
void xcorrCombineResultWCC<float >(      float*  const resData,
                                   const float*  const fftTmp1,
                                   const float*  const fftTmp2,
                                   float  tplSum,
                                   size_t resSize);
template
void xcorrCombineResultWCC<double>(      double* const resData,
                                   const double* const fftTmp1,
                                   const double* const fftTmp2,
                                   double tplSum,
                                   size_t resSize);

// combine correlations
template <typename T>
void xcorrCombineResultNReg(      T* const resCorr,
                            const T* const resOlap,
                            const T* const numTpl,
                            const T* const numRef,
                            const T* const denTpl,
                            const T* const denRef,
                            size_t resSize,
                            T      overlapRatio)
{
    assert(resCorr != NULL);
    assert(resOlap != NULL);
    assert(numTpl  != NULL);
    assert(numRef  != NULL);
    assert(denTpl  != NULL);
    assert(denRef  != NULL);
    assert(resSize > 0);

    T     threshVal = overlapRatio * array_reduce_max(resOlap, resSize);

    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        if (resOlap[i] >= threshVal) {
            resCorr[i] = (resCorr[i] - numTpl[i]*numRef[i]/resOlap[i])
                         / std::sqrt((denTpl[i] - pow2(numTpl[i])/resOlap[i]) *
                                     (denRef[i] - pow2(numRef[i])/resOlap[i]));
        }
        else {
            resCorr[i] = 0;
        }
    }
}

// instantiation
template
void xcorrCombineResultNReg<float >(      float*  const resCorr,
                                    const float*  const resOlap,
                                    const float*  const numTpl,
                                    const float*  const numRef,
                                    const float*  const denTpl,
                                    const float*  const denRef,
                                    size_t resSize,
                                    float  overlapRatio);
template
void xcorrCombineResultNReg<double>(      double* const resCorr,
                                    const double* const resOlap,
                                    const double* const numTpl,
                                    const double* const numRef,
                                    const double* const denTpl,
                                    const double* const denRef,
                                    size_t resSize,
                                    double overlapRatio);

// merge result
template <typename T>
void xcorrMergeResult(const T*      const resData,
                            T*      const resDataAbsMax,
                            size_t* const resDataMaxInd,
                      size_t        resSize,
                      size_t        indx,
                      eXCorrMerge   bAbs)
{
    assert(resData       != NULL);
    assert(resDataAbsMax != NULL);
    assert(resDataMaxInd != NULL);
    assert(resSize > 0);

    switch (bAbs) {
        case XCORR_MERGE_NEGATIVE:
            #pragma omp parallel for
            for (size_t i = 0; i < resSize; i++) {
                if (std::abs(resData[i]) > std::abs(resDataAbsMax[i])) {
                    resDataAbsMax[i] = resData[i];
                    resDataMaxInd[i] = indx;
                }
            }
            break;
        case XCORR_MERGE_POSITIVE:
            #pragma omp parallel for
            for (size_t i = 0; i < resSize; i++) {
                if (resData[i] > resDataAbsMax[i]) {
                    resDataAbsMax[i] = resData[i];
                    resDataMaxInd[i] = indx;
                }
            }
            break;
        default:
            ERROR("xcorrMergeResult", "unsupported merging mode");
    }
}

// instantiation
template
void xcorrMergeResult<float >(const float*  const resData,
                                    float*  const resDataAbsMax,
                                    size_t* const resDataMaxInd,
                              size_t        resSize,
                              size_t        indx,
                              eXCorrMerge   bAbs);
template
void xcorrMergeResult<double>(const double* const resData,
                                    double* const resDataAbsMax,
                                    size_t* const resDataMaxInd,
                              size_t        resSize,
                              size_t        indx,
                              eXCorrMerge   bAbs);

// merge result
template <typename T>
void xcorrMergeResultGlobal(const T*      const resDataAbsMax,
                            const size_t* const resDataMaxInd,
                                  T*      const resDataAbsMaxGlobal,
                                  size_t* const resDataMaxIndGlobal,
                            size_t        resSize,
                            eXCorrMerge   bAbs)
{
    assert(resDataAbsMax       != NULL);
    assert(resDataMaxInd       != NULL);
    assert(resDataAbsMaxGlobal != NULL);
    assert(resDataMaxIndGlobal != NULL);
    assert(resSize > 0);

    mutexMergeResult.lock();

    switch (bAbs) {
        case XCORR_MERGE_NEGATIVE:
            #pragma omp parallel for
            for (size_t i = 0; i < resSize; i++) {
                if (std::abs(resDataAbsMax[i]) > std::abs(resDataAbsMaxGlobal[i])) {
                    resDataAbsMaxGlobal[i] = resDataAbsMax[i];
                    resDataMaxIndGlobal[i] = resDataMaxInd[i];
                }
            }
            break;
        case XCORR_MERGE_POSITIVE:
            #pragma omp parallel for
            for (size_t i = 0; i < resSize; i++) {
                if (resDataAbsMax[i] > resDataAbsMaxGlobal[i]) {
                    resDataAbsMaxGlobal[i] = resDataAbsMax[i];
                    resDataMaxIndGlobal[i] = resDataMaxInd[i];
                }
            }
            break;
        default:
            ERROR("xcorrMergeResultGlobal", "unsupported merging mode");
    }

    mutexMergeResult.unlock();
}

// instantiation
template
void xcorrMergeResultGlobal<float >(const float*  const resDataAbsMax,
                                    const size_t* const resDataMaxInd,
                                          float*  const resDataAbsMaxGlobal,
                                          size_t* const resDataMaxIndGlobal,
                                    size_t        resSize,
                                    eXCorrMerge   bAbs);
template
void xcorrMergeResultGlobal<double>(const double* const resDataAbsMax,
                                    const size_t* const resDataMaxInd,
                                          double* const resDataAbsMaxGlobal,
                                          size_t* const resDataMaxIndGlobal,
                                    size_t        resSize,
                                    eXCorrMerge   bAbs);

} // namespace gem
