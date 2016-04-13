/***********************************************************************
 *  File:       xcorr_helper.cu
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

void printErrMsgCUFFT(cufftResult_t err)
{
    switch (err) {
        case CUFFT_SUCCESS:
            std::cout << "CUFFT operation is successful" << std::endl;
            break;
        case CUFFT_INVALID_PLAN:
            std::cout << "CUFFT is passed an invalid plan handle" << std::endl;
            break;
        case CUFFT_ALLOC_FAILED:
            std::cout << "CUFFT failed to allocate GPU memory" << std::endl;
            break;
        case CUFFT_INVALID_TYPE:
            std::cout << "CUFFT: an unsupported type" << std::endl;
            break;
        case CUFFT_INVALID_VALUE:
            std::cout << "CUFFT: a bad memory pointer" << std::endl;
            break;
        case CUFFT_INTERNAL_ERROR:
            std::cout << "CUFFT internal driver errors" << std::endl;
            break;
        case CUFFT_EXEC_FAILED:
            std::cout << "CUFFT failed to execute an FFT on the GPU" << std::endl;
            break;
        case CUFFT_SETUP_FAILED:
            std::cout << "CUFFT library failed to initialize" << std::endl;
            break;
        case CUFFT_INVALID_SIZE:
            std::cout << "CUFFT: an unsupported FFT size" << std::endl;
            break;
        case CUFFT_UNALIGNED_DATA:
            std::cout << "CUFFT: unaligned data" << std::endl;
            break;
        case CUFFT_INCOMPLETE_PARAMETER_LIST:
            std::cout << "CUFFT: incomplete parameter list" << std::endl;
            break;
        case CUFFT_INVALID_DEVICE:
            std::cout << "CUFFT: invalid device" << std::endl;
            break;
        case CUFFT_PARSE_ERROR:
            std::cout << "CUFFT: parse error" << std::endl;
            break;
        case CUFFT_NO_WORKSPACE:
            std::cout << "CUFFT: no workspace" << std::endl;
            break;
        /*case CUFFT_NOT_IMPLEMENTED:
            std::cout << "CUFFT: not implemented" << std::endl;
            break;
        case CUFFT_LICENSE_ERROR:
            std::cout << "CUFFT: license error" << std::endl;
            break;*/
        default:
            std::cout << "CUFFT: unrecognized error" << std::endl;
    }
}

/***********************************************************************
 * Helper functions for computing correlation
 **********************************************************************/

// modulation
template <typename T1, typename T2>
void cuda_pcorrModulateAndNormalize(const T1* const tplFFT,
                                    const T1* const refFFT,
                                          T1* const resFFT,
                                    size_t fftSize,
                                    T2     nElements);

// modulation
template <typename T1, typename T2> __global__
void dev_xcorrModulateAndNormalize(const T1* const tplFFT,
                                   const T1* const refFFT,
                                         T1* const resFFT,
                                   size_t fftSize,
                                   T2     nElements)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < fftSize) {
        resFFT[tid].x = (refFFT[tid].x * tplFFT[tid].x +
                         refFFT[tid].y * tplFFT[tid].y) / nElements;
        resFFT[tid].y = (refFFT[tid].y * tplFFT[tid].x -
                         refFFT[tid].x * tplFFT[tid].y) / nElements;
    }
}

template <typename T1, typename T2>
void cuda_xcorrModulateAndNormalize(const T1* const tplFFT,
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

    dev_xcorrModulateAndNormalize
        <<<iDivUp(fftSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplFFT, refFFT, resFFT, fftSize, nElements);
}

// instantiation
template
void cuda_xcorrModulateAndNormalize<cuFloatComplex,float >
        (const cuFloatComplex*  const tplFFT,
         const cuFloatComplex*  const refFFT,
               cuFloatComplex*  const resFFT,
         size_t          fftSize,
         float           nElements);

template
void cuda_xcorrModulateAndNormalize<cuDoubleComplex,double>
        (const cuDoubleComplex* const tplFFT,
         const cuDoubleComplex* const refFFT,
               cuDoubleComplex* const resFFT,
         size_t          fftSize,
         double          nElements);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateECC(T* const tplData, size_t tplSize)
{
    assert(tplData != NULL);
    assert(tplSize > 0);

    cuda_array_math_div(tplData, std::sqrt(cuda_array_reduce_sum2(tplData, tplSize)), tplSize);
}

// instantiation
template
void cuda_xcorrNormTemplateECC<float >(float*  const tplData, size_t tplSize);
template
void cuda_xcorrNormTemplateECC<double>(double* const tplData, size_t tplSize);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateECC(      T* const tplData, size_t tplSize,
                                 const T* const mskData)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplSize > 0);

    cuda_array_math_mul(tplData, mskData, tplSize);
    cuda_array_math_div(tplData, std::sqrt(cuda_array_reduce_sum2(tplData, tplSize)), tplSize);
}

// instantiation
template
void cuda_xcorrNormTemplateECC<float >(      float*  const tplData, size_t tplSize,
                                       const float*  const mskData);
template
void cuda_xcorrNormTemplateECC<double>(      double* const tplData, size_t tplSize,
                                       const double* const mskData);

// template normalization
template <typename T>
size_t cuda_xcorrNormTemplateNCC(T* const tplData, size_t tplSize)
{
    assert(tplData != NULL);
    assert(tplSize > 0);

    cuda_array_math_norm(tplData, tplSize);

    return tplSize;
}

// instantiation
template
size_t cuda_xcorrNormTemplateNCC<float >(float*  const tplData, size_t tplSize);
template
size_t cuda_xcorrNormTemplateNCC<double>(double* const tplData, size_t tplSize);

// template normalization
template <typename T>
size_t cuda_xcorrNormTemplateNCC(      T* const tplData, size_t tplSize,
                                 const T* const mskData)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplSize > 0);

    size_t    mskArea;
    T         sum, mean, std;

    // mean
    mskArea = (size_t) cuda_array_reduce_sum(mskData, tplSize);
    sum     = cuda_array_mask_sum(tplData, mskData, tplSize);
    mean    = sum / (T) mskArea;

    // std
    sum = cuda_array_mask_subsqrsum(tplData, mean, mskData, tplSize);
    std = std::sqrt(sum / (T) mskArea);

    // normalize
    cuda_array_mask_subdivmul(tplData, mean, std, mskData, tplSize);

    return mskArea;
}

// instantiation
template
size_t cuda_xcorrNormTemplateNCC<float >(      float*  const tplData, size_t tplSize,
                                         const float*  const mskData);
template
size_t cuda_xcorrNormTemplateNCC<double>(      double* const tplData, size_t tplSize,
                                         const double* const mskData);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateWCC1(T* const tplData, size_t tplSize,
                                T* const tplWght)
{
    assert(tplData != NULL);
    assert(tplWght != NULL);
    assert(tplSize > 0);

    T    *dataTmp = NULL;

    cuda_arrayDev_new(dataTmp, tplSize);

    // preprocess the weight
    cuda_array_math_div(tplWght, cuda_array_reduce_sum(tplWght, tplSize), tplSize);

    // normalize the template
    cuda_array_math_mul(dataTmp, tplData, tplWght, tplSize);
    cuda_array_math_sub(tplData, cuda_array_reduce_sum(dataTmp, tplSize), tplSize);

    cuda_array_math_sqr(dataTmp, tplData, tplSize);
    cuda_array_math_mul(dataTmp, tplWght, tplSize);

    cuda_array_math_sqrt(tplWght, tplSize);
    cuda_array_math_mul(tplData, tplWght, tplSize);
    cuda_array_math_div(tplData, std::sqrt(cuda_array_reduce_sum(dataTmp, tplSize)), tplSize);
    //cuda_array_math_sqr(tplWght, tplSize);

    cuda_arrayDev_delete(dataTmp);
}

// instantiation
template
void cuda_xcorrNormTemplateWCC1<float >(float*  const tplData, size_t tplSize,
                                        float*  const tplWght);
template
void cuda_xcorrNormTemplateWCC1<double>(double* const tplData, size_t tplSize,
                                        double* const tplWght);

// template normalization
template <typename T> __global__
void dev_array_weight_zero(      T* const tplData, size_t tplSize,
                           const T* const mskData,
                                 T* const tplWght,
                                 T* const refWght)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < tplSize && ((int) mskData[i]) ==  0) {
        tplData[i] = tplWght[i] = refWght[i] = 0;
    }
}

template <typename T> __global__
void dev_array_weight_nom(T* const tplData, size_t tplSize,
                          T* const tplWght,
                          T* const refWght,
                          T* const dataNom)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < tplSize) {
        dataNom[i] = sqrtf(tplWght[i]*refWght[i]) * tplData[i];
    }
}

template <> __global__
void dev_array_weight_nom<double>(double* const tplData, size_t tplSize,
                                  double* const tplWght,
                                  double* const refWght,
                                  double* const dataNom)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < tplSize) {
        dataNom[i] = sqrt(tplWght[i]*refWght[i]) * tplData[i];
    }
}

template <typename T> __global__
void dev_array_weight_tmp(T* const tplData, size_t tplSize,
                          T* const tplWght,
                          T* const dataTmp)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < tplSize) {
        dataTmp[i] = tplData[i] * tplData[i] * tplWght[i];
    }
}

template <typename T>
void cuda_xcorrNormTemplateWCC2(T* const tplData, size_t tplSize,
                                T* const tplWght,
                                T* const refWght)
{
    assert(tplData != NULL);
    assert(tplWght != NULL);
    assert(refWght != NULL);
    assert(tplSize > 0);

    // preprocess the weights
    cuda_array_math_div(tplWght, cuda_array_reduce_sum(tplWght, tplSize), tplSize);
    cuda_array_math_div(refWght, cuda_array_reduce_sum(refWght, tplSize), tplSize);

    // normalize the template
    T    *dataTmp = NULL, *dataNom = NULL;

    cuda_arrayDev_new(dataTmp, tplSize);
    cuda_arrayDev_new(dataNom, tplSize);

    cuda_array_math_mul(dataTmp, tplWght, tplData, tplSize);
    cuda_array_math_sub(tplData, cuda_array_reduce_sum(dataTmp, tplSize), tplSize);

    dev_array_weight_nom
        <<<iDivUp(tplSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplData, tplSize, tplWght, refWght, dataNom);

    dev_array_weight_tmp
        <<<iDivUp(tplSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplData, tplSize, tplWght, dataTmp);

    cuda_array_math_div(tplData, dataNom,
                        std::sqrt(cuda_array_reduce_sum(dataTmp, tplSize)),
                        tplSize);

    cuda_arrayDev_delete(dataTmp);
    cuda_arrayDev_delete(dataNom);
}

// instantiation
template
void cuda_xcorrNormTemplateWCC2<float >(float*  const tplData, size_t tplSize,
                                        float*  const tplWght,
                                        float*  const refWght);
template
void cuda_xcorrNormTemplateWCC2<double>(double* const tplData, size_t tplSize,
                                        double* const tplWght,
                                        double* const refWght);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateWCC1(      T* const tplData, size_t tplSize,
                                const T* const mskData,
                                      T* const tplWght)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplWght != NULL);
    assert(tplSize > 0);

    T    *dataTmp = NULL;

    cuda_arrayDev_new(dataTmp, tplSize);

    // preprocess the weight
    cuda_array_math_mul(tplWght, mskData, tplSize);
    cuda_array_math_div(tplWght, cuda_array_reduce_sum(tplWght, tplSize), tplSize);

    // normalize the template
    cuda_array_math_mul(dataTmp, tplData, tplWght, tplSize);
    cuda_array_math_sub(tplData, cuda_array_reduce_sum(dataTmp, tplSize), tplSize);
    cuda_array_math_mul(tplData, mskData, tplSize);

    cuda_array_math_sqr(dataTmp, tplData, tplSize);
    cuda_array_math_mul(dataTmp, tplWght, tplSize);

    cuda_array_math_sqrt(tplWght, tplSize);
    cuda_array_math_mul(tplData, tplWght, tplSize);
    cuda_array_math_div(tplData, std::sqrt(cuda_array_reduce_sum(dataTmp, tplSize)), tplSize);
    //cuda_array_math_sqr(tplWght, tplSize);

    cuda_arrayDev_delete(dataTmp);
}

// instantiation
template
void cuda_xcorrNormTemplateWCC1<float >(      float*  const tplData, size_t tplSize,
                                        const float*  const mskData,
                                              float*  const tplWght);
template
void cuda_xcorrNormTemplateWCC1<double>(      double* const tplData, size_t tplSize,
                                        const double* const mskData,
                                              double* const tplWght);

// template normalization
template <typename T>
void cuda_xcorrNormTemplateWCC2(      T* const tplData, size_t tplSize,
                                const T* const mskData,
                                      T* const tplWght,
                                      T* const refWght)
{
    assert(tplData != NULL);
    assert(mskData != NULL);
    assert(tplWght != NULL);
    assert(refWght != NULL);
    assert(tplSize > 0);

    dev_array_weight_zero
        <<<iDivUp(tplSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplData, tplSize, mskData, tplWght, refWght);

    // preprocess the weights
    cuda_array_math_div(tplWght, cuda_array_reduce_sum(tplWght, tplSize), tplSize);
    cuda_array_math_div(refWght, cuda_array_reduce_sum(refWght, tplSize), tplSize);

    // normalize the template
    T    *dataTmp = NULL, *dataNom = NULL;

    cuda_arrayDev_new(dataTmp, tplSize);
    cuda_arrayDev_new(dataNom, tplSize);

    cuda_array_math_mul(dataTmp, tplWght, tplData, tplSize);
    cuda_array_math_sub(tplData, cuda_array_reduce_sum(dataTmp, tplSize), tplSize);

    dev_array_weight_nom
        <<<iDivUp(tplSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplData, tplSize, tplWght, refWght, dataNom);

    dev_array_weight_tmp
        <<<iDivUp(tplSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (tplData, tplSize, tplWght, dataTmp);

    cuda_array_math_div(tplData, dataNom,
                        std::sqrt(cuda_array_reduce_sum(dataTmp, tplSize)),
                        tplSize);

    cuda_arrayDev_delete(dataTmp);
    cuda_arrayDev_delete(dataNom);
}

// instantiation
template
void cuda_xcorrNormTemplateWCC2<float >(      float*  const tplData, size_t tplSize,
                                        const float*  const mskData,
                                              float*  const tplWght,
                                              float*  const refWght);
template
void cuda_xcorrNormTemplateWCC2<double>(      double* const tplData, size_t tplSize,
                                        const double* const mskData,
                                              double* const tplWght,
                                              double* const refWght);

// combine correlations
template <typename T> __global__
void dev_xcorrCombineResultCCC(      T* const resData,
                               const T* const corrData1, T weight1,
                               const T* const corrData2, T weight2,
                               const T* const corrData3, T weight3,
                               size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] = weight1 * corrData1[tid] +
                       weight2 * corrData2[tid] -
                       weight3 * corrData3[tid];
    }
}

template <typename T>
void cuda_xcorrCombineResultCCC(      T* const resData,
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

    dev_xcorrCombineResultCCC
        <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (resData, corrData1, weight1, corrData2, weight2, corrData3, weight3, resSize);
}

// instantiation
template
void cuda_xcorrCombineResultCCC<float >(      float*  const resData,
                                        const float*  const corrData1, float  weight1,
                                        const float*  const corrData2, float  weight2,
                                        const float*  const corrData3, float  weight3,
                                        size_t resSize);
template
void cuda_xcorrCombineResultCCC<double>(      double* const resData,
                                        const double* const corrData1, double weight1,
                                        const double* const corrData2, double weight2,
                                        const double* const corrData3, double weight3,
                                        size_t resSize);

// combine correlations
template <typename T> __global__
void dev_xcorrCombineResultECC(      T* const resData,
                               const T* const fftTmp1,
                               size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] /= std::sqrt(fftTmp1[tid]);
    }
}

template <typename T>
void cuda_xcorrCombineResultECC(      T* const resData,
                                const T* const fftTmp1,
                                size_t resSize)
{
    assert(resData != NULL);
    assert(fftTmp1 != NULL);
    assert(resSize > 0);

    dev_xcorrCombineResultECC
        <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (resData, fftTmp1, resSize);
}

// instantiation
template
void cuda_xcorrCombineResultECC<float >(      float*  const resData,
                                        const float*  const fftTmp1,
                                        size_t resSize);
template
void cuda_xcorrCombineResultECC<double>(      double* const resData,
                                        const double* const fftTmp1,
                                        size_t resSize);

// combine correlations
template <typename T> __global__
void dev_xcorrCombineResultNCC(      T* const resData,
                               const T* const fftTmp1,
                               const T* const fftTmp2,
                               size_t mskArea,
                               size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] /= sqrtf(fftTmp2[tid] * (T) mskArea -
                              fftTmp1[tid] * fftTmp1[tid]);
    }
}

template <> __global__
void dev_xcorrCombineResultNCC<double>(      double* const resData,
                                       const double* const fftTmp1,
                                       const double* const fftTmp2,
                                       size_t mskArea,
                                       size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] /= sqrt(fftTmp2[tid] * (double) mskArea -
                             fftTmp1[tid] * fftTmp1[tid]);
    }
}

template <typename T>
void cuda_xcorrCombineResultNCC(      T* const resData,
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

    dev_xcorrCombineResultNCC
        <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (resData, fftTmp1, fftTmp2, mskArea, resSize);
}

// instantiation
template
void cuda_xcorrCombineResultNCC<float >(      float*  const resData,
                                        const float*  const fftTmp1,
                                        const float*  const fftTmp2,
                                        size_t mskArea,
                                        size_t resSize);
template
void cuda_xcorrCombineResultNCC<double>(      double* const resData,
                                        const double* const fftTmp1,
                                        const double* const fftTmp2,
                                        size_t mskArea,
                                        size_t resSize);

// combine correlations
template <typename T> __global__
void dev_xcorrCombineResultWCC(      T* const resData,
                               const T* const fftTmp1,
                               const T* const fftTmp2,
                               T      tplSum,
                               size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] -= tplSum * fftTmp1[tid];
        resData[tid] /= sqrtf(fftTmp2[tid] - fftTmp1[tid]*fftTmp1[tid]);
    }
}

template <> __global__
void dev_xcorrCombineResultWCC<double>(      double* const resData,
                                       const double* const fftTmp1,
                                       const double* const fftTmp2,
                                       double tplSum,
                                       size_t resSize)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        resData[tid] -= tplSum * fftTmp1[tid];
        resData[tid] /= sqrt(fftTmp2[tid] - fftTmp1[tid]*fftTmp1[tid]);
    }
}

template <typename T>
void cuda_xcorrCombineResultWCC(      T* const resData,
                                const T* const fftTmp1,
                                const T* const fftTmp2,
                                T      tplSum,
                                size_t resSize)

{
    dev_xcorrCombineResultWCC
        <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (resData, fftTmp1, fftTmp2, tplSum, resSize);
}

// instantiation
template
void cuda_xcorrCombineResultWCC<float >(      float*  const resData,
                                        const float*  const fftTmp1,
                                        const float*  const fftTmp2,
                                        float  tplSum,
                                        size_t resSize);
template
void cuda_xcorrCombineResultWCC<double>(      double* const resData,
                                        const double* const fftTmp1,
                                        const double* const fftTmp2,
                                        double tplSum,
                                        size_t resSize);

// combine correlations
template <typename T> __global__
void dev_xcorrCombineResultNReg(      T* const resCorr,
                                const T* const resOlap,
                                const T* const numTpl,
                                const T* const numRef,
                                const T* const denTpl,
                                const T* const denRef,
                                size_t resSize,
                                T      threshVal)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        if (resOlap[tid] >= threshVal) {
            resCorr[tid] = (resCorr[tid] - numTpl[tid]*numRef[tid]/resOlap[tid])
                           / std::sqrt((denTpl[tid] - numTpl[tid]*numTpl[tid]/resOlap[tid]) *
                                       (denRef[tid] - numRef[tid]*numRef[tid]/resOlap[tid]));
        }
        else {
            resCorr[tid] = 0;
        }
    }
}

template <typename T>
void cuda_xcorrCombineResultNReg(      T* const resCorr,
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

    T     threshVal = overlapRatio * cuda_array_reduce_max(resOlap, resSize);

    dev_xcorrCombineResultNReg
        <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (resCorr,  resOlap, numTpl, numRef, denTpl, denRef, resSize, threshVal);
}

// instantiation
template
void cuda_xcorrCombineResultNReg<float >(      float*  const resCorr,
                                         const float*  const resOlap,
                                         const float*  const numTpl,
                                         const float*  const numRef,
                                         const float*  const denTpl,
                                         const float*  const denRef,
                                         size_t resSize,
                                         float  overlapRatio);
template
void cuda_xcorrCombineResultNReg<double>(      double* const resCorr,
                                         const double* const resOlap,
                                         const double* const numTpl,
                                         const double* const numRef,
                                         const double* const denTpl,
                                         const double* const denRef,
                                         size_t resSize,
                                         double overlapRatio);

// merge result
template <typename T> __global__
void dev_xcorrMergeResultAbs(const T*      const resData,
                                   T*      const resDataAbsMax,
                                   size_t* const resDataMaxInd,
                             size_t resSize,
                             size_t indx)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        if (fabsf(resData[tid]) > fabsf(resDataAbsMax[tid])) {
            resDataAbsMax[tid] = resData[tid];
            resDataMaxInd[tid] = indx;
        }
    }
}

template <> __global__
void dev_xcorrMergeResultAbs<double>(const double* const resData,
                                           double* const resDataAbsMax,
                                           size_t* const resDataMaxInd,
                                     size_t resSize,
                                     size_t indx)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        if (fabs(resData[tid]) > fabs(resDataAbsMax[tid])) {
            resDataAbsMax[tid] = resData[tid];
            resDataMaxInd[tid] = indx;
        }
    }
}

template <typename T> __global__
void dev_xcorrMergeResult(const T*      const resData,
                                T*      const resDataAbsMax,
                                size_t* const resDataMaxInd,
                          size_t resSize,
                          size_t indx)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < resSize) {
        if (resData[tid] > resDataAbsMax[tid]) {
            resDataAbsMax[tid] = resData[tid];
            resDataMaxInd[tid] = indx;
        }
    }
}

template <typename T>
void cuda_xcorrMergeResult(const T*      const resData,
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
            dev_xcorrMergeResultAbs
                <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
                (resData,
                 resDataAbsMax,
                 resDataMaxInd,
                 resSize, indx);
            break;
        case XCORR_MERGE_POSITIVE:
            dev_xcorrMergeResult
                <<<iDivUp(resSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
                (resData,
                 resDataAbsMax,
                 resDataMaxInd,
                 resSize, indx);
            break;
        default:
            ERROR("cuda_xcorrMergeResult", "unsupported merging mode");
    }
}

// instantiation
template
void cuda_xcorrMergeResult<float >(const float*  const resData,
                                         float*  const resDataAbsMax,
                                         size_t* const resDataMaxInd,
                                   size_t        resSize,
                                   size_t        indx,
                                   eXCorrMerge   bAbs);
template
void cuda_xcorrMergeResult<double>(const double* const resData,
                                         double* const resDataAbsMax,
                                         size_t* const resDataMaxInd,
                                   size_t        resSize,
                                   size_t        indx,
                                   eXCorrMerge   bAbs);

} // namespace gem
