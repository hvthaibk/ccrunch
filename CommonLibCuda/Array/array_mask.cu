/***********************************************************************
 *  File:     array_mask.cu
 *
 *  Purpose:  Implementation of array-related functions
 *
 *  Author:   Thai V. Hoang
 *
 *  Contact:  hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.cuh"

namespace gem {

/*****************************************
 * Mask
 ****************************************/

// 1D copy
template <typename T> __global__
void dev_array_mask_copy(const T* const array,  size_t    nRow,
                         const T* const mask,   size_t    nRowMsk,
                               T* const output, ptrdiff_t nRowOff)
{
    ptrdiff_t    iRowMsk = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowMsk < nRowMsk) {
        if (mask[iRowMsk] > 0) {
            assert(iRowMsk + nRowOff >= 0);
            assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);

            output[iRowMsk] = array[iRowMsk+nRowOff];
        }
        else {
            output[iRowMsk] = 0;
        }
    }
}

template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,
                          const T* const mask,   size_t    nRowMsk,
                                T* const output, ptrdiff_t nRowOff)
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);

    dev_array_mask_copy
        <<<iDivUp(nRowMsk, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array,  nRow,
         mask,   nRowMsk,
         output, nRowOff);
}

// instantiation
template
void cuda_array_mask_copy<float >(const float*  const array,  size_t    nRow,
                                  const float*  const mask,   size_t    nRowMsk,
                                        float*  const output, ptrdiff_t nRowOff);
template
void cuda_array_mask_copy<double>(const double* const array,  size_t    nRow,
                                  const double* const mask,   size_t    nRowMsk,
                                        double* const output, ptrdiff_t nRowOff);

// 2D copy
template <typename T> __global__
void dev_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,
                         const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                               T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    ptrdiff_t    iRowMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iColMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iRowMsk < nRowMsk && iColMsk < nColMsk) {
        iMsk = iRowMsk*nColMsk + iColMsk;

        if (mask[iMsk] > 0) {
            assert(iRowMsk + nRowOff >= 0);
            assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
            assert(iColMsk + nColOff >= 0);
            assert(iColMsk + nColOff <  (ptrdiff_t) nCol);

            iSrc = (iRowMsk+nRowOff)*nCol + (iColMsk+nColOff);
            output[iMsk] = array[iSrc];
        }
        else {
            output[iMsk] = 0;
        }
    }
}

template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,
                          const T* const mask,   size_t    nRowMsk, size_t    nColMsk,
                                T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff)
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_mask_copy<<<dimGrid, dimBlock>>>
        (array,  nRow,    nCol,
         mask,   nRowMsk, nColMsk,
         output, nRowOff, nColOff);
}

// instantiation
template
void cuda_array_mask_copy<float >(const float*  const array,  size_t    nRow,    size_t    nCol,
                                  const float*  const mask,   size_t    nRowMsk, size_t    nColMsk,
                                        float*  const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);
template
void cuda_array_mask_copy<double>(const double* const array,  size_t    nRow,    size_t    nCol,
                                  const double* const mask,   size_t    nRowMsk, size_t    nColMsk,
                                        double* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff);

// 3D copy
template <typename T> __global__
void dev_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                         const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                               T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    ptrdiff_t    iRowMsk;
    ptrdiff_t    iColMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iColMsk < nColMsk && iSecMsk < nSecMsk) {
        for (iRowMsk = 0; iRowMsk < nRowMsk; iRowMsk++) {
            iMsk = iRowMsk*nColMsk*nSecMsk + iColMsk*nSecMsk + iSecMsk;

            if (mask[iMsk] > 0) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                assert(iColMsk + nColOff >= 0);
                assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
                assert(iSecMsk + nSecOff >= 0);
                assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

                iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
                output[iMsk] = array[iSrc];
            }
            else {
                output[iMsk] = 0;
            }
        }
    }
#else
    ptrdiff_t    iRowMsk = blockDim.z * blockIdx.z + threadIdx.z;
    ptrdiff_t    iColMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iRowMsk < nRowMsk && iColMsk < nColMsk && iSecMsk < nSecMsk) {
        iMsk = iRowMsk*nColMsk*nSecMsk + iColMsk*nSecMsk + iSecMsk;

        if (mask[iMsk] > 0) {
            assert(iRowMsk + nRowOff >= 0);
            assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
            assert(iColMsk + nColOff >= 0);
            assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
            assert(iSecMsk + nSecOff >= 0);
            assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

            iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
            output[iMsk] = array[iSrc];
        }
        else {
            output[iMsk] = 0;
        }
    }
#endif
}

template <typename T>
void cuda_array_mask_copy(const T* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                          const T* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                T* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff)
{
    assert(array  != NULL);
    assert(mask   != NULL);
    assert(output != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);
    assert(nSec > 0 && nSecMsk > 0);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y),
                     iDivUp(nRow, dimBlock.z));
#endif

    dev_array_mask_copy<<<dimGrid, dimBlock>>>
        (array,  nRow,    nCol,    nSec,
         mask,   nRowMsk, nColMsk, nSecMsk,
         output, nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_array_mask_copy<float >(const float*  const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                                  const float*  const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                        float*  const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);
template
void cuda_array_mask_copy<double>(const double* const array,  size_t    nRow,    size_t    nCol,    size_t    nSec,
                                  const double* const mask,   size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                        double* const output, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff);

// 1D replace
template <typename T> __global__
void dev_array_mask_replace(      T* const array, size_t    nRow,
                            const T* const mask,  size_t    nRowMsk,
                                  T        value, ptrdiff_t nRowOff,
                                  bool     maskInside)
{
    ptrdiff_t    iRowMsk = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRowMsk < nRowMsk) {
        if (mask[iRowMsk] > 0) {
            if (maskInside) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);

                array[iRowMsk+nRowOff] = value;
            }
            else {
                if (iRowMsk + nRowOff >= 0 &&
                    iRowMsk + nRowOff < (ptrdiff_t) nRow) {
                    array[iRowMsk+nRowOff] = value;
                }
            }
        }
    }
}

template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,
                             const T* const mask,  size_t    nRowMsk,
                                   T        value, ptrdiff_t nRowOff,
                                   bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);

    dev_array_mask_replace
        <<<iDivUp(nRowMsk, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, nRow,
         mask,  nRowMsk,
         value, nRowOff,
         maskInside);
}

// instantiation
template
void cuda_array_mask_replace<float >(      float*  const array, size_t    nRow,
                                     const float*  const mask,  size_t    nRowMsk,
                                           float         value, ptrdiff_t nRowOff,
                                           bool          maskInside);
template
void cuda_array_mask_replace<double>(      double* const array, size_t    nRow,
                                     const double* const mask,  size_t    nRowMsk,
                                           double        value, ptrdiff_t nRowOff,
                                           bool          maskInside);

// 2D replace
template <typename T> __global__
void dev_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,
                            const T* const mask,  size_t    nRowMsk, size_t    nColMsk,
                                  T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                  bool     maskInside)
{
    ptrdiff_t    iRowMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iColMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iRowMsk < nRowMsk && iColMsk < nColMsk) {
        iMsk = iRowMsk*nColMsk + iColMsk;

        if (mask[iMsk] > 0) {
            if (maskInside) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                assert(iColMsk + nColOff >= 0);
                assert(iColMsk + nColOff <  (ptrdiff_t) nCol);

                iSrc = (iRowMsk+nRowOff)*nCol + (iColMsk+nColOff);
                array[iSrc] = value;
            }
            else {
                if (iRowMsk + nRowOff >= 0 &&
                    iRowMsk + nRowOff < (ptrdiff_t) nRow &&
                    iColMsk + nColOff >= 0 &&
                    iColMsk + nColOff < (ptrdiff_t) nCol) {
                    iSrc = (iRowMsk+nRowOff)*nCol + (iColMsk+nColOff);
                    array[iSrc] = value;
                }
            }
        }
    }
}

template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,
                             const T* const mask,  size_t    nRowMsk, size_t    nColMsk,
                                   T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                   bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_mask_replace<<<dimGrid, dimBlock>>>
        (array, nRow,    nCol,
         mask,  nRowMsk, nColMsk,
         value, nRowOff, nColOff,
         maskInside);
}

// instantiation
template
void cuda_array_mask_replace<float >(      float*  const array, size_t    nRow,    size_t    nCol,
                                     const float*  const mask,  size_t    nRowMsk, size_t    nColMsk,
                                           float         value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                           bool          maskInside);
template
void cuda_array_mask_replace<double>(      double* const array, size_t    nRow,    size_t    nCol,
                                     const double* const mask,  size_t    nRowMsk, size_t    nColMsk,
                                           double        value, ptrdiff_t nRowOff, ptrdiff_t nColOff,
                                           bool          maskInside);

// 3D replace
template <typename T> __global__
void dev_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                            const T* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                  T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                  bool     maskInside = true)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    ptrdiff_t    iRowMsk;
    ptrdiff_t    iColMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iColMsk < nColMsk && iSecMsk < nSecMsk) {
        for (iRowMsk = 0; iRowMsk < nRowMsk; iRowMsk++) {
            iMsk = iRowMsk*nColMsk*nSecMsk + iColMsk*nSecMsk + iSecMsk;

            if (mask[iMsk] > 0) {
                if (maskInside) {
                    assert(iRowMsk + nRowOff >= 0);
                    assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                    assert(iColMsk + nColOff >= 0);
                    assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
                    assert(iSecMsk + nSecOff >= 0);
                    assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

                    iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
                    array[iSrc] = value;
                }
                else {
                    if (iRowMsk + nRowOff >= 0 &&
                        iRowMsk + nRowOff < (ptrdiff_t) nRow &&
                        iColMsk + nColOff >= 0 &&
                        iColMsk + nColOff < (ptrdiff_t) nCol &&
                        iSecMsk + nSecOff >= 0 &&
                        iSecMsk + nSecOff < (ptrdiff_t) nSec) {
                        iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
                        array[iSrc] = value;
                    }
                }
            }
        }
    }
#else
    ptrdiff_t    iRowMsk = blockDim.z * blockIdx.z + threadIdx.z;
    ptrdiff_t    iColMsk = blockDim.y * blockIdx.y + threadIdx.y;
    ptrdiff_t    iSecMsk = blockDim.x * blockIdx.x + threadIdx.x;
    ptrdiff_t    iSrc, iMsk;

    if (iRowMsk < nRowMsk && iColMsk < nColMsk && iSecMsk < nSecMsk) {
        iMsk = iRowMsk*nColMsk*nSecMsk + iColMsk*nSecMsk + iSecMsk;

        if (mask[iMsk] > 0) {
            if (maskInside) {
                assert(iRowMsk + nRowOff >= 0);
                assert(iRowMsk + nRowOff <  (ptrdiff_t) nRow);
                assert(iColMsk + nColOff >= 0);
                assert(iColMsk + nColOff <  (ptrdiff_t) nCol);
                assert(iSecMsk + nSecOff >= 0);
                assert(iSecMsk + nSecOff <  (ptrdiff_t) nSec);

                iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
                array[iSrc] = value;
            }
            else {
                if (iRowMsk + nRowOff >= 0 &&
                    iRowMsk + nRowOff < (ptrdiff_t) nRow &&
                    iColMsk + nColOff >= 0 &&
                    iColMsk + nColOff < (ptrdiff_t) nCol &&
                    iSecMsk + nSecOff >= 0 &&
                    iSecMsk + nSecOff < (ptrdiff_t) nSec) {
                    iSrc = (iRowMsk+nRowOff)*nCol*nSec + (iColMsk+nColOff)*nSec + (iSecMsk+nSecOff);
                    array[iSrc] = value;
                }
            }
        }
    }
#endif
}

template <typename T>
void cuda_array_mask_replace(      T* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                             const T* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                   T        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                   bool     maskInside)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0 && nRowMsk > 0);
    assert(nCol > 0 && nColMsk > 0);
    assert(nSec > 0 && nSecMsk > 0);

#ifdef __GEM_CUDA_ARCH_HOST_130__
    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y));
#else
    dim3    dimBlock(BLOCK_3D_NSEC, BLOCK_3D_NCOL, BLOCK_3D_NROW);
    dim3    dimGrid (iDivUp(nSec, dimBlock.x),
                     iDivUp(nCol, dimBlock.y),
                     iDivUp(nRow, dimBlock.z));
#endif

    dev_array_mask_replace<<<dimGrid, dimBlock>>>
        (array, nRow,    nCol,    nSec,
         mask,  nRowMsk, nColMsk, nSecMsk,
         value, nRowOff, nColOff, nSecOff,
         maskInside);
}

// instantiation
template
void cuda_array_mask_replace<float >(      float*  const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                                     const float*  const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                           float         value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                           bool          maskInside);
template
void cuda_array_mask_replace<double>(      double* const array, size_t    nRow,    size_t    nCol,    size_t    nSec,
                                     const double* const mask,  size_t    nRowMsk, size_t    nColMsk, size_t    nSecMsk,
                                           double        value, ptrdiff_t nRowOff, ptrdiff_t nColOff, ptrdiff_t nSecOff,
                                           bool          maskInside);

// sum
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_mask_sum(const T* const array, const T* const mask, size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        if (mask[i]) {
            myVal += array[i];
        }
        if ((i + blockSize < nRow) && mask[i+blockSize]) {
            myVal += array[i+blockSize];
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { sData[tid] = myVal = myVal + sData[tid+512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sData[tid] = myVal = myVal + sData[tid+256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sData[tid] = myVal = myVal + sData[tid+128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sData[tid] = myVal = myVal + sData[tid+ 64]; } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;

        if (blockSize >=  64) { sDataVol[tid] = myVal = myVal + sDataVol[tid+32]; }
        if (blockSize >=  32) { sDataVol[tid] = myVal = myVal + sDataVol[tid+16]; }
        if (blockSize >=  16) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 8]; }
        if (blockSize >=   8) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 4]; }
        if (blockSize >=   4) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 2]; } __syncthreads();
        if (blockSize >=   2) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 1]; }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
    }
}

template <typename T>
T cuda_array_mask_sum(const T* const array, const T* const mask, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    unsigned int    smemSize = (BLOCK_1D_NROW <= 32) ?
                               (2*BLOCK_1D_NROW * (unsigned int) sizeof(T)) :
                               (  BLOCK_1D_NROW * (unsigned int) sizeof(T));

    T         *arrayRedc = NULL;
    T         *dev_arrayRedc = NULL;

    array_new(arrayRedc, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayRedc, nblocksCUDA);

    if (isPow2(nRow)) {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_mask_sum<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_mask_sum<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_mask_sum<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_mask_sum<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_mask_sum<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_mask_sum<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_mask_sum<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_mask_sum<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_mask_sum<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_mask_sum<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_mask_sum<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_mask_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_mask_sum<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_mask_sum<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_mask_sum<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_mask_sum<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_mask_sum<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_mask_sum<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_mask_sum<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_mask_sum<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_mask_sum<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_mask_sum<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_mask_sum<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, mask, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_mask_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);

    T         sum = array_reduce_sum(arrayRedc, nblocksCUDA);

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);

    return sum;
}

// instantiation
template
float  cuda_array_mask_sum<float >(const float*  const array, const float*  const mask, size_t nRow, unsigned int nblocksCUDA);
template
double cuda_array_mask_sum<double>(const double* const array, const double* const mask, size_t nRow, unsigned int nblocksCUDA);

// substract
template <typename T> __global__
void dev_array_mask_sub(T* const arrayDst, const T* const arraySrc, const T* const mask, T value, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow && mask[iRow] > 0) {
        arrayDst[iRow] = arraySrc[iRow] - value;
    }
}

template <typename T>
void cuda_array_mask_sub(T* const arrayDst, const T* const arraySrc, const T* const mask, T value, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(mask     != NULL);
    assert(nRow > 0);

    dev_array_mask_sub
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayDst, arraySrc, mask, value, nRow);
}

// instantiation
template
void cuda_array_mask_sub<float >(float*  const arrayDst, const float*  const arraySrc, const float  *mask, float  value, size_t nRow);
template
void cuda_array_mask_sub<double>(double* const arrayDst, const double* const arraySrc, const double *mask, double value, size_t nRow);

// substract - square - sum
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_mask_subsqrsum(const T* const array, T value, const T* const mask, size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        if (mask[i]) {
            myVal += (array[i]-value)*(array[i]-value);
        }
        if ((i + blockSize < nRow) && mask[i+blockSize]) {
            myVal += (array[i+blockSize]-value)*(array[i+blockSize]-value);
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { sData[tid] = myVal = myVal + sData[tid+512]; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sData[tid] = myVal = myVal + sData[tid+256]; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sData[tid] = myVal = myVal + sData[tid+128]; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sData[tid] = myVal = myVal + sData[tid+ 64]; } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;

        if (blockSize >=  64) { sDataVol[tid] = myVal = myVal + sDataVol[tid+32]; }
        if (blockSize >=  32) { sDataVol[tid] = myVal = myVal + sDataVol[tid+16]; }
        if (blockSize >=  16) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 8]; }
        if (blockSize >=   8) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 4]; }
        if (blockSize >=   4) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 2]; }
        if (blockSize >=   2) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 1]; }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
    }
}

template <typename T>
T cuda_array_mask_subsqrsum(const T* const array, T value, const T* const mask, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    unsigned int    smemSize = (BLOCK_1D_NROW <= 32) ?
                               (2*BLOCK_1D_NROW * (unsigned int) sizeof(T)) :
                               (  BLOCK_1D_NROW * (unsigned int) sizeof(T));

    T         *arrayRedc = NULL;
    T         *dev_arrayRedc = NULL;

    array_new(arrayRedc, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayRedc, nblocksCUDA);

    if (isPow2(nRow)) {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_mask_subsqrsum<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_mask_subsqrsum<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_mask_subsqrsum<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_mask_subsqrsum<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_mask_subsqrsum<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_mask_subsqrsum<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_mask_subsqrsum<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_mask_subsqrsum<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_mask_subsqrsum<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_mask_subsqrsum<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_mask_subsqrsum<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_mask_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_mask_subsqrsum<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_mask_subsqrsum<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_mask_subsqrsum<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_mask_subsqrsum<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_mask_subsqrsum<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_mask_subsqrsum<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_mask_subsqrsum<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_mask_subsqrsum<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_mask_subsqrsum<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_mask_subsqrsum<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_mask_subsqrsum<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, value, mask, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_mask_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);

    T         sum = array_reduce_sum(arrayRedc, nblocksCUDA);

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);

    return sum;
}

// instantiation
template
float  cuda_array_mask_subsqrsum(const float*  const array, float  value, const float*  const mask, size_t nRow, unsigned int nblocksCUDA);
template
double cuda_array_mask_subsqrsum(const double* const array, double value, const double* const mask, size_t nRow, unsigned int nblocksCUDA);

// substract - divide - multiply
template <typename T> __global__
void dev_array_mask_subdivmul(T* const array, T valsub, T valdiv, const T* const mask, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        array[iRow] = (array[iRow] - valsub) / valdiv * mask[iRow];
    }
}

template <typename T>
void cuda_array_mask_subdivmul(T* const array, T valsub, T valdiv, const T* const mask, size_t nRow)
{
    assert(array != NULL);
    assert(mask  != NULL);
    assert(nRow > 0);

    dev_array_mask_subdivmul
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (array, valsub, valdiv, mask, nRow);
}

// instantiation
template
void cuda_array_mask_subdivmul<float >(float*  const array, float  valsub, float  valdiv, const float*  const mask, size_t nRow);
template
void cuda_array_mask_subdivmul<double>(double* const array, double valsub, double valdiv, const double* const mask, size_t nRow);

} // namespace gem
