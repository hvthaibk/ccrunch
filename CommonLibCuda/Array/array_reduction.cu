/***********************************************************************
 *  File:     array_reduction.cu
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

#define MAX_INDX1(myMax,myIndx,data,i) \
{ \
    if (myMax < data[i]) { \
        myMax  = data[i]; \
        myIndx = i; \
    } \
}

#define MAX_INDX2(myMax,myIndx,data,indx,i) \
{ \
    if (myMax < data[i]) { \
        myMax  = data[i]; \
        myIndx = indx[i]; \
    } \
}

#define MIN_INDX1(myMin,myIndx,data,i) \
{ \
    if (myMin > data[i]) { \
        myMin  = data[i]; \
        myIndx = i; \
    } \
}

#define MIN_INDX2(myMax,myIndx,data,indx,i) \
{ \
    if (myMax > data[i]) { \
        myMax  = data[i]; \
        myIndx = indx[i]; \
    } \
}

#ifndef max
  #define max(a,b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
  #define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

/*****************************************
 * Reduction
 ****************************************/

// max
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_max(const T* const array, size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = array[0];

    while (i < nRow) {
        myVal = max(myVal,array[i]);
        if (i + blockSize < nRow) {
            myVal = max(myVal,array[i+blockSize]);
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { sData[tid] = myVal = max(myVal,sData[tid+512]); } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sData[tid] = myVal = max(myVal,sData[tid+256]); } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sData[tid] = myVal = max(myVal,sData[tid+128]); } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sData[tid] = myVal = max(myVal,sData[tid+ 64]); } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;

        if (blockSize >=  64) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+32]); }
        if (blockSize >=  32) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+16]); }
        if (blockSize >=  16) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+ 8]); }
        if (blockSize >=   8) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+ 4]); }
        if (blockSize >=   4) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+ 2]); }
        if (blockSize >=   2) { sDataVol[tid] = myVal = max(myVal,sDataVol[tid+ 1]); }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
    }
}

template <typename T>
T cuda_array_reduce_max(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
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
            //case 1024:  dev_array_reduce_max<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_max<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_max<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_max<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_max<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_max<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_max<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_max<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_max<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_max<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_max<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_max(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_max<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_max<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_max<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_max<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_max<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_max<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_max<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_max<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_max<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_max<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_max<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_max(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);

    T         max = array_reduce_max(arrayRedc, nblocksCUDA);

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);

    return max;
}

// instantiation
template
int32_t  cuda_array_reduce_max<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_max<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_max<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_max<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_max<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_max<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// min
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_min(const T* const array, size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = array[0];

    while (i < nRow) {
        myVal = min(myVal,array[i]);
        if (i + blockSize < nRow) {
            myVal = min(myVal,array[i+blockSize]);
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { sData[tid] = myVal = min(myVal,sData[tid+512]); } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { sData[tid] = myVal = min(myVal,sData[tid+256]); } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { sData[tid] = myVal = min(myVal,sData[tid+128]); } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { sData[tid] = myVal = min(myVal,sData[tid+ 64]); } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;

        if (blockSize >=  64) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+32]); }
        if (blockSize >=  32) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+16]); }
        if (blockSize >=  16) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+ 8]); }
        if (blockSize >=   8) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+ 4]); }
        if (blockSize >=   4) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+ 2]); }
        if (blockSize >=   2) { sDataVol[tid] = myVal = min(myVal,sDataVol[tid+ 1]); }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
    }
}

template <typename T>
T cuda_array_reduce_min(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
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
            //case 1024:  dev_array_reduce_min<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_min<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_min<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_min<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_min<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_min<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_min<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_min<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_min<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_min<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_min<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_min(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_min<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_min<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_min<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_min<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_min<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_min<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_min<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_min<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_min<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_min<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_min<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_min(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);

    T         min = array_reduce_min(arrayRedc, nblocksCUDA);

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);

    return min;
}

// instantiation
template
int32_t  cuda_array_reduce_min<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_min<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_min<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_min<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_min<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_min<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// sum
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_sum(const T* const array, size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        myVal += array[i];
        if (i + blockSize < nRow) {
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
        if (blockSize >=   4) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 2]; }
        if (blockSize >=   2) { sDataVol[tid] = myVal = myVal + sDataVol[tid+ 1]; }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
    }
}

template <typename T>
T cuda_array_reduce_sum(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
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
            //case 1024:  dev_array_reduce_sum<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_sum<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_sum<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_sum<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_sum<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_sum<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_sum<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_sum<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_sum<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_sum<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_sum<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_sum<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_sum<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_sum<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_sum<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_sum<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_sum<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_sum<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_sum<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_sum<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_sum<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_sum<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_sum(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
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
int32_t  cuda_array_reduce_sum<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_sum<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_sum<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_sum<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_sum<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_sum<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// sum of squares
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_sum2(const T* const array, size_t nRow, T* const arrayRedc)
{

    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        myVal += array[i]*array[i];
        if (i + blockSize < nRow) {
            myVal += array[i+blockSize]*array[i+blockSize];
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
T cuda_array_reduce_sum2(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
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
            //case 1024:  dev_array_reduce_sum2<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_sum2<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_sum2<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_sum2<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_sum2<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_sum2<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_sum2<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_sum2<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_sum2<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_sum2<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_sum2<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_sum2(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_sum2<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_sum2<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_sum2<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_sum2<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_sum2<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_sum2<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_sum2<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_sum2<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_sum2<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_sum2<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_sum2<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_sum2(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);

    T         sum2 = array_reduce_sum(arrayRedc, nblocksCUDA);

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);

    return sum2;
}

// instantiation
template
int32_t  cuda_array_reduce_sum2<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_sum2<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_sum2<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_sum2<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_sum2<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_sum2<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// mean
template <typename T>
T cuda_array_reduce_mean(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    return cuda_array_reduce_sum(array, nRow, nblocksCUDA) / (T) nRow;
}

// instantiation
template
int32_t  cuda_array_reduce_mean<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_mean<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_mean<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_mean<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_mean<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_mean<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// std
template <typename T>
T cuda_array_reduce_std(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    T         mean, sum2;
    T         *arrayTmp = NULL;

    cuda_arrayDev_new(arrayTmp, nRow);

    mean = cuda_array_reduce_sum(array, nRow, nblocksCUDA) / (T) nRow;

    cuda_array_math_sub(arrayTmp, array, mean, nRow);

    sum2 = cuda_array_reduce_sum2(arrayTmp, nRow, nblocksCUDA);

    cuda_arrayDev_delete(arrayTmp);

    return (T) std::sqrt(sum2 / (T) nRow);
}

// instantiation
template
int32_t  cuda_array_reduce_std<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_std<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_std<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_std<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_std<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_std<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// index max
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_maxindx(const T* const array, size_t nRow, T* const arrayRedc, size_t* const arrayIndx)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    __shared__ T    sData[blockSize];
    __shared__ T    sIndx[blockSize];

    T         myVal  = array[0];
    T         myIndx = 0;

    while (i < nRow) {
        MAX_INDX1(myVal,myIndx,array,i);
        if (i + blockSize < nRow) {
            MAX_INDX1(myVal,myIndx,array,i+blockSize);
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    sIndx[tid] = myIndx;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { MAX_INDX2(myVal,myIndx,sData,sIndx,tid+512); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { MAX_INDX2(myVal,myIndx,sData,sIndx,tid+256); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { MAX_INDX2(myVal,myIndx,sData,sIndx,tid+128); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { MAX_INDX2(myVal,myIndx,sData,sIndx,tid+ 64); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;
        volatile T    *sIndxVol = sIndx;

        if (blockSize >=  64) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+32); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=  32) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+16); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=  16) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 8); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   8) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 4); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   4) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 2); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   2) { MAX_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 1); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
        arrayIndx[blockIdx.x] = myIndx;
    }
}

template <typename T>
size_t cuda_array_reduce_maxindx(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    unsigned int    smemSize = (BLOCK_1D_NROW <= 32) ?
                               (2*BLOCK_1D_NROW * (unsigned int) sizeof(T)) :
                               (  BLOCK_1D_NROW * (unsigned int) sizeof(T));

    T         *arrayRedc = NULL;
    T         *dev_arrayRedc = NULL;
    size_t    *arrayIndx = NULL;
    size_t    *dev_arrayIndx = NULL;

    array_new(arrayRedc, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayRedc, nblocksCUDA);
    array_new(arrayIndx, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayIndx, nblocksCUDA);


    if (isPow2(nRow)) {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_maxindx<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  512:  dev_array_reduce_maxindx<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  256:  dev_array_reduce_maxindx<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  128:  dev_array_reduce_maxindx<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   64:  dev_array_reduce_maxindx<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   32:  dev_array_reduce_maxindx<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   16:  dev_array_reduce_maxindx<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    8:  dev_array_reduce_maxindx<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    4:  dev_array_reduce_maxindx<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    2:  dev_array_reduce_maxindx<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    1:  dev_array_reduce_maxindx<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            default  :  std::cout << "cuda_array_reduce_maxindx(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_maxindx<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  512:  dev_array_reduce_maxindx<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  256:  dev_array_reduce_maxindx<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  128:  dev_array_reduce_maxindx<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   64:  dev_array_reduce_maxindx<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   32:  dev_array_reduce_maxindx<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   16:  dev_array_reduce_maxindx<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    8:  dev_array_reduce_maxindx<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    4:  dev_array_reduce_maxindx<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    2:  dev_array_reduce_maxindx<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    1:  dev_array_reduce_maxindx<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            default  :  std::cout << "cuda_array_reduce_maxindx(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);
    cuda_array_memcpy_d2h(arrayIndx, dev_arrayIndx, nblocksCUDA);

    T         maxVal  = arrayRedc[0];
    size_t    maxIndx = arrayIndx[0];

    for (size_t i = 1; i < nblocksCUDA; i++) {
        if (arrayRedc[i] > maxVal) {
            maxVal  = arrayRedc[i];
            maxIndx = arrayIndx[i];
        }
    }

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);
    array_delete(arrayIndx);
    cuda_arrayDev_delete(dev_arrayIndx);

    return maxIndx;
}

// instantiation
template
size_t cuda_array_reduce_maxindx<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_maxindx<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_maxindx<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_maxindx<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_maxindx<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_maxindx<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// index min
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_minindx(const T* const array, size_t nRow, T* const arrayRedc, size_t* const arrayIndx)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    __shared__ T    sData[blockSize];
    __shared__ T    sIndx[blockSize];

    T         myVal  = array[0];
    T         myIndx = 0;

    while (i < nRow) {
        MIN_INDX1(myVal,myIndx,array,i);
        if (i + blockSize < nRow) {
            MIN_INDX1(myVal,myIndx,array,i+blockSize);
        }
        i += gridSize;
    }

    sData[tid] = myVal;
    sIndx[tid] = myIndx;
    __syncthreads();

    // reduction in shared memory
    if (blockSize >= 1024) { if (tid < 512) { MIN_INDX2(myVal,myIndx,sData,sIndx,tid+512); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  512) { if (tid < 256) { MIN_INDX2(myVal,myIndx,sData,sIndx,tid+256); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  256) { if (tid < 128) { MIN_INDX2(myVal,myIndx,sData,sIndx,tid+128); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }
    if (blockSize >=  128) { if (tid <  64) { MIN_INDX2(myVal,myIndx,sData,sIndx,tid+ 64); sData[tid] = myVal; sIndx[tid] = myIndx; } __syncthreads(); }

    if (tid < 32) {
        volatile T    *sDataVol = sData;
        volatile T    *sIndxVol = sIndx;

        if (blockSize >=  64) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+32); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=  32) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+16); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=  16) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 8); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   8) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 4); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   4) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 2); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
        if (blockSize >=   2) { MIN_INDX2(myVal,myIndx,sDataVol,sIndxVol,tid+ 1); sDataVol[tid] = myVal; sIndxVol[tid] = myIndx; }
    }

    // result from this block
    if (tid == 0) {
        arrayRedc[blockIdx.x] = sData[0];
        arrayIndx[blockIdx.x] = myIndx;
    }
}

template <typename T>
size_t cuda_array_reduce_minindx(const T* const array, size_t nRow, unsigned int nblocksCUDA)
{
    assert(array != NULL);
    assert(nRow > 0);
    assert(isPow2(nblocksCUDA));

    unsigned int    smemSize = (BLOCK_1D_NROW <= 32) ?
                               (2*BLOCK_1D_NROW * (unsigned int) sizeof(T)) :
                               (  BLOCK_1D_NROW * (unsigned int) sizeof(T));

    T         *arrayRedc = NULL;
    T         *dev_arrayRedc = NULL;
    size_t    *arrayIndx = NULL;
    size_t    *dev_arrayIndx = NULL;

    array_new(arrayRedc, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayRedc, nblocksCUDA);
    array_new(arrayIndx, nblocksCUDA);
    cuda_arrayDev_new(dev_arrayIndx, nblocksCUDA);


    if (isPow2(nRow)) {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_minindx<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  512:  dev_array_reduce_minindx<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  256:  dev_array_reduce_minindx<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  128:  dev_array_reduce_minindx<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   64:  dev_array_reduce_minindx<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   32:  dev_array_reduce_minindx<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   16:  dev_array_reduce_minindx<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    8:  dev_array_reduce_minindx<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    4:  dev_array_reduce_minindx<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    2:  dev_array_reduce_minindx<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    1:  dev_array_reduce_minindx<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            default  :  std::cout << "cuda_array_reduce_minindx(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_minindx<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  512:  dev_array_reduce_minindx<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  256:  dev_array_reduce_minindx<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case  128:  dev_array_reduce_minindx<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   64:  dev_array_reduce_minindx<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   32:  dev_array_reduce_minindx<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case   16:  dev_array_reduce_minindx<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    8:  dev_array_reduce_minindx<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    4:  dev_array_reduce_minindx<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    2:  dev_array_reduce_minindx<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            case    1:  dev_array_reduce_minindx<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(array, nRow, dev_arrayRedc, dev_arrayIndx);  break;
            default  :  std::cout << "cuda_array_reduce_minindx(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }

    cuda_array_memcpy_d2h(arrayRedc, dev_arrayRedc, nblocksCUDA);
    cuda_array_memcpy_d2h(arrayIndx, dev_arrayIndx, nblocksCUDA);

    T         minVal  = arrayRedc[0];
    size_t    minIndx = arrayIndx[0];

    for (size_t i = 1; i < nblocksCUDA; i++) {
        if (arrayRedc[i] < minVal) {
            minVal  = arrayRedc[i];
            minIndx = arrayIndx[i];
        }
    }

    array_delete(arrayRedc);
    cuda_arrayDev_delete(dev_arrayRedc);
    array_delete(arrayIndx);
    cuda_arrayDev_delete(dev_arrayIndx);

    return minIndx;
}

// instantiation
template
size_t cuda_array_reduce_minindx<int32_t >(const int32_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_minindx<uint32_t>(const uint32_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_minindx<int64_t >(const int64_t*  const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_minindx<uint64_t>(const uint64_t* const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_minindx<float   >(const float*    const array, size_t nRow, unsigned int nblocksCUDA);
template
size_t cuda_array_reduce_minindx<double  >(const double*   const array, size_t nRow, unsigned int nblocksCUDA);

// corr
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_corr(const T* const arraySrc1,
                           const T* const arraySrc2,
                           size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        myVal += arraySrc1[i]*arraySrc2[i];
        if (i + blockSize < nRow) {
            myVal += arraySrc1[i+blockSize]*arraySrc2[i+blockSize];
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
T cuda_array_reduce_corr(const T* const arraySrc1,
                         const T* const arraySrc2,
                         size_t nRow, unsigned int nblocksCUDA)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
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
            //case 1024:  dev_array_reduce_corr<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_corr<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_corr<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_corr<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_corr<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_corr<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_corr<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_corr<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_corr<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_corr<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_corr<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_corr(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_corr<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_corr<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_corr<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_corr<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_corr<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_corr<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_corr<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_corr<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_corr<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_corr<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_corr<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_corr(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
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
int32_t  cuda_array_reduce_corr(const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_corr(const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_corr(const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_corr(const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_corr(const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_corr(const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow, unsigned int nblocksCUDA);

// conv
template <typename T, size_t blockSize, bool nRowIsPow2> __global__
void dev_array_reduce_conv(const T* const arraySrc1,
                           const T* const arraySrc2,
                           size_t nRow, T* const arrayRedc)
{
    size_t    tid      = threadIdx.x;
    size_t    i        = blockIdx.x*blockSize*2 + threadIdx.x;
    size_t    gridSize = blockSize*2*gridDim.x;

    T         *sData = SharedMemory<T>();
    T         myVal  = 0;

    while (i < nRow) {
        myVal += arraySrc1[i]*arraySrc2[nRow-1-i];
        if (i + blockSize < nRow) {
            myVal += arraySrc1[nRow-1-i-blockSize]*arraySrc2[nRow-1-i-blockSize];
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
T cuda_array_reduce_conv(const T* const arraySrc1,
                         const T* const arraySrc2,
                         size_t nRow, unsigned int nblocksCUDA)
{
    assert(arraySrc1 != NULL);
    assert(arraySrc2 != NULL);
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
            //case 1024:  dev_array_reduce_conv<T, 1024, true ><<<nblocksCUDA, 1024, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_conv<T,  512, true ><<<nblocksCUDA,  512, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_conv<T,  256, true ><<<nblocksCUDA,  256, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_conv<T,  128, true ><<<nblocksCUDA,  128, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_conv<T,   64, true ><<<nblocksCUDA,   64, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_conv<T,   32, true ><<<nblocksCUDA,   32, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_conv<T,   16, true ><<<nblocksCUDA,   16, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_conv<T,    8, true ><<<nblocksCUDA,    8, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_conv<T,    4, true ><<<nblocksCUDA,    4, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_conv<T,    2, true ><<<nblocksCUDA,    2, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_conv<T,    1, true ><<<nblocksCUDA,    1, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_conv(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
        }
    }
    else {
        switch (BLOCK_1D_NROW) {
            //case 1024:  dev_array_reduce_conv<T, 1024, false><<<nblocksCUDA, 1024, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  512:  dev_array_reduce_conv<T,  512, false><<<nblocksCUDA,  512, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  256:  dev_array_reduce_conv<T,  256, false><<<nblocksCUDA,  256, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case  128:  dev_array_reduce_conv<T,  128, false><<<nblocksCUDA,  128, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   64:  dev_array_reduce_conv<T,   64, false><<<nblocksCUDA,   64, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   32:  dev_array_reduce_conv<T,   32, false><<<nblocksCUDA,   32, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case   16:  dev_array_reduce_conv<T,   16, false><<<nblocksCUDA,   16, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    8:  dev_array_reduce_conv<T,    8, false><<<nblocksCUDA,    8, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    4:  dev_array_reduce_conv<T,    4, false><<<nblocksCUDA,    4, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    2:  dev_array_reduce_conv<T,    2, false><<<nblocksCUDA,    2, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            case    1:  dev_array_reduce_conv<T,    1, false><<<nblocksCUDA,    1, smemSize>>>(arraySrc1, arraySrc2, nRow, dev_arrayRedc);  break;
            default  :  std::cout << "cuda_array_reduce_conv(): unsupported value of BLOCK_1D_NROW = " << BLOCK_1D_NROW << std::endl;    exit(EXIT_FAILURE);
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
int32_t  cuda_array_reduce_conv(const int32_t*  const arraySrc1, const int32_t*  const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
uint32_t cuda_array_reduce_conv(const uint32_t* const arraySrc1, const uint32_t* const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
int64_t  cuda_array_reduce_conv(const int64_t*  const arraySrc1, const int64_t*  const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
uint64_t cuda_array_reduce_conv(const uint64_t* const arraySrc1, const uint64_t* const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
float    cuda_array_reduce_conv(const float*    const arraySrc1, const float*    const arraySrc2, size_t nRow, unsigned int nblocksCUDA);
template
double   cuda_array_reduce_conv(const double*   const arraySrc1, const double*   const arraySrc2, size_t nRow, unsigned int nblocksCUDA);

} // namespace gem
