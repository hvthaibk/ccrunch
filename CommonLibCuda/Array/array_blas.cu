/***********************************************************************
 *  File:       array_blas.cu
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.cuh"

namespace gem {

/*****************************************
 * BLAS
 ****************************************/

// blas axpy
template <typename T> __global__
void dev_array_blas_axpy(      T* const arrayY,
                         const T* const arrayX,
                         T alpha, size_t nRow)
{
    size_t    iRow = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow) {
        arrayY[iRow] += alpha*arrayX[iRow];
    }
}

template <typename T>
void cuda_array_blas_axpy(      T* const arrayY,
                          const T* const arrayX,
                          T alpha, size_t nRow)
{
    assert(arrayX != NULL);
    assert(arrayY != NULL);
    assert(nRow > 0);

    dev_array_blas_axpy
        <<<iDivUp(nRow, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayY, arrayX, alpha, nRow);
}

// instantiation
template
void cuda_array_blas_axpy<float >(      float*  const arrayY,
                                  const float*  const arrayX,
                                  float  alpha, size_t nRow);
template
void cuda_array_blas_axpy<double>(      double* const arrayY,
                                  const double* const arrayX,
                                  double alpha, size_t nRow);

// blas gemm: C = A*B, nColA = nRowB
template <typename T, unsigned int BLOCK_SIZE> __global__
void dev_array_blas_gemm2(      T* const arrayC,
                          const T* const arrayA,
                          const T* const arrayB,
                          size_t nRowA, size_t nColA, size_t nColB)
{
    size_t    iBRow = blockIdx.y;
    size_t    iBCol = blockIdx.x;
    size_t    iTRow = threadIdx.y;
    size_t    iTCol = threadIdx.x;
    size_t    iRow  = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol  = blockDim.x * blockIdx.x + threadIdx.x;

    __shared__ T    As[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ T    Bs[BLOCK_SIZE][BLOCK_SIZE];
               T    Ctmp = 0;
    size_t          a = iBRow * BLOCK_SIZE * nColA ;
    size_t          b = iBCol * BLOCK_SIZE;

    for (size_t i = 0; i < (nColA-1)/BLOCK_SIZE+1; i++) {

        As[iTRow][iTCol] = arrayA[a + nColA * iTRow + iTCol];
        Bs[iTRow][iTCol] = arrayB[b + nColB * iTRow + iTCol];
        a +=         BLOCK_SIZE;
        b += nColB * BLOCK_SIZE;

        __syncthreads();

        for (size_t k = 0; k < BLOCK_SIZE; ++k) {
            Ctmp += As[iTRow][k] * Bs[k][iTCol];
        }

        __syncthreads();
    }

    arrayC[iRow*nColB+iCol] = Ctmp;
}

template <typename T, unsigned int BLOCK_SIZE> __global__
void dev_array_blas_gemm1(      T* const arrayC,
                          const T* const arrayA,
                          const T* const arrayB,
                          size_t nRowA, size_t nColA, size_t nColB)
{
    size_t    iBRow = blockIdx.y;
    size_t    iBCol = blockIdx.x;
    size_t    iTRow = threadIdx.y;
    size_t    iTCol = threadIdx.x;
    size_t    iRow  = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol  = blockDim.x * blockIdx.x + threadIdx.x;

    __shared__ T    As[BLOCK_SIZE][BLOCK_SIZE];
    __shared__ T    Bs[BLOCK_SIZE][BLOCK_SIZE];
               T    Ctmp = 0;
    size_t          a = iBRow * BLOCK_SIZE * nColA ;
    size_t          b = iBCol * BLOCK_SIZE;

    for (size_t i = 0; i < (nColA-1)/BLOCK_SIZE+1; i++) {

        if (iBRow*BLOCK_SIZE+iTRow < nRowA && i*BLOCK_SIZE+iTCol < nColA) {
               As[iTRow][iTCol] = arrayA[a + nColA * iTRow + iTCol];      }
        else { As[iTRow][iTCol] = 0; }
        a += BLOCK_SIZE;

        if (i*BLOCK_SIZE+iTRow < nColA && iBCol*BLOCK_SIZE+iTCol < nColB) {
               Bs[iTRow][iTCol] = arrayB[b + nColB * iTRow + iTCol];      }
        else { Bs[iTRow][iTCol] = 0; }
        b += nColB * BLOCK_SIZE;

        __syncthreads();

        for (size_t k = 0; k < BLOCK_SIZE; ++k) {
            Ctmp += As[iTRow][k] * Bs[k][iTCol];
        }

        __syncthreads();
    }

    if (iRow < nRowA && iCol < nColB) {
        arrayC[iRow*nColB+iCol] = Ctmp;
    }
}

template <typename T> __global__
void dev_array_blas_gemm0(      T* const arrayC,
                         const T* const arrayA,
                         const T* const arrayB,
                         size_t nRowA, size_t nColA, size_t nColB)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    T         Ctmp = 0;

    if (iRow < nRowA && iCol < nColB) {
        for (size_t i = 0; i < nColA; i++) {
            Ctmp += arrayA[iRow*nColA+i] * arrayB[i*nColB+iCol];
        }

        arrayC[iRow*nColB+iCol] = Ctmp;
    }
}

template <typename T>
void cuda_array_blas_gemm(      T* const arrayC,
                          const T* const arrayA,
                          const T* const arrayB,
                          size_t nRowA, size_t nColA, size_t nColB,
                          unsigned int mode)
{
    assert(arrayA != NULL);
    assert(arrayB != NULL);
    assert(arrayC != NULL);
    assert(nRowA > 0);
    assert(nColA > 0);
    assert(nColB > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nColB, dimBlock.x),
                     iDivUp(nRowA, dimBlock.y));

    switch (mode) {
        case 0:
            dev_array_blas_gemm0<<<dimGrid, dimBlock>>>
                (arrayC, arrayA, arrayB, nRowA, nColA, nColB);
            break;
        case 1:
            require(BLOCK_2D_NROW == BLOCK_2D_NCOL,
                    "cuda_array_blas_gemm: BLOCK_2D_NROW != BLOCK_2D_NCOL");

            dev_array_blas_gemm1<T,BLOCK_2D_NROW><<<dimGrid, dimBlock>>>
                (arrayC, arrayA, arrayB, nRowA, nColA, nColB);
            break;
        case 2:
            require(BLOCK_2D_NROW == BLOCK_2D_NCOL,
                    "cuda_array_blas_gemm: BLOCK_2D_NROW != BLOCK_2D_NCOL");
            require(nRowA % BLOCK_2D_NROW == 0 &&
                    nColA % BLOCK_2D_NROW == 0 &&
                    nColB % BLOCK_2D_NROW == 0,
                    "cuda_array_blas_gemm: conditions on matrix sizes are not satisfied");

            dev_array_blas_gemm2<T,BLOCK_2D_NROW><<<dimGrid, dimBlock>>>
                (arrayC, arrayA, arrayB, nRowA, nColA, nColB);
            break;
        default:
            ERROR("cuda_array_blas_gemm", "unsupported multiplication mode")
    }
}

// instantiation
template
void cuda_array_blas_gemm<float >(      float*  const arrayC,
                                  const float*  const arrayA,
                                  const float*  const arrayB,
                                  size_t nRowA, size_t nColA, size_t nColB,
                                  unsigned int mode);
template
void cuda_array_blas_gemm<double>(      double* const arrayC,
                                  const double* const arrayA,
                                  const double* const arrayB,
                                  size_t nRowA, size_t nColA, size_t nColB,
                                  unsigned int mode);

} // namespace gem
