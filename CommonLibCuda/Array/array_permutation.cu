/***********************************************************************
 *  File:       array_permutation.cu
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
 * Permutation
 ****************************************/

// 2D
template <typename T> __global__
void dev_array_permute(const T* const arraySrc, T* const arrayDst,
                       size_t nRow, size_t nCol)
{
    __shared__ T    tile[TILE_DIM][TILE_DIM+1];

    size_t    iRow = TILE_DIM * blockIdx.y + threadIdx.y;
    size_t    iCol = TILE_DIM * blockIdx.x + threadIdx.x;

    if (iRow < nRow && iCol < nCol) {
        size_t    iIn  = iRow * nCol + iCol;

        for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
            tile[threadIdx.y+i][threadIdx.x] = arraySrc[iIn+i*nCol];
        }
    }

    __syncthreads();

    // arrayDst: #rows = nCol, #cols = nRow
    iRow = TILE_DIM * blockIdx.x + threadIdx.y;     // iRow of arrayDst
    iCol = TILE_DIM * blockIdx.y + threadIdx.x;     // iCol of arrayDst

    if (iRow < nCol && iCol < nRow) {
        size_t    iOut = iRow * nRow + iCol;

        for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
            arrayDst[iOut+i*nRow] = tile[threadIdx.x][threadIdx.y+i];
        }
    }
}

template <typename T>
void cuda_array_permute(const T* const arraySrc, T* const arrayDst,
                        size_t nRow, size_t nCol)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0);

    dim3    dimBlock(TILE_DIM, BLOCK_ROWS);
    dim3    dimGrid (iDivUp(nCol, TILE_DIM),
                     iDivUp(nRow, TILE_DIM));

    dev_array_permute<<<dimGrid, dimBlock>>>
        (arraySrc, arrayDst, nRow, nCol);
}

// instantiation
template
void cuda_array_permute<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<float   >(const float*  const arraySrc, float*  const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<double  >(const double* const arraySrc, double* const arrayDst,
                                  size_t nRow, size_t nCol);
template
void cuda_array_permute<cuFloatComplex >(const cuFloatComplex*  const arraySrc, cuFloatComplex*  const arrayDst,
                                         size_t nRow, size_t nCol);
template
void cuda_array_permute<cuDoubleComplex>(const cuDoubleComplex* const arraySrc, cuDoubleComplex* const arrayDst,
                                         size_t nRow, size_t nCol);

// 3D
template <typename T> __global__
void dev_array_permute_231(const T* const arraySrc, T* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec)
{
    __shared__ T    tile[TILE_DIM][TILE_DIM+1];

    size_t    iRow = TILE_DIM * blockIdx.y + threadIdx.y;
    size_t    iSec = TILE_DIM * blockIdx.x + threadIdx.x;

    size_t    iRowDst;
    size_t    iColDst = TILE_DIM * blockIdx.x + threadIdx.y;
    size_t    iSecDst = TILE_DIM * blockIdx.y + threadIdx.x;
    size_t    iSrc, iDst;

    for (size_t iCol = 0; iCol < nCol; iCol++) {
        if (iRow < nRow && iSec < nSec) {
            iSrc  = iRow*nCol*nSec + iCol*nSec + iSec;

            for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
                tile[threadIdx.y+i][threadIdx.x] = arraySrc[iSrc+i*nCol*nSec];
            }
        }

        __syncthreads();

        // arrayDst: #rows = nCol, #cols = nSec, #secs = nRow
        iRowDst = iCol;

        if (iColDst < nSec && iSecDst < nRow) {
            iDst = iRowDst*nSec*nRow + iColDst*nRow + iSecDst;

            for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
                arrayDst[iDst+i*nRow] = tile[threadIdx.x][threadIdx.y+i];
            }
        }

        __syncthreads();
    }
}

template <typename T> __global__
void dev_array_permute_312(const T* const arraySrc, T* const arrayDst,
                           size_t nRow, size_t nCol, size_t nSec)
{
    __shared__ T    tile[TILE_DIM][TILE_DIM+1];

    size_t    iCol = TILE_DIM * blockIdx.y + threadIdx.y;
    size_t    iSec = TILE_DIM * blockIdx.x + threadIdx.x;

    size_t    iRowDst = TILE_DIM * blockIdx.x + threadIdx.y;
    size_t    iColDst;
    size_t    iSecDst = TILE_DIM * blockIdx.y + threadIdx.x;
    size_t    iSrc, iDst;

    for (size_t iRow = 0; iRow < nRow; iRow++) {
        if (iCol < nCol && iSec < nSec) {
            iSrc  = iRow*nCol*nSec + iCol*nSec + iSec;

            for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
                tile[threadIdx.y+i][threadIdx.x] = arraySrc[iSrc+i*nSec];
            }
        }

        __syncthreads();

        // arrayDst: #rows = nSec, #cols = nRow, #secs = nCol
        iColDst = iRow;

        if (iRowDst < nSec && iSecDst < nCol) {
            iDst = iRowDst*nRow*nCol + iColDst*nCol + iSecDst;

            for (size_t i = 0; i < TILE_DIM; i += BLOCK_ROWS) {
                arrayDst[iDst+i*nRow*nCol] = tile[threadIdx.x][threadIdx.y+i];
            }
        }

        __syncthreads();
    }
}

template <typename T>
void cuda_array_permute(const T* const arraySrc, T* const arrayDst,
                        size_t nRow, size_t nCol, size_t nSec,
                        ePermute permute)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    dim3    dimBlock(TILE_DIM, BLOCK_ROWS), dimGrid;

    switch (permute) {
        case PERMUTE3D_123:
            ERROR("cuda_array_permute", "no need to permute, memcpy is faster");
            break;
        case PERMUTE3D_132:
        case PERMUTE3D_213:
        case PERMUTE3D_321:
            ERROR("cuda_array_permute", "unsupported permutation mode");
            break;
        case PERMUTE3D_231:
            dimGrid.x = iDivUp(nSec, TILE_DIM);
            dimGrid.y = iDivUp(nRow, TILE_DIM);
            dimGrid.z = 1;
            dev_array_permute_231<<<dimGrid, dimBlock>>>
                (arraySrc, arrayDst, nRow, nCol, nSec);
            break;
        case PERMUTE3D_312:
            dimGrid.x = iDivUp(nSec, TILE_DIM);
            dimGrid.y = iDivUp(nCol, TILE_DIM);
            dimGrid.z = 1;
            dev_array_permute_312<<<dimGrid, dimBlock>>>
                (arraySrc, arrayDst, nRow, nCol, nSec);
            break;
        default:
            ERROR("cuda_array_permute", "unsupported permutation mode");
    }
}

// instantiation
template
void cuda_array_permute<int32_t >(const int32_t*  const arraySrc, int32_t*  const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<uint32_t>(const uint32_t* const arraySrc, uint32_t* const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<int64_t >(const int64_t*  const arraySrc, int64_t*  const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<uint64_t>(const uint64_t* const arraySrc, uint64_t* const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<float   >(const float*  const arraySrc, float*  const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<double  >(const double* const arraySrc, double* const arrayDst,
                                  size_t nRow, size_t nCol, size_t nSec,
                                  ePermute permute);
template
void cuda_array_permute<cuFloatComplex >(const cuFloatComplex*  const arraySrc, cuFloatComplex*  const arrayDst,
                                         size_t nRow, size_t nCol, size_t nSec,
                                         ePermute permute);
template
void cuda_array_permute<cuDoubleComplex>(const cuDoubleComplex* const arraySrc, cuDoubleComplex* const arrayDst,
                                         size_t nRow, size_t nCol, size_t nSec,
                                         ePermute permute);

} // namespace gem
