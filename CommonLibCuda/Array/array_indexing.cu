/***********************************************************************
 *  File:       array_indexing.cu
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
 * Indexing
 ****************************************/

// 1D index -> 2D sub-indices
template <typename T> __global__
void dev_array_index_ind2sub(const T* const arrayInd,
                                   T* const arrayRow,
                                   T* const arrayCol,
                                   size_t   length,
                                   size_t   nRow, size_t nCol)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        assert(arrayInd[i] < nRow*nCol);

        arrayRow[i] = arrayInd[i] / nCol;
        arrayCol[i] = arrayInd[i] % nCol;
    }
}

template <typename T>
void cuda_array_index_ind2sub(const T* const arrayInd,
                                    T* const arrayRow,
                                    T* const arrayCol,
                                    size_t   length,
                                    size_t   nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    dev_array_index_ind2sub<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayInd,
         arrayRow, arrayCol,
         length,
         nRow, nCol);
}

// instantiation
template
void cuda_array_index_ind2sub<size_t>(const size_t* const arrayInd,
                                            size_t* const arrayRow,
                                            size_t* const arrayCol,
                                            size_t        length,
                                            size_t        nRow, size_t nCol);

// 1D index -> 3D sub-indices
template <typename T> __global__
void dev_array_index_ind2sub(const T* const arrayInd,
                                   T* const arrayRow,
                                   T* const arrayCol,
                                   T* const arraySec,
                                   size_t   length,
                                   size_t   nRow, size_t nCol, size_t nSec)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;
    size_t    nColSec = nCol*nSec;

    if (i < length) {
        assert(arrayInd[i] < nRow*nCol*nSec);

        arrayRow[i] = arrayInd[i] / nColSec;
        arrayCol[i] = (arrayInd[i] - arrayRow[i]*nColSec) / nSec;
        arraySec[i] = arrayInd[i] % nSec;
    }
}

template <typename T>
void cuda_array_index_ind2sub(const T* const arrayInd,
                                    T* const arrayRow,
                                    T* const arrayCol,
                                    T* const arraySec,
                                    size_t   length,
                                    size_t   nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arraySec != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    dev_array_index_ind2sub<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayInd,
         arrayRow, arrayCol, arraySec,
         length,
         nRow, nCol, nSec);
}

// instantiation
template
void cuda_array_index_ind2sub<size_t>(const size_t* const arrayInd,
                                            size_t* const arrayRow,
                                            size_t* const arrayCol,
                                            size_t* const arraySec,
                                            size_t        length,
                                            size_t        nRow, size_t nCol, size_t nSec);

// 2D sub-indices -> 1D index
template <typename T> __global__
void dev_array_index_sub2ind(const T* const arrayRow,
                             const T* const arrayCol,
                                   T* const arrayInd,
                                   size_t   length,
                                   size_t   nRow, size_t nCol)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        assert(arrayRow[i] < nRow);
        assert(arrayCol[i] < nCol);

        arrayInd[i] = arrayRow[i]*nCol + arrayCol[i];
    }
}

template <typename T>
void cuda_array_index_sub2ind(const T* const arrayRow,
                              const T* const arrayCol,
                                    T* const arrayInd,
                                    size_t   length,
                                    size_t   nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0);

    dev_array_index_sub2ind<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRow, arrayCol,
         arrayInd,
         length,
         nRow, nCol);
}

// instantiation
template
void cuda_array_index_sub2ind<size_t>(const size_t* const arrayRow,
                                      const size_t* const arrayCol,
                                            size_t* const arrayInd,
                                            size_t        length,
                                            size_t        nRow, size_t nCol);

// 3D sub-indices -> 1D index
template <typename T> __global__
void dev_array_index_sub2ind(const T* const arrayRow,
                             const T* const arrayCol,
                             const T* const arraySec,
                                   T* const arrayInd,
                                   size_t   length,
                                   size_t   nRow, size_t nCol, size_t nSec)
{
    size_t    i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < length) {
        assert(arrayRow[i] < nRow);
        assert(arrayCol[i] < nCol);
        assert(arraySec[i] < nSec);

        arrayInd[i] = (arrayRow[i]*nCol + arrayCol[i])*nSec + arraySec[i];
    }
}

template <typename T>
void cuda_array_index_sub2ind(const T* const arrayRow,
                              const T* const arrayCol,
                              const T* const arraySec,
                                    T* const arrayInd,
                                    size_t   length,
                                    size_t   nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(arraySec != NULL);
    assert(arrayInd != NULL);
    assert(length > 0);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    dev_array_index_sub2ind<<<iDivUp(length, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (arrayRow, arrayCol, arraySec,
         arrayInd,
         length,
         nRow, nCol, nSec);
}

// instantiation
template
void cuda_array_index_sub2ind<size_t>(const size_t* const arrayRow,
                                      const size_t* const arrayCol,
                                      const size_t* const arraySec,
                                            size_t* const arrayInd,
                                            size_t        length,
                                            size_t        nRow, size_t nCol, size_t nSec);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2> __global__
void dev_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                             size_t nRow, size_t nCol)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow && iCol < nCol) {
        arrayRow[iRow*nCol+iCol] = (T2) arrayCol[iCol*nRow+iRow];
    }
}

template <typename T1, typename T2>
void cuda_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                              size_t nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_index_col2row<<<dimGrid, dimBlock>>>
        (arrayCol, arrayRow, nRow, nCol);
}

// instantiation
template
void cuda_array_index_col2row<int32_t ,int32_t >(const int32_t*  const arrayCol, int32_t*  const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<uint32_t,uint32_t>(const uint32_t* const arrayCol, uint32_t* const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<int64_t ,int64_t >(const int64_t*  const arrayCol, int64_t*  const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<uint64_t,uint64_t>(const uint64_t* const arrayCol, uint64_t* const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<float   ,float   >(const float*    const arrayCol, float*    const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<float   ,double  >(const float*    const arrayCol, double*   const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<double  ,float   >(const double*   const arrayCol, float*    const arrayRow,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_col2row<double  ,double  >(const double*   const arrayCol, double*   const arrayRow,
                                                 size_t nRow, size_t nCol);

// column-major -> row-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2> __global__
void dev_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                             size_t nRow, size_t nCol, size_t nSec)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRow;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    if (iCol < nCol && iSec < nSec) {
        for (iRow = 0; iRow < nRow; iRow++) {
            arrayRow[(iRow*nCol+iCol)*nSec+iSec] =
            (T2) arrayCol[(iSec*nCol+iCol)*nRow+iRow];
        }
    }
#else
    size_t    iRow = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow && iCol < nCol && iSec < nSec) {
        arrayRow[(iRow*nCol+iCol)*nSec+iSec] =
        (T2) arrayCol[(iSec*nCol+iCol)*nRow+iRow];
    }
#endif
}

template <typename T1, typename T2>
void cuda_array_index_col2row(const T1* const arrayCol, T2* const arrayRow,
                              size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

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

    dev_array_index_col2row<<<dimGrid, dimBlock>>>
        (arrayCol, arrayRow, nRow, nCol, nSec);
}

// instantiation
template
void cuda_array_index_col2row<int32_t ,int32_t >(const int32_t*  const arrayCol, int32_t*  const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<uint32_t,uint32_t>(const uint32_t* const arrayCol, uint32_t* const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<int64_t ,int64_t >(const int64_t*  const arrayCol, int64_t*  const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<uint64_t,uint64_t>(const uint64_t* const arrayCol, uint64_t* const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<float   ,float   >(const float*    const arrayCol, float*    const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<float   ,double  >(const float*    const arrayCol, double*   const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<double  ,float   >(const double*   const arrayCol, float*    const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_col2row<double  ,double  >(const double*   const arrayCol, double*   const arrayRow,
                                                 size_t nRow, size_t nCol, size_t nSec);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow)
template <typename T1, typename T2> __global__
void dev_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                             size_t nRow, size_t nCol)
{
    size_t    iRow = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iCol = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow && iCol < nCol) {
        arrayCol[iCol*nRow+iRow] = (T2) arrayRow[iRow*nCol+iCol];
    }
}

template <typename T1, typename T2>
void cuda_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                              size_t nRow, size_t nCol)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0);

    dim3    dimBlock(BLOCK_2D_NCOL, BLOCK_2D_NROW);
    dim3    dimGrid (iDivUp(nCol, dimBlock.x),
                     iDivUp(nRow, dimBlock.y));

    dev_array_index_row2col<<<dimGrid, dimBlock>>>
        (arrayRow, arrayCol, nRow, nCol);
}

// instantiation
template
void cuda_array_index_row2col<int32_t ,int32_t >(const int32_t*  const arrayRow, int32_t*  const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<uint32_t,uint32_t>(const uint32_t* const arrayRow, uint32_t* const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<int64_t ,int64_t >(const int64_t*  const arrayRow, int64_t*  const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<uint64_t,uint64_t>(const uint64_t* const arrayRow, uint64_t* const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<float   ,float   >(const float*    const arrayRow, float*    const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<float   ,double  >(const float*    const arrayRow, double*   const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<double  ,float   >(const double*   const arrayRow, float*    const arrayCol,
                                                 size_t nRow, size_t nCol);
template
void cuda_array_index_row2col<double  ,double  >(const double*   const arrayRow, double*   const arrayCol,
                                                 size_t nRow, size_t nCol);

// row-major -> column-major
// nRow = dim1(arrayRow), nCol = dim2(arrayRow), nSec = dim3(arrayRow)
template <typename T1, typename T2> __global__
void dev_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                             size_t nRow, size_t nCol, size_t nSec)
{
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
    size_t    iRow;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    if (iCol < nCol && iSec < nSec) {
        for (iRow = 0; iRow < nRow; iRow++) {
            arrayCol[(iSec*nCol+iCol)*nRow+iRow] =
            (T2) arrayRow[(iRow*nCol+iCol)*nSec+iSec];
        }
    }
#else
    size_t    iRow = blockDim.z * blockIdx.z + threadIdx.z;
    size_t    iCol = blockDim.y * blockIdx.y + threadIdx.y;
    size_t    iSec = blockDim.x * blockIdx.x + threadIdx.x;

    if (iRow < nRow && iCol < nCol && iSec < nSec) {
        arrayCol[(iSec*nCol+iCol)*nRow+iRow] =
        (T2) arrayRow[(iRow*nCol+iCol)*nSec+iSec];
    }
#endif
}

template <typename T1, typename T2>
void cuda_array_index_row2col(const T1* const arrayRow, T2* const arrayCol,
                              size_t nRow, size_t nCol, size_t nSec)
{
    assert(arrayRow != NULL);
    assert(arrayCol != NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

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

    dev_array_index_row2col<<<dimGrid, dimBlock>>>
        (arrayRow, arrayCol, nRow, nCol, nSec);
}

// instantiation
template
void cuda_array_index_row2col<int32_t ,int32_t >(const int32_t*  const arrayRow, int32_t*  const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<uint32_t,uint32_t>(const uint32_t* const arrayRow, uint32_t* const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<int64_t ,int64_t >(const int64_t*  const arrayRow, int64_t*  const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<uint64_t,uint64_t>(const uint64_t* const arrayRow, uint64_t* const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<float   ,float   >(const float*    const arrayRow, float*    const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<float   ,double  >(const float*    const arrayRow, double*   const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<double  ,float   >(const double*   const arrayRow, float*    const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);
template
void cuda_array_index_row2col<double  ,double  >(const double*   const arrayRow, double*   const arrayCol,
                                                 size_t nRow, size_t nCol, size_t nSec);

} // namespace gem
