/***********************************************************************
 *  File:       transform_expand_rotation.cu
 *
 *  Purpose:    Implementation of expansion functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.cuh"

namespace gem {

/**************************
 * Expansion
 *************************/

// 2D
template <typename T>
void cuda_transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff;

    diagLength = (size_t) std::ceil(dist((T)0,(T)0,(T)nRowSrc,(T)nColSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;

    cuda_arrayDev_delete(arrayDst);
    cuda_arrayDev_new_zero(arrayDst, nRowDst*nColDst);

    cuda_array_replace<T,T>(arrayDst, nRowDst, nColDst,
                            arraySrc, nRowSrc, nColSrc,
                                      nRowOff, nColOff);
}

// instantiation
template
void cuda_transform_expand_rotate_full<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                     float*&       arrayDst, size_t& nRowDst, size_t& nColDst);
template
void cuda_transform_expand_rotate_full<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc,
                                                     double*&      arrayDst, size_t& nRowDst, size_t& nColDst);

// 3D
template <typename T>
void cuda_transform_expand_rotate_full(const T* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                             T*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst)
{
    assert(arraySrc != NULL && arrayDst == NULL);
    assert(nRowSrc > 0 && nColSrc > 0 && nSecSrc > 0);

    size_t    diagLength;
    size_t    nRowOff, nColOff, nSecOff;

    diagLength = (size_t) std::ceil(dist((T)0,(T)0,(T)0,(T)nRowSrc,(T)nColSrc,(T)nSecSrc));

    nRowDst = ((diagLength-nRowSrc) % 2) ? diagLength+1 : diagLength;
    nColDst = ((diagLength-nColSrc) % 2) ? diagLength+1 : diagLength;
    nSecDst = ((diagLength-nSecSrc) % 2) ? diagLength+1 : diagLength;

    nRowOff = (nRowDst-nRowSrc) / 2;
    nColOff = (nColDst-nColSrc) / 2;
    nSecOff = (nSecDst-nSecSrc) / 2;

    cuda_arrayDev_delete(arrayDst);
    cuda_arrayDev_new_zero(arrayDst, nRowDst*nColDst*nSecDst);

    cuda_array_replace<T,T>(arrayDst, nRowDst, nColDst, nSecDst,
                            arraySrc, nRowSrc, nColSrc, nSecSrc,
                                      nRowOff, nColOff, nSecOff);
}

// instantiation
template
void cuda_transform_expand_rotate_full<float >(const float*  const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                     float*&       arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);
template
void cuda_transform_expand_rotate_full<double>(const double* const arraySrc, size_t  nRowSrc, size_t  nColSrc, size_t  nSecSrc,
                                                     double*&      arrayDst, size_t& nRowDst, size_t& nColDst, size_t& nSecDst);

} // namespace gem
