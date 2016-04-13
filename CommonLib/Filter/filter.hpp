/***********************************************************************
 *  File:       filter.hpp
 *
 *  Purpose:    Header file for filtering functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_FILTER_HPP__
#define __GEM_LIB_FILTER_HPP__

#include "array.hpp"

namespace gem {

/**************************
 * Average
 *************************/

// 1D
template <typename T>
void filter_average(T* const array, size_t nRow);

// 2D
template <typename T>
void filter_average(T* const array, size_t nRow, size_t nCol);

// 3D
template <typename T>
void filter_average(T* const array, size_t nRow, size_t nCol, size_t nSec);

/**************************
 * Disk
 *************************/

// 1D
template <typename T>
void filter_disk(T* const array, size_t nRow, T radius, bool norm = 1);

// 2D
template <typename T>
void filter_disk(T* const array, size_t nRow, size_t nCol, T radius, bool norm = 1);

// 3D
template <typename T>
void filter_disk(T* const array, size_t nRow, size_t nCol, size_t nSec, T radius, bool norm = 1);

/**************************
 * Gaussian
 *************************/

// 1D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, T sigma = 0.5);

// 2D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, size_t nCol, T sigma = 0.5);

// 3D
template <typename T>
void filter_gaussian(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma = 0.5);

/**************************
 * Difference of Gaussian
 *************************/

// 1D
template <typename T>
void filter_dog(T* const array, size_t nRow, T sigma1, T sigma2);

// 2D
template <typename T>
void filter_dog(T* const array, size_t nRow, size_t nCol, T sigma1, T sigma2);

// 3D
template <typename T>
void filter_dog(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma1, T sigma2);

/**************************
 * Laplacian of Gaussian
 *************************/

// 1D
template <typename T>
void filter_log(T* const array, size_t nRow, T sigma = 0.5);

// 2D
template <typename T>
void filter_log(T* const array, size_t nRow, size_t nCol, T sigma = 0.5);

// 3D
template <typename T>
void filter_log(T* const array, size_t nRow, size_t nCol, size_t nSec, T sigma = 0.5);

/**************************
 * Laplacian
 *************************/

// 1D
template <typename T>
void filter_laplacian(T* const array, size_t nRow);

// 2D
template <typename T>
void filter_laplacian(T* const array, size_t nRow, size_t nCol, T alpha = 0.2);

// 3D
template <typename T>
void filter_laplacian(T* const array, size_t nRow, size_t nCol, size_t nSec, T alpha = 3.0/16, T beta = 0);

} // namespace gem

#endif
