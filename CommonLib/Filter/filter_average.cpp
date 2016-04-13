/***********************************************************************
 *  File:       filter_average.cpp
 *
 *  Purpose:    Implementation of filtering functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "filter.hpp"

namespace gem {

/**************************
 * Average
 *************************/

// 1D
template <typename T>
void filter_average(T* const array, size_t nRow)
{
    array_setval(array, (T) 1.0/(T) (nRow), nRow);
}

// instantiation
template
void filter_average<float >(float*  const array, size_t nRow);
template
void filter_average<double>(double* const array, size_t nRow);

// 2D
template <typename T>
void filter_average(T* const array, size_t nRow, size_t nCol)
{
    array_setval(array, (T) 1.0/(T) (nRow*nCol), nRow*nCol);
}

// instantiation
template
void filter_average<float >(float*  const array, size_t nRow, size_t nCol);
template
void filter_average<double>(double* const array, size_t nRow, size_t nCol);

// 3D
template <typename T>
void filter_average(T* const array, size_t nRow, size_t nCol, size_t nSec)
{
    array_setval(array, (T) 1.0/(T) (nRow*nCol*nSec), nRow*nCol*nSec);
}

// instantiation
template
void filter_average<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec);
template
void filter_average<double>(double* const array, size_t nRow, size_t nCol, size_t nSec);

} // namespace gem
