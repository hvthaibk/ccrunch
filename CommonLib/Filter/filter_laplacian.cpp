/***********************************************************************
 *  File:       filter_laplacian.cpp
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
 * Laplacian
 *************************/

// 1D
template <typename T>
void filter_laplacian(T* const array, size_t nRow)
{
    assert(array != NULL);

    if (nRow != 3) {
        std::cout << "filter_laplacian(): "
                  << "the filter size should be 3"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    array[0] = 1;    array[1] = -2;    array[2] = 1;
}

// instantiation
template
void filter_laplacian<float >(float*  const array, size_t nRow);
template
void filter_laplacian<double>(double* const array, size_t nRow);

// 2D
template <typename T>
void filter_laplacian(T* const array, size_t nRow, size_t nCol, T alpha)
{
    assert(array != NULL);

    if ((nRow != 3) || (nCol != 3)) {
        std::cout << "filter_laplacian(): "
                  << "the filter size should be 3x3"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((alpha < 0.0) || (alpha > 1.0)) {
        std::cout << "filter_laplacian(): "
                  << "alpha shoule be in the range [0,1]"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    array[0] = alpha;    array[1] = 1-alpha;    array[2] = alpha;
    array[3] = 1-alpha;  array[4] = -4;         array[5] = 1-alpha;
    array[6] = alpha;    array[7] = 1-alpha;    array[8] = alpha;

    array_math_div(array, 1+alpha, nRow*nCol);
}

// instantiation
template
void filter_laplacian<float >(float*  const array, size_t nRow, size_t nCol, float  alpha);
template
void filter_laplacian<double>(double* const array, size_t nRow, size_t nCol, double alpha);

// 3D
template <typename T>
void filter_laplacian(T* const array, size_t nRow, size_t nCol, size_t nSec, T alpha, T beta)
{
    assert(array != NULL);

    T    gamma = 1 - 2*alpha - 8*beta/6;

    if ((nRow != 3) || (nCol != 3) || (nSec != 3)) {
        std::cout << "filter_laplacian(): "
                  << "the filter size should be 3x3x3"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((alpha < 0.0) || (alpha > 0.5)) {
        std::cout << "filter_laplacian(): "
                  << "alpha shoule be in the range [0,0.5]"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((beta < 0.0) || (beta > 6/8)) {
        std::cout << "filter_laplacian(): "
                  << "beta shoule be in the range [0,6/8]"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    if (gamma < 0) {
        std::cout << "filter_laplacian(): "
                  << "beta shoule be in the range [0,6/8]"
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    array[0]  = beta;    array[1]  = alpha;    array[2]  = beta;
    array[3]  = alpha;   array[4]  = gamma;    array[5]  = alpha;
    array[6]  = beta;    array[7]  = alpha;    array[8]  = beta;

    array[9]  = alpha;   array[10] = gamma;    array[11] = alpha;
    array[12] = gamma;   array[13] = -6;       array[14] = gamma;
    array[15] = alpha;   array[16] = gamma;    array[17] = alpha;

    array[18] = beta;    array[19] = alpha;    array[20] = beta;
    array[21] = alpha;   array[22] = gamma;    array[23] = alpha;
    array[24] = beta;    array[25] = alpha;    array[26] = beta;
}

// instantiation
template
void filter_laplacian<float >(float*  const array, size_t nRow, size_t nCol, size_t nSec, float  alpha, float  beta);
template
void filter_laplacian<double>(double* const array, size_t nRow, size_t nCol, size_t nSec, double alpha, double beta);

} // namespace gem
