/***********************************************************************
 *  File:       cufft_wrapper.cuh
 *
 *  Purpose:    Header file for cufft-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_CUFFT_WRAPPER_HPP__
#define __GEM_LIB_CUFFT_WRAPPER_HPP__

#include "xcorr.hpp"
#include "array.cuh"

#include <cufft.h>
#include <cuComplex.h>

namespace gem {

/*****************************************
 * Forward and inverse FFTs (slow)
 ****************************************/

// 2D float
void cufft_dft_r2c(float          *dataReal, size_t nRowReal, size_t nColReal,
                   cuFloatComplex *dataCplx);

void cufft_dft_c2r(cuFloatComplex *dataCplx,
                   float          *dataReal, size_t nRowReal, size_t nColReal);

// 2D double
void cufft_dft_r2c(double          *dataReal, size_t nRowReal, size_t nColReal,
                   cuDoubleComplex *dataCplx);

void cufft_dft_c2r(cuDoubleComplex *dataCplx,
                   double          *dataReal, size_t nRowReal, size_t nColReal);

// 3D float
void cufft_dft_r2c(float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                   cuFloatComplex *dataCplx);

void cufft_dft_c2r(cuFloatComplex *dataCplx,
                   float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void cufft_dft_r2c(double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                   cuDoubleComplex *dataCplx);

void cufft_dft_c2r(cuDoubleComplex *dataCplx,
                   double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

/*****************************************
 * Forward and inverse FFTs (1d1d1d)
 ****************************************/

// 2D float
void cufft_dft_r2c_fast(float          *dataReal, size_t nRowReal, size_t nColReal,
                        cuFloatComplex *dataCplx);

void cufft_dft_c2r_fast(cuFloatComplex *dataCplx,
                        float          *dataReal, size_t nRowReal, size_t nColReal);

// 2D double
void cufft_dft_r2c_fast(double          *dataReal, size_t nRowReal, size_t nColReal,
                        cuDoubleComplex *dataCplx);

void cufft_dft_c2r_fast(cuDoubleComplex *dataCplx,
                        double          *dataReal, size_t nRowReal, size_t nColReal);

// 3D float
void cufft_dft_r2c_fast(float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuFloatComplex *dataCplx);

void cufft_dft_c2r_fast(cuFloatComplex *dataCplx,
                        float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void cufft_dft_r2c_fast(double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuDoubleComplex *dataCplx);

void cufft_dft_c2r_fast(cuDoubleComplex *dataCplx,
                        double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

/*****************************************
 * Forward and inverse FFTs (2d1d)
 ****************************************/

// 3D float
void cufft_dft_r2c_2d1d(float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuFloatComplex *dataCplx);

void cufft_dft_c2r_2d1d(cuFloatComplex *dataCplx,
                        float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void cufft_dft_r2c_2d1d(double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuDoubleComplex *dataCplx);

void cufft_dft_c2r_2d1d(cuDoubleComplex *dataCplx,
                        double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

} // namespace gem

#endif
