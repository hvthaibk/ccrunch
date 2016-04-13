/***********************************************************************
 *  File:       fftw_wrapper.hpp
 *
 *  Purpose:    Header file for fftw-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_FFTW_WRAPPER_HPP__
#define __GEM_LIB_FFTW_WRAPPER_HPP__

#include <fftw3.h>

typedef     fftwf_complex   fftwComplexSingle;
typedef     fftw_complex    fftwComplexDouble;

typedef     fftwf_plan      fftwPlanSingle;
typedef     fftw_plan       fftwPlanDouble;

namespace gem {

/*****************************************
 * malloc / memcheck
 ****************************************/

// malloc
void fftw_new(fftwComplexSingle*& array1D, size_t nRow);
void fftw_new(fftwComplexDouble*& array1D, size_t nRow);

// free
void fftw_delete(fftwComplexSingle*& array1D);
void fftw_delete(fftwComplexDouble*& array1D);

// check memory leak
void fftw_memcheck(void);

/*****************************************
 * fft plans
 ****************************************/

void fftw_create_plan(fftwPlanSingle* fftwPlan);
void fftw_create_plan(fftwPlanDouble* fftwPlan);

/*****************************************
 * Forward and inverse FFTs (slow)
 ****************************************/

// 1D float
void fftw_dft_r2c(float*             dataReal, size_t nRowReal,
                  fftwComplexSingle* dataCplx,
                  fftwPlanSingle*    fftwPlanFwd);

void fftw_dft_c2r(fftwComplexSingle* dataCplx,
                  float*             dataReal, size_t nRowReal);

// 1D double
void fftw_dft_r2c(double*            dataReal, size_t nRowReal,
                  fftwComplexDouble* dataCplx);

void fftw_dft_c2r(fftwComplexDouble* dataCplx,
                  double*            dataReal, size_t nRowReal);

// 2D float
void fftw_dft_r2c(float         *dataReal, size_t nRowReal, size_t nColReal,
                  fftwf_complex *dataCplx);

void fftw_dft_c2r(fftwf_complex *dataCplx,
                  float         *dataReal, size_t nRowReal, size_t nColReal);

// 2D double
void fftw_dft_r2c(double       *dataReal, size_t nRowReal, size_t nColReal,
                  fftw_complex *dataCplx);

void fftw_dft_c2r(fftw_complex *dataCplx,
                  double       *dataReal, size_t nRowReal, size_t nColReal);

// 3D float
void fftw_dft_r2c(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                  fftwf_complex *dataCplx);

void fftw_dft_c2r(fftwf_complex *dataCplx,
                  float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void fftw_dft_r2c(double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                  fftw_complex *dataCplx);

void fftw_dft_c2r(fftw_complex *dataCplx,
                  double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

/*****************************************
 * Forward and inverse FFTs (1d1d1d)
 ****************************************/

// 2D float
void fftw_dft_r2c_fast(float         *dataReal, size_t nRowReal, size_t nColReal,
                       fftwf_complex *dataCplx);

void fftw_dft_c2r_fast(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal);

// 2D double
void fftw_dft_r2c_fast(double       *dataReal, size_t nRowReal, size_t nColReal,
                       fftw_complex *dataCplx);

void fftw_dft_c2r_fast(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal);

// 3D float
void fftw_dft_r2c_fast(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftwf_complex *dataCplx);

void fftw_dft_c2r_fast(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void fftw_dft_r2c_fast(double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftw_complex *dataCplx);

void fftw_dft_c2r_fast(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

/*****************************************
 * Forward and inverse FFTs (2d1d)
 ****************************************/

// 3D float
void fftw_dft_r2c_2d1d(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftwf_complex *dataCplx);

void fftw_dft_c2r_2d1d(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double
void fftw_dft_r2c_2d1d(double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftw_complex *dataCplx);

void fftw_dft_c2r_2d1d(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

} // namespace gem

#endif
