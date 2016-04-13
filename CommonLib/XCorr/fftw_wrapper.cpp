/***********************************************************************
 *  File:       fftw_wrapper.cpp
 *
 *  Purpose:    Implementation of fftw-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.hpp"
#include "fftw_lock.hpp"
#include "fftw_wrapper.hpp"

namespace gem {

#ifndef NDEBUG
static int      countAllocFFTW = 0;
#endif

/***************************************************
 * Mutex for FFTW
 **************************************************/

inline void mutexFFTW_lock(void)   { mutexFFTW.lock();   }
inline void mutexFFTW_unlock(void) { mutexFFTW.unlock(); }

/*****************************************
 * malloc / memcheck
 ****************************************/

// malloc float
void fftw_new(fftwComplexSingle*& array1D, size_t nRow)
{
    assert(array1D == NULL);
    assert(nRow > 0);

    FFTW_LOCK( array1D = (fftwComplexSingle*) fftwf_malloc(sizeof(fftwComplexSingle) * nRow) );

#ifndef NDEBUG
    countAllocFFTW++;
#endif
}

// malloc double
void fftw_new(fftwComplexDouble*& array1D, size_t nRow)
{
    assert(array1D == NULL);
    assert(nRow > 0);

    FFTW_LOCK( array1D = (fftwComplexDouble*) fftw_malloc(sizeof(fftwComplexDouble) * nRow) );

#ifndef NDEBUG
    countAllocFFTW++;
#endif
}

// free float
void fftw_delete(fftwComplexSingle*& array1D)
{
    if (array1D != NULL) {
        FFTW_LOCK( fftwf_free(array1D) );
        array1D = NULL;

#ifndef NDEBUG
        countAllocFFTW--;
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("fftw_delete", "attempt to delete a NULL object");
#endif
    }
}

// free double
void fftw_delete(fftwComplexDouble*& array1D)
{
    if (array1D != NULL) {
        FFTW_LOCK( fftw_free(array1D) );
        array1D = NULL;

#ifndef NDEBUG
        countAllocFFTW--;
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("fftw_delete", "attempt to delete a NULL object");
#endif
    }
}

// check memory leak
void fftw_memcheck(void)
{
#ifndef NDEBUG
    assert(countAllocFFTW == 0);
#else
    WARNING("fftw_memcheck", "only available in DEBUG mode");
#endif
}

/*****************************************
 * Forward and inverse FFTs (slow)
 ****************************************/

// 1D float
void fftw_dft_r2c(float*             dataReal, size_t nRowReal,
                  fftwComplexSingle* dataCplx);

void fftw_dft_c2r(fftwComplexSingle* dataCplx,
                  float*             dataReal, size_t nRowReal);

// 1D double
void fftw_dft_r2c(double*            dataReal, size_t nRowReal,
                  fftwComplexDouble* dataCplx);

void fftw_dft_c2r(fftwComplexDouble* dataCplx,
                  double*            dataReal, size_t nRowReal);

// 2D float: forward
void fftw_dft_r2c(float         *dataReal, size_t nRowReal, size_t nColReal,
                  fftwf_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftwf_plan fftPlan;
    FFTW_LOCK( fftPlan = fftwf_plan_dft_r2c_2d((int) nRowReal, (int) nColReal,
                                               dataReal, dataCplx,
                                               FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftwf_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftwf_destroy_plan(fftPlan) );
}

// 2D float: inverse
void fftw_dft_c2r(fftwf_complex *dataCplx,
                  float         *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftwf_plan fftPlan;
    FFTW_LOCK( fftPlan = fftwf_plan_dft_c2r_2d((int) nRowReal, (int) nColReal,
                                               dataCplx, dataReal,
                                               FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftwf_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftwf_destroy_plan(fftPlan) );

    // normalize the output
    array_math_div(dataReal, (float) (nRowReal*nColReal), nRowReal*nColReal);
}

// 2D double: forward
void fftw_dft_r2c(double       *dataReal, size_t nRowReal, size_t nColReal,
                  fftw_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftw_plan fftPlan;
    FFTW_LOCK( fftPlan = fftw_plan_dft_r2c_2d((int) nRowReal, (int) nColReal,
                                              dataReal, dataCplx,
                                              FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftw_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftw_destroy_plan(fftPlan) );
}

// 2D double: inverse
void fftw_dft_c2r(fftw_complex *dataCplx,
                  double       *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftw_plan fftPlan;
    FFTW_LOCK( fftPlan = fftw_plan_dft_c2r_2d((int) nRowReal, (int) nColReal,
                                              dataCplx, dataReal,
                                              FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftw_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftw_destroy_plan(fftPlan) );

    // normalize the output
    array_math_div(dataReal, (double) (nRowReal*nColReal), nRowReal*nColReal);
}

// 3D float: forward
void fftw_dft_r2c(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                  fftwf_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    // setup FFTW plan
    fftwf_plan fftPlan;
    FFTW_LOCK( fftPlan = fftwf_plan_dft_r2c_3d((int) nRowReal, (int) nColReal, (int) nSecReal,
                                               dataReal, dataCplx,
                                               FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftwf_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftwf_destroy_plan(fftPlan) );
}

// 3D float: inverse
void fftw_dft_c2r(fftwf_complex *dataCplx,
                  float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftwf_plan fftPlan;
    FFTW_LOCK( fftPlan = fftwf_plan_dft_c2r_3d((int) nRowReal, (int) nColReal, (int) nSecReal,
                                               dataCplx, dataReal,
                                               FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftwf_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftwf_destroy_plan(fftPlan) );

    // normalize the output
    array_math_div(dataReal, (float) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);
}

// 3D double: forward
void fftw_dft_r2c(double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                  fftw_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftw_plan fftPlan;
    FFTW_LOCK( fftPlan = fftw_plan_dft_r2c_3d((int) nRowReal, (int) nColReal, (int) nSecReal,
                                              dataReal, dataCplx,
                                              FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftw_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftw_destroy_plan(fftPlan) );
}

// 3D double: inverse
void fftw_dft_c2r(fftw_complex *dataCplx,
                  double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup FFTW plan
    fftw_plan fftPlan;
    FFTW_LOCK( fftPlan = fftw_plan_dft_c2r_3d((int) nRowReal, (int) nColReal, (int) nSecReal,
                                              dataCplx, dataReal,
                                              FFTW_ESTIMATE) );
    assert(fftPlan != NULL);

    // obtain the FFT
    fftw_execute(fftPlan);

    // deallocate FFTW plan
    FFTW_LOCK( fftw_destroy_plan(fftPlan) );

    // normalize the output
    array_math_div(dataReal, (double) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);
}

/*****************************************
 * Forward and inverse FFTs (1d1d1d)
 ****************************************/

// 2D float: forward
void fftw_dft_r2c_fast(float         *dataReal, size_t nRowReal, size_t nColReal,
                       fftwf_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal};
    int         n2[] = {(int) nRowReal};
    fftwf_plan  fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft_r2c((int) 1,  n1, (int) nRowReal,
                                                  dataReal, n1, (int) 1, (int) nColReal,
                                                  dataTmp1, n1, (int) 1, (int) nColCplx,
                                                  FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft((int) 1,  n2, (int) nColCplx,
                                              dataCplx, n2, (int) 1, (int) nRowReal,
                                              dataTmp1, n2, (int) 1, (int) nRowReal,
                                              FFTW_FORWARD, FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nRowReal, nColCplx);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nColCplx, nRowReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 2D float: inverse
void fftw_dft_c2r_fast(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nRowReal};
    int         n2[] = {(int) nColReal};
    fftwf_plan  fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft((int) 1,  n1, (int) nColCplx,
                                              dataTmp1, n1, (int) 1, (int) nRowReal,
                                              dataCplx, n1, (int) 1, (int) nRowReal,
                                              FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft_c2r((int) 1,  n2, (int) nRowReal,
                                                  dataTmp1, n2, (int) 1, (int) nColCplx,
                                                  dataReal, n2, (int) 1, (int) nColReal,
                                                  FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nRowReal, nColCplx);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nColCplx, nRowReal);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // normalize the output
    array_math_div(dataReal, (float) (nRowReal*nColReal), nRowReal*nColReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 2D double: forward
void fftw_dft_r2c_fast(double       *dataReal, size_t nRowReal, size_t nColReal,
                       fftw_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate FFTW input and output arrays
    fftw_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal};
    int         n2[] = {(int) nRowReal};
    fftw_plan   fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft_r2c((int) 1,  n1, (int) nRowReal,
                                                 dataReal, n1, (int) 1, (int) nColReal,
                                                 dataTmp1, n1, (int) 1, (int) nColCplx,
                                                 FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft((int) 1,  n2, (int) nColCplx,
                                             dataCplx, n2, (int) 1, (int) nRowReal,
                                             dataTmp1, n2, (int) 1, (int) nRowReal,
                                             FFTW_FORWARD, FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nRowReal, nColCplx);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nColCplx, nRowReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 2D double: inverse
void fftw_dft_c2r_fast(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate FFTW input and output arrays
    fftw_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nRowReal};
    int         n2[] = {(int) nColReal};
    fftw_plan   fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft((int) 1,  n1, (int) nColCplx,
                                             dataTmp1, n1, (int) 1, (int) nRowReal,
                                             dataCplx, n1, (int) 1, (int) nRowReal,
                                             FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft_c2r((int) 1,  n2, (int) nRowReal,
                                                 dataTmp1, n2, (int) 1, (int) nColCplx,
                                                 dataReal, n2, (int) 1, (int) nColReal,
                                                 FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nRowReal, nColCplx);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nColCplx, nRowReal);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // normalize the output
    array_math_div(dataReal, (double) (nRowReal*nColReal), nRowReal*nColReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 3D float: forward
void fftw_dft_r2c_fast(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftwf_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nSecReal};
    int         n2[] = {(int) nColReal};
    int         n3[] = {(int) nRowReal};
    fftwf_plan  fftPlan1, fftPlan2, fftPlan3;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft_r2c((int) 1,  n1, (int) (nRowReal*nColReal),
                                                  dataReal, n1, (int) 1, (int) nSecReal,
                                                  dataTmp1, n1, (int) 1, (int) nSecCplx,
                                                  FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft((int) 1,  n2, (int) (nSecCplx*nRowReal),
                                              dataCplx, n2, (int) 1, (int) nColReal,
                                              dataTmp1, n2, (int) 1, (int) nColReal,
                                              FFTW_FORWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan3 = fftwf_plan_many_dft((int) 1,  n3, (int) (nColReal*nSecCplx),
                                              dataCplx, n3, (int) 1, (int) nRowReal,
                                              dataTmp1, n3, (int) 1, (int) nRowReal,
                                              FFTW_FORWARD, FFTW_ESTIMATE) );

    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);
    assert(fftPlan3 != NULL);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_312);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nSecCplx, nRowReal, nColReal,
                  PERMUTE3D_312);

    // batch FFT 3
    fftwf_execute(fftPlan3);

    // permute 3
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan3) );
    fftw_delete(dataTmp1);
}

// 3D float: inverse
void fftw_dft_c2r_fast(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal};
    int         n2[] = {(int) nRowReal};
    int         n3[] = {(int) nSecReal};
    fftwf_plan  fftPlan1, fftPlan2, fftPlan3;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft((int) 1,  n1, (int) (nSecCplx*nRowReal),
                                              dataTmp1, n1, (int) 1, (int) nColReal,
                                              dataCplx, n1, (int) 1, (int) nColReal,
                                              FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft((int) 1,  n2, (int) (nColReal*nSecCplx),
                                              dataTmp1, n2, (int) 1, (int) nRowReal,
                                              dataCplx, n2, (int) 1, (int) nRowReal,
                                              FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan3 = fftwf_plan_many_dft_c2r((int) 1,  n3, (int) (nRowReal*nColReal),
                                                  dataTmp1, n3, (int) 1, (int) nSecCplx,
                                                  dataReal, n3, (int) 1, (int) nSecReal,
                                                  FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);
    assert(fftPlan3 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_312);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nSecCplx, nRowReal, nColReal,
                  PERMUTE3D_312);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // permute 3
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // batch FFT 3
    fftwf_execute(fftPlan3);

    // normalize the output
    array_math_div(dataReal, (float) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan3) );
    fftw_delete(dataTmp1);
}

// 3D double: forward
void fftw_dft_r2c_fast(double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftw_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftw_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nSecReal};
    int         n2[] = {(int) nColReal};
    int         n3[] = {(int) nRowReal};
    fftw_plan   fftPlan1, fftPlan2, fftPlan3;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft_r2c((int) 1,  n1, (int) (nRowReal*nColReal),
                                                 dataReal, n1, (int) 1, (int) nSecReal,
                                                 dataTmp1, n1, (int) 1, (int) nSecCplx,
                                                 FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft((int) 1,  n2, (int) (nSecCplx*nRowReal),
                                             dataCplx, n2, (int) 1, (int) nColReal,
                                             dataTmp1, n2, (int) 1, (int) nColReal,
                                             FFTW_FORWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan3 = fftw_plan_many_dft((int) 1,  n3, (int) (nColReal*nSecCplx),
                                             dataCplx, n3, (int) 1, (int) nRowReal,
                                             dataTmp1, n3, (int) 1, (int) nRowReal,
                                             FFTW_FORWARD, FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);
    assert(fftPlan3 != NULL);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_312);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nSecCplx, nRowReal, nColReal,
                  PERMUTE3D_312);

    // batch FFT 3
    fftw_execute(fftPlan3);

    // permute 3
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan3) );
    fftw_delete(dataTmp1);
}

// 3D double: inverse
void fftw_dft_c2r_fast(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftw_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal};
    int         n2[] = {(int) nRowReal};
    int         n3[] = {(int) nSecReal};
    fftw_plan   fftPlan1, fftPlan2, fftPlan3;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft((int) 1,  n1, (int) (nSecCplx*nRowReal),
                                             dataTmp1, n1, (int) 1, (int) nColReal,
                                             dataCplx, n1, (int) 1, (int) nColReal,
                                             FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft((int) 1,  n2, (int) (nColReal*nSecCplx),
                                             dataTmp1, n2, (int) 1, (int) nRowReal,
                                             dataCplx, n2, (int) 1, (int) nRowReal,
                                             FFTW_BACKWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan3 = fftw_plan_many_dft_c2r((int) 1,  n3, (int) (nRowReal*nColReal),
                                                 dataTmp1, n3, (int) 1, (int) nSecCplx,
                                                 dataReal, n3, (int) 1, (int) nSecReal,
                                                 FFTW_ESTIMATE) );
    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);
    assert(fftPlan3 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_312);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nSecCplx, nRowReal, nColReal,
                  PERMUTE3D_312);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // permute 3
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // batch FFT 3
    fftw_execute(fftPlan3);

    // normalize the output
    array_math_div(dataReal, (double) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan3) );
    fftw_delete(dataTmp1);
}

/*****************************************
 * Forward and inverse FFTs (2d1d)
 ****************************************/

// 3D float: forward
void fftw_dft_r2c_2d1d(float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftwf_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal, (int) nSecReal};
    int         n2[] = {(int) nRowReal};
    fftwf_plan  fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft_r2c((int) 2,  n1,   (int) nRowReal,
                                                  dataReal, NULL, (int) 1, (int) (nColReal*nSecReal),
                                                  dataTmp1, NULL, (int) 1, (int) (nColReal*nSecCplx),
                                                  FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft((int) 1,  n2,   (int) (nColReal*nSecCplx),
                                              dataCplx, NULL, (int) 1, (int) nRowReal,
                                              dataTmp1, NULL, (int) 1, (int) nRowReal,
                                              FFTW_FORWARD, FFTW_ESTIMATE) );

    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_231);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataTmp1),
                  reinterpret_cast<std::complex<float>*>(dataCplx),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 3D float: inverse
void fftw_dft_c2r_2d1d(fftwf_complex *dataCplx,
                       float         *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftwf_complex   *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nRowReal};
    int         n2[] = {(int) nColReal, (int) nSecReal};
    fftwf_plan  fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftwf_plan_many_dft((int) 1,  n1,   (int) (nColReal*nSecCplx),
                                              dataTmp1, NULL, (int) 1, (int) nRowReal,
                                              dataCplx, NULL, (int) 1, (int) nRowReal,
                                              FFTW_FORWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftwf_plan_many_dft_c2r((int) 2,  n2,   (int) nRowReal,
                                                  dataTmp1, NULL, (int) 1, (int) (nColReal*nSecCplx),
                                                  dataReal, NULL, (int) 1, (int) (nColReal*nSecReal),
                                                  FFTW_ESTIMATE) );

    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_231);

    // batch FFT 1
    fftwf_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<float>*>(dataCplx),
                  reinterpret_cast<std::complex<float>*>(dataTmp1),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // batch FFT 2
    fftwf_execute(fftPlan2);

    // normalize the output
    array_math_div(dataReal, (float) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftwf_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftwf_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 3D double: forward
void fftw_dft_r2c_2d1d(double        *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                       fftw_complex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftw_complex    *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nColReal, (int) nSecReal};
    int         n2[] = {(int) nRowReal};
    fftw_plan   fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft_r2c((int) 2,  n1,   (int) nRowReal,
                                                 dataReal, NULL, (int) 1, (int) (nColReal*nSecReal),
                                                 dataTmp1, NULL, (int) 1, (int) (nColReal*nSecCplx),
                                                 FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft((int) 1,  n2,   (int) (nColReal*nSecCplx),
                                             dataCplx, NULL, (int) 1, (int) nRowReal,
                                             dataTmp1, NULL, (int) 1, (int) nRowReal,
                                             FFTW_FORWARD, FFTW_ESTIMATE) );

    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_231);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataTmp1),
                  reinterpret_cast<std::complex<double>*>(dataCplx),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

// 3D double: inverse
void fftw_dft_c2r_2d1d(fftw_complex *dataCplx,
                       double       *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    size_t    nSecCplx = nSecReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColReal * nSecCplx;

    // allocate FFTW input and output arrays
    fftw_complex    *dataTmp1 = NULL;
    fftw_new(dataTmp1, fftSize);

    // setup FFTW plans
    int         n1[] = {(int) nRowReal};
    int         n2[] = {(int) nColReal, (int) nSecReal};
    fftw_plan   fftPlan1, fftPlan2;
    FFTW_LOCK( fftPlan1 = fftw_plan_many_dft((int) 1,  n1,   (int) (nColReal*nSecCplx),
                                             dataTmp1, NULL, (int) 1, (int) nRowReal,
                                             dataCplx, NULL, (int) 1, (int) nRowReal,
                                             FFTW_FORWARD, FFTW_ESTIMATE) );
    FFTW_LOCK( fftPlan2 = fftw_plan_many_dft_c2r((int) 2,  n2,   (int) nRowReal,
                                                 dataTmp1, NULL, (int) 1, (int) (nColReal*nSecCplx),
                                                 dataReal, NULL, (int) 1, (int) (nColReal*nSecReal),
                                                 FFTW_ESTIMATE) );

    assert(fftPlan1 != NULL);
    assert(fftPlan2 != NULL);

    // permute 1
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nRowReal, nColReal, nSecCplx,
                  PERMUTE3D_231);

    // batch FFT 1
    fftw_execute(fftPlan1);

    // permute 2
    array_permute(reinterpret_cast<std::complex<double>*>(dataCplx),
                  reinterpret_cast<std::complex<double>*>(dataTmp1),
                  nColReal, nSecCplx, nRowReal,
                  PERMUTE3D_312);

    // batch FFT 2
    fftw_execute(fftPlan2);

    // normalize the output
    array_math_div(dataReal, (double) (nRowReal*nColReal*nSecReal), nRowReal*nColReal*nSecReal);

    // deallocate FFTW arrays and plans
    FFTW_LOCK( fftw_destroy_plan(fftPlan1) );
    FFTW_LOCK( fftw_destroy_plan(fftPlan2) );
    fftw_delete(dataTmp1);
}

} // namespace gem
