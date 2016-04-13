/***********************************************************************
 *  File:       cufft_wrapper.cu
 *
 *  Purpose:    Implementation of cufft-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cufft_wrapper.cuh"

namespace gem {

/*****************************************
 * Forward and inverse FFTs (slow)
 ****************************************/

// 2D float: forward
void cufft_dft_r2c(float          *dataReal, size_t nRowReal, size_t nColReal,
                   cuFloatComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan2d(&fftPlan, (int) nRowReal, (int) nColReal, CUFFT_R2C) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecR2C(fftPlan, (cufftReal *)    dataReal,
                                           (cufftComplex *) dataCplx) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 2D float: inverse
void cufft_dft_c2r(cuFloatComplex *dataCplx,
                   float          *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan2d(&fftPlan, (int) nRowReal, (int) nColReal, CUFFT_C2R) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecC2R(fftPlan, (cufftComplex *) dataCplx,
                                           (cufftReal *)    dataReal) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 2D double: forward
void cufft_dft_r2c(double          *dataReal, size_t nRowReal, size_t nColReal,
                   cuDoubleComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan2d(&fftPlan, (int) nRowReal, (int) nColReal, CUFFT_D2Z) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecD2Z(fftPlan, (cufftDoubleReal *)    dataReal,
                                           (cufftDoubleComplex *) dataCplx) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 2D double: inverse
void cufft_dft_c2r(cuDoubleComplex *dataCplx,
                   double          *dataReal, size_t nRowReal, size_t nColReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan2d(&fftPlan, (int) nRowReal, (int) nColReal, CUFFT_D2Z) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecZ2D(fftPlan, (cufftDoubleComplex *) dataCplx,
                                           (cufftDoubleReal *)    dataReal) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 3D float: forward
void cufft_dft_r2c(float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                   cuFloatComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan3d(&fftPlan, (int) nRowReal, (int) nColReal, (int) nSecReal, CUFFT_R2C) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecR2C(fftPlan, (cufftReal *)    dataReal,
                                           (cufftComplex *) dataCplx) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 3D float: inverse
void cufft_dft_c2r(cuFloatComplex *dataCplx,
                   float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan3d(&fftPlan, (int) nRowReal, (int) nColReal, (int) nSecReal, CUFFT_C2R) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecC2R(fftPlan, (cufftComplex *) dataCplx,
                                           (cufftReal *)    dataReal) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 3D double: forward
void cufft_dft_r2c(double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                   cuDoubleComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan3d(&fftPlan, (int) nRowReal, (int) nColReal, (int) nSecReal, CUFFT_D2Z) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecD2Z(fftPlan, (cufftDoubleReal *)    dataReal,
                                           (cufftDoubleComplex *) dataCplx) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

// 3D double: inverse
void cufft_dft_c2r(cuDoubleComplex *dataCplx,
                   double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0 && nSecReal > 0);

    // setup CUFFT plan
    cufftHandle     fftPlan;
    CUFFT_SAFE_CALL( cufftPlan3d(&fftPlan, (int) nRowReal, (int) nColReal, (int) nSecReal, CUFFT_D2Z) );

    // obtain the FFT
    CUFFT_SAFE_CALL( cufftExecZ2D(fftPlan, (cufftDoubleComplex *) dataCplx,
                                           (cufftDoubleReal *)    dataReal) );

    // deallocate CUFFT plan
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan) );
}

/*****************************************
 * Forward and inverse FFTs (1d1d1d)
 ****************************************/

// 2D float: forward
void cufft_dft_r2c_fast(float          *dataReal, size_t nRowReal, size_t nColReal,
                        cuFloatComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate CUFFT input and output arrays
    cuFloatComplex    *dataTmp1 = NULL;
    cuda_arrayDev_new(dataTmp1, fftSize);

    // setup CUFFT plans
    cufftHandle fftPlan1, fftPlan2;
    CUFFT_SAFE_CALL( cufftPlan1d(&fftPlan1, (int) nColReal, CUFFT_R2C, (int) nRowReal) );
    CUFFT_SAFE_CALL( cufftPlan1d(&fftPlan2, (int) nRowReal, CUFFT_C2C, (int) nColCplx) );

    // batch FFT 1
    CUFFT_SAFE_CALL( cufftExecR2C(fftPlan1, (cufftReal *)    dataReal,
                                            (cufftComplex *) dataTmp1) );

    // permute 1
    cuda_array_permute(reinterpret_cast<cuFloatComplex*>(dataTmp1),
                       reinterpret_cast<cuFloatComplex*>(dataCplx),
                       nRowReal, nColCplx);

    // batch FFT 2
    CUFFT_SAFE_CALL( cufftExecC2C(fftPlan2, (cufftComplex *) dataCplx,
                                            (cufftComplex *) dataTmp1,
                                            CUFFT_FORWARD) );

    // permute 2
    cuda_array_permute(reinterpret_cast<cuFloatComplex*>(dataTmp1),
                       reinterpret_cast<cuFloatComplex*>(dataCplx),
                       nColCplx, nRowReal);

    // deallocate CUFFT arrays and plans
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan1) );
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan2) );
    cuda_arrayDev_delete(dataTmp1);
}

// 2D float: forward
/*void cufft_dft_r2c_fast(float          *dataReal, size_t nRowReal, size_t nColReal,
                        cuFloatComplex *dataCplx)
{
    assert(dataReal != NULL);
    assert(dataCplx != NULL);
    assert(nRowReal > 0 && nColReal > 0);

    size_t    nColCplx = nColReal / 2 + 1;
    size_t    fftSize  = nRowReal * nColCplx;

    // allocate CUFFT input and output arrays
    cuFloatComplex    *dataTmp1 = NULL;
    cuda_arrayDev_new(dataTmp1, fftSize);

    // setup CUFFT plans
    int         n1[] = {(int) nColReal};
    int         n2[] = {(int) nRowReal};
    cufftHandle fftPlan1, fftPlan2;
    CUFFT_SAFE_CALL( cufftPlanMany(&fftPlan1, (int) 1,  n1,
                                   n1, (int) 1, (int) nColReal,
                                   n1, (int) 1, (int) nColCplx,
                                   CUFFT_R2C, (int) nRowReal) );
    CUFFT_SAFE_CALL( cufftPlanMany(&fftPlan2, (int) 1,  n2,
                                   n2, (int) nColCplx, (int) 1,
                                   n2, (int) nColCplx, (int) 1,
                                   CUFFT_C2C, (int) nColCplx) );
    //CUFFT_SAFE_CALL( cufftPlanMany(&fftPlan2, (int) 1,  n2,
    //                               n2, (int) 1, (int) nRowReal,
    //                               n2, (int) 1, (int) nRowReal,
    //                               CUFFT_C2C, (int) nColCplx) );

    // batch FFT 1
    CUFFT_SAFE_CALL( cufftExecR2C(fftPlan1, (cufftReal *)    dataReal,
                                            (cufftComplex *) dataTmp1) );

    // batch FFT 2
    CUFFT_SAFE_CALL( cufftExecC2C(fftPlan2, (cufftComplex *) dataCplx,
                                            (cufftComplex *) dataTmp1,
                                            CUFFT_FORWARD) );

    // deallocate CUFFT arrays and plans
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan1) );
    CUFFT_SAFE_CALL( cufftDestroy(fftPlan2) );
    cuda_arrayDev_delete(dataTmp1);
}*/

// 2D float: inverse
void cufft_dft_c2r_fast(cuFloatComplex *dataCplx,
                        float          *dataReal, size_t nRowReal, size_t nColReal);

// 2D double: forward
void cufft_dft_r2c_fast(double          *dataReal, size_t nRowReal, size_t nColReal,
                        cuDoubleComplex *dataCplx);

// 2D double: inverse
void cufft_dft_c2r_fast(cuDoubleComplex *dataCplx,
                        double          *dataReal, size_t nRowReal, size_t nColReal);

// 3D float: forward
void cufft_dft_r2c_fast(float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuFloatComplex *dataCplx);

// 3D float: inverse
void cufft_dft_c2r_fast(cuFloatComplex *dataCplx,
                        float          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal);

// 3D double: forward
void cufft_dft_r2c_fast(double          *dataReal, size_t nRowReal, size_t nColReal, size_t nSecReal,
                        cuDoubleComplex *dataCplx);

// 3D double: inverse
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
