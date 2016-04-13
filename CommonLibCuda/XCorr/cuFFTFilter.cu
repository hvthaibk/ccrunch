/***********************************************************************
 *  File:       cuCUFFTFilter.cu
 *
 *  Purpose:    Implementation of an CUFFT-filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuFFTFilter.cuh"

namespace gem {

template <class T1, class T2> __global__
void dev_filterModulateAndNormalize(      T1* const dataFFT,
                                    const T2* const kernel,
                                    size_t fftSize,
                                    T2     nElements)
{
    const size_t    tid = blockDim.x * blockIdx.x + threadIdx.x;

    if (tid < fftSize) {
        dataFFT[tid].x = dataFFT[tid].x * kernel[tid] / nElements;
        dataFFT[tid].y = dataFFT[tid].y * kernel[tid] / nElements;
    }
}

cuFFTFilter::cuFFTFilter(void)
{
    _kernelSingle = NULL;
    _kernelDouble = NULL;

    _dataSingle = NULL;
    _dataDouble = NULL;

    _dataSingleFFT = NULL;
    _dataDoubleFFT = NULL;
}

void cuFFTFilter::prepare(const float*  const kernel, size_t nrow, size_t ncol)
{
    _nrow    = nrow;
    _ncol    = ncol;
    _fftSize = _nrow * (_ncol / 2 + 1);
    _resSize = _nrow * _ncol;

    // memory
    cuda_arrayDev_new(_kernelSingle, _fftSize);
    cuda_array_memcpy_d2d(_kernelSingle, kernel, _fftSize);

    cuda_arrayDev_new(_dataSingle, _resSize);
    cuda_arrayDev_new(_dataSingleFFT, _fftSize);

    // plans
    CUFFT_SAFE_CALL(cufftPlan2d(&_fftPlanFwd, (int) nrow, (int) ncol, CUFFT_R2C));
    CUFFT_SAFE_CALL(cufftPlan2d(&_fftPlanInv, (int) nrow, (int) ncol, CUFFT_C2R));
}

void cuFFTFilter::prepare(const double* const kernel, size_t nrow, size_t ncol)
{
    _nrow    = nrow;
    _ncol    = ncol;
    _fftSize = _nrow * (_ncol / 2 + 1);
    _resSize = _nrow * _ncol;

    // memory
    cuda_arrayDev_new(_kernelDouble, _fftSize);
    cuda_array_memcpy_d2d(_kernelDouble, kernel, _fftSize);

    cuda_arrayDev_new(_dataDouble, _resSize);
    cuda_arrayDev_new(_dataDoubleFFT, _fftSize);

    // plans
    CUFFT_SAFE_CALL(cufftPlan2d(&_fftPlanFwd, (int) nrow, (int) ncol, CUFFT_D2Z));
    CUFFT_SAFE_CALL(cufftPlan2d(&_fftPlanInv, (int) nrow, (int) ncol, CUFFT_Z2D));
}

void cuFFTFilter::perform(float*  const data)
{
    cuda_array_memcpy_d2d(_dataSingle, data, _resSize);

    cufftExecR2C(_fftPlanFwd, (cufftReal *)    _dataSingle,
                                    (cufftComplex *) _dataSingleFFT);

    dev_filterModulateAndNormalize
        <<<iDivUp(_fftSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (_dataSingleFFT, _kernelSingle, _fftSize, (float) _resSize);

    cufftExecC2R(_fftPlanInv, (cufftComplex *) _dataSingleFFT,
                                    (cufftReal *)    _dataSingle);

    cuda_array_memcpy_d2d(data, _dataSingle, _resSize);
}

void cuFFTFilter::perform(double* const data)
{
    cuda_array_memcpy_d2d(_dataDouble, data, _resSize);

    cufftExecD2Z(_fftPlanFwd, (cufftDoubleReal *)    _dataDouble,
                                    (cufftDoubleComplex *) _dataDoubleFFT);

    dev_filterModulateAndNormalize
        <<<iDivUp(_fftSize, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (_dataDoubleFFT, _kernelDouble, _fftSize, (double) _resSize);

    cufftExecZ2D(_fftPlanInv, (cufftDoubleComplex *) _dataDoubleFFT,
                                    (cufftDoubleReal *)    _dataDouble);

    cuda_array_memcpy_d2d(data, _dataDouble, _resSize);
}

void cuFFTFilter::clean(const float* const kernel)
{
    require(kernel != 0);

    // single-precision
    cuda_arrayDev_delete(_kernelSingle);
    cuda_arrayDev_delete(_dataSingle);
    cuda_arrayDev_delete(_dataSingleFFT);

    CUFFT_SAFE_CALL(cufftDestroy(_fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(_fftPlanInv));
}

void cuFFTFilter::clean(const double* const kernel)
{
    require(kernel != 0);

    // double-precision
    cuda_arrayDev_delete(_kernelDouble);
    cuda_arrayDev_delete(_dataDouble);
    cuda_arrayDev_delete(_dataDoubleFFT);

    CUFFT_SAFE_CALL(cufftDestroy(_fftPlanFwd));
    CUFFT_SAFE_CALL(cufftDestroy(_fftPlanInv));
}

} // namespace gem
