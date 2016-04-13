/***********************************************************************
 *  File:       cFilterFrequency.cpp
 *
 *  Purpose:    Implementation of a frequency filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFilterFrequency.hpp"

#include "cNDgrid.hpp"
#include "cCoordinate.hpp"
#include "cTransThreshold.hpp"
#include "cFFTFilter.hpp"

namespace gem {

template <typename T>
void cFilterFrequency<T>::memFree(void)
{
    _kernel.memFree();
}

template <typename T>
void cFilterFrequency<T>::kernelReverse(void)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelReverse(void)");

    _kernel.requireNonEmpty(funcName);

    #pragma omp parallel for
    for (size_t i = 0; i < _kernel.getNelement(); i++) {
        _kernel[i] = 1 - _kernel[i];
    }
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2Distance(const cSize3& size)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2Distance("
                                    "const cSize3& size)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    cData3<T>   R, C;

    cNDgrid<T>().genRCfreq(size[0], size[1], R, C);
    cCoordinate<T>().computeDistance(C, R, _kernel);
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2Ring(const cSize3& size, T radius)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2Ring("
                                    "const cSize3& size, T radius)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    kernelFreq2Distance(size);

    cTransThreshold<T>().threshold(_kernel, radius);

    kernelReverse();
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2NotchLineHor(const cSize3& size, T radius, T sigma)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2NotchLineHor("
                                    "const cSize3& size, T radius, T sigma)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    T           fExp = -1 / (2 * pow2(sigma));
    cData3<T>   R, C, kernel1, kernel2;

    kernelFreq2Ring(size, radius);
    kernel1 = _kernel;
    kernelFreq2HighPass(size, 3*radius);
    kernel2 = _kernel;

    cNDgrid<T>().genRCfreq(size[0], size[1], R, C);

    #pragma omp parallel for
    for (size_t i = 0 ; i < _kernel.getNelement(); i++) {
        _kernel[i] = 1 - (std::exp(fExp * pow2(R[i]))) *
                         (1 - kernel1[i]) * kernel2[i];
    }
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2NotchLineVer(const cSize3& size, T radius, T sigma)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2NotchLineVer("
                                    "const cSize3& size, T radius, T sigma)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    T           fExp = -1 / (2 * pow2(sigma));
    cData3<T>   R, C, kernel1, kernel2;

    kernelFreq2Ring(size, radius);
    kernel1 = _kernel;
    kernelFreq2HighPass(size, 3*radius);
    kernel2 = _kernel;

    cNDgrid<T>().genRCfreq(size[0], size[1], R, C);

    #pragma omp parallel for
    for (size_t i = 0 ; i < _kernel.getNelement(); i++) {
        _kernel[i] = 1 - (std::exp(fExp * pow2(C[i]))) *
                         (1 - kernel1[i]) * kernel2[i];
    }
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2LowPass(const cSize3& size, T sigma)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2LowPass("
                                    "const cSize3& size, T sigma)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    T           fExp = -1 / (2 * pow2(sigma));
    cData3<T>   R, C;

    cNDgrid<T>().genRCfreq(size[0], size[1], R, C);

    #pragma omp parallel for
    for (size_t i = 0 ; i < _kernel.getNelement(); i++) {
        _kernel[i] = std::exp(fExp * (pow2(R[i]) + pow2(C[i])));
    }
}

template <typename T>
void cFilterFrequency<T>::kernelFreq2HighPass(const cSize3& size, T sigma)
{
    const std::string   funcName("void cFilterFrequency<T>::kernelFreq2HighPass("
                                    "const cSize3& size, T sigma)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1,
            funcName + ": incorrect size");

    kernelFreq2LowPass(size, sigma);

    kernelReverse();
}

template <typename T>
cData3<T> cFilterFrequency<T>::getKernelFull2(const cSize3& size)
{
    const std::string   funcName("void cFilterFrequency<T>::getKernelFull2("
                                        "const cSize3& size, "
                                        "cData3<T>& kernelFull)");

    _kernel.requireNonEmpty(funcName);

    cData3<T>   kernelTmp, kernelFull;
    size_t      nrowHalf = _kernel.getNcol();

    kernelTmp.opCircShift(_kernel, cOffset3(size[0]/2,0,0));

    kernelFull.memReAllocZero(size);

    // the right half
    if (((size[1] >> 1) << 1) == size[1]) {
        kernelFull.opReplace(kernelTmp, cSize3(0,size[1]/2-1,0));
    }
    else {
        kernelFull.opReplace(kernelTmp, cSize3(0,size[1]/2  ,0));
    }

    // the left half
    if (((size[1] >> 1) << 1) == size[1]) {
        for (size_t irow = 0; irow < size[0]; irow++) {
            for (size_t icol = 0; icol < size[1]/2-1; icol++) {
                kernelFull[irow*size[1]+icol] = kernelTmp[irow*nrowHalf+nrowHalf-2-icol];
            }
        }
    }
    else {
        for (size_t irow = 0; irow < size[0]; irow++) {
            for (size_t icol = 0; icol < size[1]/2; icol++) {
                kernelFull[irow*size[1]+icol] = kernelTmp[irow*nrowHalf+nrowHalf-1-icol];
            }
        }
    }

    if (((size[1] >> 1) << 1) == size[1]) {
        kernelFull.opCircShift(cOffset3(0,1,0));
    }

    return kernelFull;
}

template <typename T>
void cFilterFrequency<T>::filterFreq2(cData3<T>& data)
{
    const std::string   funcName("void cFilterFrequency<T>::filterFreq("
                                    "cData3<T>& data)");

    _kernel.requireNonEmpty(funcName);
       data.requireNonEmpty(funcName);

    cFFTFilter      objFilter;

    objFilter.prepare(_kernel.getAddrData(), data.getNrow(), data.getNcol());
    objFilter.perform(   data.getAddrData());
    objFilter.clean  (_kernel.getAddrData());
}

// instantiation
template class cFilterFrequency<float >;
template class cFilterFrequency<double>;

} // namespace gem
