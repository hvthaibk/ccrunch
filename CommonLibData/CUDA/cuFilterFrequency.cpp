/***********************************************************************
 *  File:       cuFilterFrequency.cpp
 *
 *  Purpose:    Implementation of a frequency filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuFilterFrequency.hpp"

#include "cuFFTFilter.cuh"

namespace gem {

template <typename T>
void cuFilterFrequency<T>::memFree(void)
{
    _kernel.memFree();
}

template <typename T>
void cuFilterFrequency<T>::copyKernelToGPU(const cFilterFrequency<T>& other)
{
    const std::string   funcName("void cuFilterFrequency<T>::copyKernelToGPU("
                                    "const cFilterFrequency<T>& other)");

    other._kernel.requireNonEmpty(funcName);

    _kernel = other._kernel;
}

template <typename T>
void cuFilterFrequency<T>::filterFreq2(cuData3<T>& data)
{
    const std::string   funcName("void cuFilterFrequency<T>::filterFreq2("
                                    "cuData3<T>& data);");

    _kernel.requireNonEmpty(funcName);
       data.requireNonEmpty(funcName);

    cuFFTFilter     objFilter;

    objFilter.prepare(_kernel.getAddrData(), data.getNrow(), data.getNcol());
    objFilter.perform(   data.getAddrData());
    objFilter.clean  (_kernel.getAddrData());
}

// instantiation
template class cuFilterFrequency<float >;
template class cuFilterFrequency<double>;

} // namespace gem
