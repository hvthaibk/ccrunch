/***********************************************************************
 *  File:       cData3.cpp
 *
 *  Purpose:    Implementation of a floating-point data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cData3.hpp"

#include "filter.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuData3.hpp"
#include "array.cuh"
#endif  // __GEM_USE_CUDA__

namespace gem {

template <typename T>
cData3<T>::cData3(const cSize3& size)
    : cArray3<T>(size)
{}

template <typename T>
cData3<T>::cData3(const cData3& other)
    : cArray3<T>(other)
{}


#ifdef __GEM_USE_CUDA__
template <typename T>
cData3<T>& cData3<T>::operator=(const cuData3<T>& other)
{
    cArray3<T>::memReAlloc(other.getSize());

    if (cArray3<T>::getNelement() > 0) {
        cuda_array_memcpy_d2h(cArray3<T>::_data, other.getAddrData(), cArray3<T>::getNelement());
    }

    return *this;
}
#endif  // __GEM_USE_CUDA__

template <typename T>
bool cData3<T>::isFinite(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::isFinite(void)");

    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        if (!std::isfinite(cArray3<T>::_data[i])) {
            return false;
        }
    }

    return true;
}

template <typename T>
bool cData3<T>::isInf(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::isInf(void)");

    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        if (std::isinf(cArray3<T>::_data[i])) {
            return true;
        }
    }

    return false;
}

template <typename T>
bool cData3<T>::isNaN(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::isNaN(void)");

    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        if (std::isnan(cArray3<T>::_data[i])) {
            return true;
        }
    }

    return false;
}

template <typename T>
void cData3<T>::ceil(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::ceil(void)");

    #pragma omp parallel for
    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        cArray3<T>::_data[i] = std::ceil(cArray3<T>::_data[i]);
    }
}

template <typename T>
void cData3<T>::floor(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::floor(void)");

    #pragma omp parallel for
    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        cArray3<T>::_data[i] = std::floor(cArray3<T>::_data[i]);
    }
}

template <typename T>
void cData3<T>::round(void)
{
    cArray3<T>::requireNonEmpty("void cData3<T>::round(void)");

    #pragma omp parallel for
    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        cArray3<T>::_data[i] = gem::round(cArray3<T>::_data[i]);
    }
}

template <typename T>
void cData3<T>::linspace(T low, T high, size_t npoint)
{
    const std::string   funcName("void cData3<T>::linspace("
                                    "T low, T high, size_t npoint)");

    require(npoint > 1, funcName + ": only one point is requested");

    T       step = (high-low) / (T) (npoint-1);

    cArray3<T>::memReAlloc(cSize3(npoint,1,1));

    for (size_t i = 0; i < npoint; i++) {
        cArray3<T>::_data[i] = low + (T) i * step;
    }
}

// instantiation
template class cData3<float >;
template class cData3<double>;

} // namespace gem
