/***********************************************************************
 *  File:       cuData3.cpp
 *
 *  Purpose:    Implementation of a floating-point data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuData3.hpp"

namespace gem {

template <typename T>
cuData3<T>::cuData3(const cSize3& size)
    : cuArray3<T>(size)
{}

template <typename T>
cuData3<T>::cuData3(const cuData3& other)
    : cuArray3<T>(other)
{}

template <typename T>
cuData3<T>& cuData3<T>::operator=(const cData3<T>& other)
{
    cuArray3<T>::memReAlloc(other.getSize());

    if (cuArray3<T>::getNelement() > 0) {
        cuda_array_memcpy_h2d(cuArray3<T>::_data, other.getAddrData(), cuArray3<T>::getNelement());
    }

    return *this;
}

// instantiation
template class cuData3<float >;
template class cuData3<double>;

} // namespace gem
