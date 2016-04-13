/***********************************************************************
 *  File:       cuTransThreshold.cpp
 *
 *  Purpose:    Implementation of a thresholding transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuTransThreshold.hpp"

namespace gem {

template <typename T>
void cuTransThreshold<T>::threshold(cuData3<T>& dataSrc, T cutoff)
{
    const std::string   funcName("void cuTransThreshold<T>::threshold(cuData3<T>& dataSrc, T cutoff)");

    dataSrc.requireNonEmpty(funcName);

    cuda_array_threshold(dataSrc.getAddrData(), dataSrc.getNelement(), cutoff);
}

template <typename T>
void cuTransThreshold<T>::threshold(const cuData3<T>& dataSrc, cuData3<T>& dataDst, T cutoff)
{
    const std::string   funcName("void cuTransThreshold<T>::threshold(const cuData3<T>& dataSrc, cuData3<T>& dataDst, T cutoff)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    cuda_array_threshold(dataSrc.getAddrData(), dataSrc.getNelement(), cutoff,
                         dataDst.getAddrData());
}

// instantiation
template class cuTransThreshold<float >;
template class cuTransThreshold<double>;

} // namespace gem
