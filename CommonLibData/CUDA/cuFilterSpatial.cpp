/***********************************************************************
 *  File:       cuFilterSpatial.cpp
 *
 *  Purpose:    Implementation of a spatial filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuFilterSpatial.hpp"

#include "xcorr.cuh"

namespace gem {

template <typename T>
void cuFilterSpatial<T>::memFree(void)
{
    _kernel.memFree();
}

template <typename T>
void cuFilterSpatial<T>::copyKernelToGPU(const cFilterSpatial<T>& other)
{
    const std::string   funcName("void cuFilterFrequency<T>::copyKernelToGPU("
                                    "const cFilterFrequency<T>& other)");

    other._kernel.requireNonEmpty(funcName);

    _kernel = other._kernel;
}

template <typename T>
void cuFilterSpatial<T>::filterSpatial(cuData3<T>& data, eXCorrRes shape)
{
    const std::string   funcName("void cuFilterSpatial<T>::filterSpatial("
                                    "cuData3<T>& data, eXCorrRes shape);");

       data.requireNonEmpty(funcName);
    _kernel.requireNonEmpty(funcName);

    require(data.getDimension() == _kernel.getDimension(),
            funcName + ": data and kernel have different dimensions");

    cuData3<T>    other(data.getSize() + _kernel.getSize() - 1);

    switch (data.getDimension()) {
        case 1:
            cuda_xcorr(_kernel.getAddrData(), _kernel.getNrow(),
                          data.getAddrData(),    data.getNrow(),
                         other.getAddrData(),
                       XCORR_RES_FULL);
            break;
        case 2:
            cuda_xcorr(_kernel.getAddrData(), _kernel.getNrow(), _kernel.getNcol(),
                          data.getAddrData(),    data.getNrow(),    data.getNcol(),
                         other.getAddrData(),
                       XCORR_RES_FULL);
            break;
        case 3:
            cuda_xcorr(_kernel.getAddrData(), _kernel.getNrow(), _kernel.getNcol(), _kernel.getNsec(),
                          data.getAddrData(),    data.getNrow(),    data.getNcol(),    data.getNsec(),
                         other.getAddrData(),
                       XCORR_RES_FULL);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    // cropping
    switch (shape) {
        case XCORR_RES_FULL:
            data = other;
            break;
        case XCORR_RES_SAME:
            cuda_array_crop(other.getAddrData(), other.getNrow(), other.getNcol(), other.getNsec(),
                             data.getAddrData(),  data.getNrow(),  data.getNcol(),  data.getNsec(),
                            (_kernel.getNrow()-1)/2,
                            (_kernel.getNcol()-1)/2,
                            (_kernel.getNsec()-1)/2);
            break;
        case XCORR_RES_VALID:
            ERROR(funcName, "unsupported xcorr mode");
            break;
        default:
            ERROR(funcName, "unsupported xcorr mode");
    }
}

// instantiation
template class cuFilterSpatial<float >;
template class cuFilterSpatial<double>;

} // namespace gem
