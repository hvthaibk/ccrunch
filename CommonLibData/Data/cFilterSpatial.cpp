/***********************************************************************
 *  File:       cFilterSpatial.cpp
 *
 *  Purpose:    Implementation of a spatial filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFilterSpatial.hpp"

#include "filter.hpp"

namespace gem {

template <typename T>
void cFilterSpatial<T>::memFree(void)
{
    _kernel.memFree();
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialAverage(size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialAverage("
                                    "size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_average(_kernel.getAddrData(), _kernel.getNrow());
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_average(_kernel.getAddrData(), _kernel.getNrow(),
                                                  _kernel.getNcol());
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_average(_kernel.getAddrData(), _kernel.getNrow(),
                                                  _kernel.getNcol(),
                                                  _kernel.getNsec());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialDisk(T radius, size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialDisk("
                                    "T radius, size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_disk(_kernel.getAddrData(), _kernel.getNrow(), radius, true);
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_disk(_kernel.getAddrData(), _kernel.getNrow(),
                                               _kernel.getNcol(), radius, true);
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_disk(_kernel.getAddrData(), _kernel.getNrow(),
                                               _kernel.getNcol(),
                                               _kernel.getNsec(), radius, true);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialGaussian(T sigma, size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialGaussian("
                                    "T sigma, size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_gaussian(_kernel.getAddrData(), _kernel.getNrow(), sigma);
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_gaussian(_kernel.getAddrData(), _kernel.getNrow(),
                                                   _kernel.getNcol(), sigma);
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_gaussian(_kernel.getAddrData(), _kernel.getNrow(),
                                                   _kernel.getNcol(),
                                                   _kernel.getNsec(), sigma);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialLaplacian(T alpha, T beta, size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialLaplacian("
                                    "T alpha, T beta, size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_laplacian(_kernel.getAddrData(), _kernel.getNrow());
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_laplacian(_kernel.getAddrData(), _kernel.getNrow(),
                                                    _kernel.getNcol(), alpha);
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_laplacian(_kernel.getAddrData(), _kernel.getNrow(),
                                                    _kernel.getNcol(),
                                                    _kernel.getNsec(), alpha, beta);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialDoG(T sigma1, T sigma2, size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialDoG("
                                    "T sigma1, T sigma2, size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_dog(_kernel.getAddrData(), _kernel.getNrow(), sigma1, sigma2);
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_dog(_kernel.getAddrData(), _kernel.getNrow(),
                                              _kernel.getNcol(), sigma1, sigma2);
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_dog(_kernel.getAddrData(), _kernel.getNrow(),
                                              _kernel.getNcol(),
                                              _kernel.getNsec(), sigma1, sigma2);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::kernelSpatialLoG(T sigma, size_t extend, size_t dim)
{
    const std::string   funcName("void cFilterSpatial<T>::kernelSpatialLoG("
                                    "T sigma, size_t extend, size_t dim)");

    switch (dim) {
        case 1:
            _kernel.memReAlloc((cSize3(extend,1,1)));

            filter_log(_kernel.getAddrData(), _kernel.getNrow(), sigma);
            break;
        case 2:
            _kernel.memReAlloc((cSize3(extend,extend,1)));

            filter_log(_kernel.getAddrData(), _kernel.getNrow(),
                                              _kernel.getNcol(), sigma);
            break;
        case 3:
            _kernel.memReAlloc((cSize3(extend,extend,extend)));

            filter_log(_kernel.getAddrData(), _kernel.getNrow(),
                                              _kernel.getNcol(),
                                              _kernel.getNsec(), sigma);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cFilterSpatial<T>::filterSpatial(cData3<T>& data, eXCorrRes shape)
{
    const std::string   funcName("void cFilterSpatial<T>::filterSpatial("
                                    "cData3<T>& data, eXCorrRes shape)");

      data.requireNonEmpty(funcName);
    _kernel.requireNonEmpty(funcName);

    require(data.getDimension() == _kernel.getDimension(),
            funcName + ": data and _kernel have different dimensions");

    cData3<T>    other(data.getSize() + _kernel.getSize() - 1);

    switch (data.getDimension()) {
        case 1:
            xcorr(_kernel.getAddrData(), _kernel.getNrow(),
                     data.getAddrData(),    data.getNrow(),
                    other.getAddrData(),
                  XCORR_RES_FULL);
            break;
        case 2:
            xcorr(_kernel.getAddrData(), _kernel.getNrow(), _kernel.getNcol(),
                     data.getAddrData(),    data.getNrow(),   data.getNcol(),
                    other.getAddrData(),
                  XCORR_RES_FULL);
            break;
        case 3:
            xcorr(_kernel.getAddrData(), _kernel.getNrow(), _kernel.getNcol(), _kernel.getNsec(),
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
            array_crop(other.getAddrData(), other.getNrow(), other.getNcol(), other.getNsec(),
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
template class cFilterSpatial<float >;
template class cFilterSpatial<double>;

} // namespace gem
