/***********************************************************************
 *  File:       cFilterSpatial.hpp
 *
 *  Purpose:    Header file for a spatial filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CFILTER_SPATIAL_HPP__
#define __GEM_CFILTER_SPATIAL_HPP__

#include "xcorr.hpp"
#include "cData3.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuFilterSpatial.hpp"
#endif

namespace gem {

template <typename T>
class cFilterSpatial
{
#ifdef __GEM_USE_CUDA__
template<typename U> friend class cuFilterSpatial;
#endif

private:
    cData3<T>       _kernel;

public:
     cFilterSpatial(void) {};
    ~cFilterSpatial()     {};

    void    memFree(void);

    // kernels
    void    kernelSpatialAverage  (                    size_t extend, size_t dim);
    void    kernelSpatialDisk     (T radius,           size_t extend, size_t dim);
    void    kernelSpatialGaussian (T sigma,            size_t extend, size_t dim);
    void    kernelSpatialLaplacian(T alpha,  T beta,   size_t extend, size_t dim);
    void    kernelSpatialDoG      (T sigma1, T sigma2, size_t extend, size_t dim);
    void    kernelSpatialLoG      (T sigma,            size_t extend, size_t dim);

    // spatial filtering
    void    filterSpatial(cData3<T>& data, eXCorrRes shape = XCORR_RES_SAME);
};

} // namespace gem

#endif
