/***********************************************************************
 *  File:       cuFilterSpatial.hpp
 *
 *  Purpose:    Header file for a spatial filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUFILTER_SPATIAL_HPP__
#define __GEM_CUFILTER_SPATIAL_HPP__

#include "xcorr.hpp"
#include "cFilterSpatial.hpp"
#include "cuData3.hpp"

namespace gem {

template <typename T> class cFilterSpatial;

template <typename T>
class cuFilterSpatial
{
private:
    cuData3<T>      _kernel;

public:
     cuFilterSpatial(void) {};
    ~cuFilterSpatial()     {};

    void    memFree(void);
    void    copyKernelToGPU(const cFilterSpatial<T>& other);

    // spatial filtering
    void    filterSpatial(cuData3<T>& data, eXCorrRes shape = XCORR_RES_SAME);

};

} // namespace gem

#endif
