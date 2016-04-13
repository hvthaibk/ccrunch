/***********************************************************************
 *  File:       cuFilterFrequency.hpp
 *
 *  Purpose:    Header file for a frequency filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUFILTER_FREQUENCY_HPP__
#define __GEM_CUFILTER_FREQUENCY_HPP__

#include "cFilterFrequency.hpp"
#include "cuData3.hpp"

namespace gem {

template <typename T> class cFilterFrequency;

template <typename T>
class cuFilterFrequency
{
private:
    cuData3<T>      _kernel;

public:
     cuFilterFrequency(void) {};
    ~cuFilterFrequency()     {};

    void    memFree(void);
    void    copyKernelToGPU(const cFilterFrequency<T>& other);

    // spatial filtering
    void    filterFreq2(cuData3<T>& data);

};

} // namespace gem

#endif
