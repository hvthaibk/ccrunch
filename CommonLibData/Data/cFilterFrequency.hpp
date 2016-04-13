/***********************************************************************
 *  File:       cFilterFrequency.hpp
 *
 *  Purpose:    Header file for a frequency filter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CFILTER_FREQUENCY_HPP__
#define __GEM_CFILTER_FREQUENCY_HPP__

#include "cData3.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuFilterFrequency.hpp"
#endif

namespace gem {

template <typename T>
class cFilterFrequency
{
#ifdef __GEM_USE_CUDA__
template<typename U> friend class cuFilterFrequency;
#endif

private:
    cData3<T>       _kernel;

public:
     cFilterFrequency(void) {};
    ~cFilterFrequency()     {};

    void    memFree(void);

    // kernels
    void    kernelReverse          (void);
    void    kernelFreq2Distance    (const cSize3& size);
    void    kernelFreq2Ring        (const cSize3& size, T radius);
    void    kernelFreq2NotchLineHor(const cSize3& size, T radius, T sigma);
    void    kernelFreq2NotchLineVer(const cSize3& size, T radius, T sigma);
    void    kernelFreq2LowPass     (const cSize3& size, T sigma);
    void    kernelFreq2HighPass    (const cSize3& size, T sigma);

    cData3<T>  getKernelFull2(const cSize3& size);
    cData3<T>& getKernel     (void)                    {  return _kernel;  };
    void       setKernel     (const cData3<T>& kernel) { _kernel = kernel; };

    // frequency filtering
    void    filterFreq2(cData3<T>& data);
};

} // namespace gem

#endif
