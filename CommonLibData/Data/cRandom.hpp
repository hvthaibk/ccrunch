/***********************************************************************
 *  File:       cRandom.hpp
 *
 *  Purpose:    Header file for a random class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CRANDOM_HPP__
#define __GEM_CRANDOM_HPP__

#include "cData3.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuData3.hpp"
#endif

namespace gem {

template <typename T>
class cRandom
{
private:

public:
     cRandom(void) {};
    ~cRandom()     {};

    // normally distributed pseudorandom numbers
    void    randn(cData3<T>& data, T mean, T std, const cSize3& size);

    // uniformly distributed pseudorandom numbers
    void    rand (cData3<T>& data, T min,  T max, const cSize3& size);

    // uniformly distributed pseudorandom numbers for zero-valued elements
    void    randZero(cData3<T>& data, T min, T max, T threshold);

#ifdef __GEM_USE_CUDA__

    // normally distributed pseudorandom numbers
    void    randn(cuData3<T>& data, T mean, T std, const cSize3& size, size_t seed = 32);

    // uniformly distributed pseudorandom numbers
    void    rand (cuData3<T>& data, T min,  T max, const cSize3& size, size_t seed = 32);

    // uniformly distributed pseudorandom numbers for zero-valued elements
    void    randZero(cuData3<T>& data, T min, T max, T threshold);

#endif
};

} // namespace gem

#endif
