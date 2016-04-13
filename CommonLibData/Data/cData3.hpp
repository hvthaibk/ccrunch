/***********************************************************************
 *  File:       cData3.hpp
 *
 *  Purpose:    Header file for a floating-point data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CDATA3_HPP__
#define __GEM_CDATA3_HPP__

#include "cArray3.hpp"

#include <iostream>
#include <fstream>

namespace gem {

#ifdef __GEM_USE_CUDA__
template <typename T>
class cuData3;
#endif

template <typename T>
class cData3 : public cArray3<T>
{
public:
     cData3(const cSize3& size = cSize3(0,0,0));
     cData3(const cData3& other);
    ~cData3() {};

    // operator overloading
#ifdef __GEM_USE_CUDA__
    cData3&     operator=(const cuData3<T>& other);
#endif

    // validity
    bool    isFinite(void);
    bool    isInf   (void);
    bool    isNaN   (void);

    // interger approximation
    void    ceil (void);
    void    floor(void);
    void    round(void);

    // linspaced
    void    linspace(T low, T high, size_t npoint);
};

typedef     cData3<float >      cData3f;
typedef     cData3<double>      cData3d;

} // namespace gem

#endif
