/***********************************************************************
 *  File:       cuData3.hpp
 *
 *  Purpose:    Header file for a floating-point data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUDATA3_HPP__
#define __GEM_CUDATA3_HPP__

#include "transform.cuh"
#include "cuArray3.hpp"
#include "cData3.hpp"

namespace gem {

template <typename T>
class cuData3 : public cuArray3<T>
{
public:
     cuData3(const cSize3& size = cSize3(0,0,0));
     cuData3(const cuData3& other);
    ~cuData3() {};

    // operator overloading
    cuData3&    operator=(const cData3<T>& other);
};

typedef     cuData3<float >     cuData3f;
typedef     cuData3<double>     cuData3d;

} // namespace gem

#endif
