/***********************************************************************
 *  File:       cBasisPolar.hpp
 *
 *  Purpose:    Header file for a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CBASIS_POLAR_HPP__
#define __GEM_CBASIS_POLAR_HPP__

#include "cBasis.hpp"

namespace gem {

template <typename T>
class cBasisPolar : public cBasis<T>
{
private:
    cData3<T>       _basisRad, _basisThe;

public:
     cBasisPolar(void);
    ~cBasisPolar()    ;
};

} // namespace gem

#endif
