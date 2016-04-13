/***********************************************************************
 *  File:       cBasisSpherical.hpp
 *
 *  Purpose:    Header file for a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CBASIS_SPHERICAL_HPP__
#define __GEM_CBASIS_SPHERICAL_HPP__

#include "cBasis.hpp"

namespace gem {

template <typename T>
class cBasisSpherical : public cBasis<T>
{
private:
    cData3<T>       _basisRad, _basisThe, _basisPhi;

public:
     cBasisSpherical(void);
    ~cBasisSpherical()    ;
};

} // namespace gem

#endif
