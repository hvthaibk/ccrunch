/***********************************************************************
 *  File:       basis_conversion.cpp
 *
 *  Purpose:    Implementation of basis-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "basis.hpp"

#include <boost/math/special_functions/factorials.hpp>

namespace gem {

namespace bmath = boost::math;
namespace bmp   = boost::multiprecision;

/*****************************************
 * Basis conversion
 ****************************************/

template <typename T>
void basis_convert_G2H1D(unsigned int nmax, T sigma,
                         T* const arrayCoef);

void basis_convert_G2H1D(unsigned int nmax, floatmp sigma,
                         floatmp* const arrayCoef);

} // namespace gem
