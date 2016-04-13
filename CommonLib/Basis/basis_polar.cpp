/***********************************************************************
 *  File:       basis_polar.cpp
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

namespace gem {

/*****************************************
 * Polar basis
 ****************************************/

template <typename T>
void basis_polar_Rn(const T* const arrayCoord,  size_t       length,
                          T* const arrayRadial, unsigned int orderMax,
                          T sigma,
                    eBasisFunction method)
{
    assert(arrayCoord  != NULL);
    assert(arrayRadial != NULL);
    assert(length > 0);
    assert(sigma  > 0);

    switch (method) {
        case BASIS_POLAR_RN:
            for (size_t i = 0; i < length; i++) {
            }
            break;
        default:
            ERROR("basis_polar_Rn", "unsupported basis function");
    }
}

// instantiation
template
void basis_polar_Rn<float >(const float*  const arrayCoord,  size_t       length,
                                  float*  const arrayRadial, unsigned int orderMax,
                                  float  sigma,
                            eBasisFunction method);
template
void basis_polar_Rn<double>(const double* const arrayCoord,  size_t       length,
                                  double* const arrayRadial, unsigned int orderMax,
                                  double sigma,
                            eBasisFunction method);

template <typename T>
void basis_spher_Am(const T*         const arrayCoord,   size_t       length,
                    std::complex<T>* const arrayAngular, unsigned int orderMax)
{
    assert(arrayCoord   != NULL);
    assert(arrayAngular != NULL);
    assert(length > 0);
}

// instantiation
template
void basis_spher_Am<float >(const float*          const arrayCoord,   size_t       length,
                            std::complex<float>*  const arrayAngular, unsigned int orderMax);
template
void basis_spher_Am<double>(const double*         const arrayCoord,   size_t       length,
                            std::complex<double>* const arrayAngular, unsigned int orderMax);

} // namespace gem
