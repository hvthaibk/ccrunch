/***********************************************************************
 *  File:       basis_spherical.cpp
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
 * Spherical basis
 ****************************************/

template <typename T>
void basis_spher_Rln(const T* const arrayCoord,  size_t       length,
                           T* const arrayRadial, unsigned int orderMax,
                           T sigma,
                     eBasisFunction method);

template <typename T>
void basis_spher_Ylm(const T*         const arrayCoord,   size_t       length,
                     std::complex<T>* const arrayAngular, unsigned int orderMax);

template <typename T>
void basis_spher_ylm(const T* const arrayCoord,   size_t       length,
                           T* const arrayAngular, unsigned int orderMax);

/*template <typename T>
void basis_Rln(T* const arrayCoord,  size_t       length,
               T* const arrayRadial, unsigned int orderMax,
               eBasisSphericalRadial type)
{
    assert(arrayCoord  != NULL);
    assert(arrayRadial != NULL);
    assert(length != 0);
}

// instantiation
template
void basis_Rln<float >(float*  const arrayCoord,  size_t       length,
                       float*  const arrayRadial, unsigned int orderMax,
                       eBasisSphericalRadial type);
template
void basis_Rln<double>(double* const arrayCoord,  size_t       length,
                       double* const arrayRadial, unsigned int orderMax,
                       eBasisSphericalRadial type);

template <typename T>
void basis_Ylm(T*               const arrayCoord,   size_t       length,
               std::complex<T>* const arrayAngular, unsigned int orderMax)
{
    assert(arrayCoord   != NULL);
    assert(arrayAngular != NULL);
    assert(length != 0);
}

// instantiation
template
void basis_Ylm<float >(float*                const arrayCoord,   size_t       length,
                       std::complex<float >* const arrayAngular, unsigned int orderMax);
template
void basis_Ylm<double>(double*               const arrayCoord,   size_t       length,
                       std::complex<double>* const arrayAngular, unsigned int orderMax);

template <typename T>
void basis_ylm(T* const arrayCoord,   size_t       length,
               T* const arrayAngular, unsigned int orderMax)
{
    assert(arrayCoord   != NULL);
    assert(arrayAngular != NULL);
    assert(length != 0);
}

// instantiation
template
void basis_ylm<float >(float*  const arrayCoord,   size_t       length,
                       float*  const arrayAngular, unsigned int orderMax);
template
void basis_ylm<double>(double* const arrayCoord,   size_t       length,
                       double* const arrayAngular, unsigned int orderMax);*/

} // namespace gem
