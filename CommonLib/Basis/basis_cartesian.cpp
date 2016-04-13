/***********************************************************************
 *  File:       basis_cartesian.cpp
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

#include <cmath>
#include <tr1/cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/factorials.hpp>

namespace gem {

namespace bmath = boost::math;
namespace bmp   = boost::multiprecision;

/*****************************************
 * Cartesian basis
 ****************************************/

template <typename T>
void basis_cart_gaussian(const T* const arrayCoord, size_t length,
                               T* const arrayBasis,
                         unsigned int nmax, T sigma,
                         eBasisAccuracy accuracy)
{
    assert(arrayCoord != NULL);
    assert(arrayBasis != NULL);
    assert(length > 0);
    assert(sigma  > 0);

    unsigned int    iOrder;
    size_t          iLen;
    T               fPre, fExp1, fExp2;
    T               x, deltaX = 0, deltaX2 = 0, valTmp1 = 0, valTmp2 = 0;

    // make sure that deltaX is constant
    if (length > 1) {
        deltaX  = arrayCoord[1] - arrayCoord[0];
        deltaX2 = deltaX / 2;
    }

    switch (accuracy) {
        case BASIS_ACCU_APPROX:
            // using definition + recursive computation
            fExp1 = -1 / (2 * pow2(sigma));

            #pragma omp parallel for private(iOrder,x,fExp2)
            for (iLen = 0; iLen < length; iLen++) {
                x     = arrayCoord[iLen];
                fExp2 = std::exp(fExp1 * pow2(x));

                // iOrder = 0
                arrayBasis[iLen] = fExp2 * deltaX;

                // iOrder > 0
                for (iOrder = 1; iOrder <= nmax; iOrder++) {
                    arrayBasis[iOrder*length+iLen] = x / sigma *
                        arrayBasis[(iOrder-1)*length+iLen];
                }
            }
            break;
        case BASIS_ACCU_ACCURATE:
            // using error function + recursive computation
            #pragma omp parallel for private(fPre,iOrder,x,fExp1,valTmp1,valTmp2)
            for (iLen = 0; iLen < length; iLen++) {
                x     = arrayCoord[iLen];

                // iOrder = 0
                fPre    = std::sqrt((T) M_PI * pow2(sigma) / 2);
                fExp1   = 1 / (std::sqrt((T) 2) * sigma);
                valTmp1 = std::tr1::erf(fExp1 * (x-deltaX2));
                valTmp2 = std::tr1::erf(fExp1 * (x+deltaX2));

                arrayBasis[iLen] = fPre * (valTmp2 - valTmp1);

                // iOrder = 1
                if (nmax >= 1) {
                    fExp1   = -1 / (2 * pow2(sigma));
                    valTmp1 = std::exp(fExp1 * pow2(x-deltaX2));
                    valTmp2 = std::exp(fExp1 * pow2(x+deltaX2));

                    arrayBasis[length+iLen] = sigma * (valTmp1 - valTmp2);

                    // iOrder >= 2
                    fExp1   = -1 / (2 * pow2(sigma));
                    valTmp1 = sigma * std::exp(fExp1 * pow2(x-deltaX2));
                    valTmp2 = sigma * std::exp(fExp1 * pow2(x+deltaX2));

                    for (iOrder = 2; iOrder <= nmax; iOrder++) {
                        valTmp1 *= (x-deltaX2)/sigma;
                        valTmp2 *= (x+deltaX2)/sigma;

                        arrayBasis[iOrder*length+iLen] = valTmp1 - valTmp2 +
                            (T) (iOrder-1) * arrayBasis[(iOrder-2)*length+iLen];
                    }
                }
            }
            break;
        default:
            ERROR("basis_cart_gaussian", "unsupported accuracy");
    }
}

void basis_cart_gaussian(const floatmp* const arrayCoord, size_t length,
                               floatmp* const arrayBasis,
                         unsigned int nmax, floatmp sigma)
{
    assert(arrayCoord != NULL);
    assert(arrayBasis != NULL);
    assert(length > 0);
    assert(sigma  > 0);

    unsigned int    iOrder;
    size_t          iLen;
    floatmp         fPre, fExp1, fExp2;
    floatmp         x, deltaX = 0, deltaX2 = 0, valTmp1 = 0, valTmp2 = 0;

    // make sure that deltaX is constant
    if (length > 1) {
        deltaX  = arrayCoord[1] - arrayCoord[0];
        deltaX2 = deltaX / 2;
    }

    // using error function + recursive computation
    #pragma omp parallel for private(fPre,iOrder,x,fExp1,valTmp1,valTmp2)
    for (iLen = 0; iLen < length; iLen++) {
        x     = arrayCoord[iLen];

        // iOrder = 0
        fPre    = bmp::sqrt(bmath::constants::pi<floatmp>() * bmp::pow(sigma,2) / 2);
        fExp1   = 1 / (bmp::sqrt((floatmp) 2) * sigma);
        valTmp1 = bmath::erf(fExp1 * (x-deltaX2));
        valTmp2 = bmath::erf(fExp1 * (x+deltaX2));

        arrayBasis[iLen] = fPre * (valTmp2 - valTmp1);

        // iOrder = 1
        if (nmax >= 1) {
            fExp1   = -1 / (2 * bmp::pow(sigma,2));
            valTmp1 = bmp::exp(fExp1 * bmp::pow(x-deltaX2,2));
            valTmp2 = bmp::exp(fExp1 * bmp::pow(x+deltaX2,2));

            arrayBasis[length+iLen] = sigma * (valTmp1 - valTmp2);

            // iOrder >= 2
            fExp1   = -1 / (2 * bmp::pow(sigma,2));
            valTmp1 = sigma * bmp::exp(fExp1 * bmp::pow(x-deltaX2,2));
            valTmp2 = sigma * bmp::exp(fExp1 * bmp::pow(x+deltaX2,2));

            for (iOrder = 2; iOrder <= nmax; iOrder++) {
                valTmp1 *= (x-deltaX2)/sigma;
                valTmp2 *= (x+deltaX2)/sigma;

                arrayBasis[iOrder*length+iLen] = valTmp1 - valTmp2 +
                    (floatmp) (iOrder-1) * arrayBasis[(iOrder-2)*length+iLen];
            }
        }
    }
}

// instantiation
template
void basis_cart_gaussian<float >(const float*  const arrayCoord, size_t length,
                                       float*  const arrayBasis,
                                 unsigned int nmax, float  sigma,
                                 eBasisAccuracy accuracy);
template
void basis_cart_gaussian<double>(const double* const arrayCoord, size_t length,
                                       double* const arrayBasis,
                                 unsigned int nmax, double sigma,
                                 eBasisAccuracy accuracy);

template <typename T>
void basis_cart_hermite(const T* const arrayCoord, size_t length,
                              T* const arrayBasis,
                        unsigned int nmax, T sigma,
                        eBasisAccuracy accuracy)
{
    assert(arrayCoord != NULL);
    assert(arrayBasis != NULL);
    assert(length > 0);
    assert(sigma  > 0);

    unsigned int    iOrder, jOrder;
    size_t          iLen;
    T               fPre, fExp1, fExp2, fCval;
    T               x, deltaX = 0;
    T*              arrayGauss = NULL;

    // make sure that deltaX is constant
    if (length > 1) deltaX = arrayCoord[1] - arrayCoord[0];

    switch (accuracy) {
        case BASIS_ACCU_APPROX:
            // using the definition of "physicists" Hermite polynomials + recursive computation
            fExp1 = -1 / (2 * pow2(sigma));

            #pragma omp parallel for private(iOrder,x,fPre,fExp2)
            for (iLen = 0; iLen < length; iLen++) {
                x     = arrayCoord[iLen];
                fPre  = sigma * std::sqrt((T)M_PI);
                fExp2 = std::exp(fExp1 * pow2(x));

                /*// using std::tr1::hermite() for all iOrder
                for (iOrder = 0; iOrder <= nmax; iOrder++) {
                    fPre *= (iOrder == 0) ? 1 : (T) (2*iOrder);

                    arrayBasis[iOrder*length+iLen] = 1 / std::sqrt(fPre) * fExp2 *
                        std::tr1::hermite(iOrder,x/sigma) * deltaX;
                }*/

                // using recurrence relation
                // iOrder = 0
                arrayBasis[iLen] = 1 / std::sqrt(fPre) * fExp2 * deltaX;

                // iOrder = 1
                if (nmax >= 1) {
                    fPre *= 2;
                    arrayBasis[length+iLen] = 1 / std::sqrt(fPre) *
                        fExp2 * 2 * x / sigma * deltaX;
                }

                // iOrder >= 2
                for (iOrder = 2; iOrder <= nmax; iOrder++) {
                    arrayBasis[iOrder*length+iLen] = std::sqrt((T)2/(T)iOrder) *
                        x / sigma * arrayBasis[(iOrder-1)*length+iLen] -
                        std::sqrt((T)(iOrder-1)/(T)iOrder) *
                        arrayBasis[(iOrder-2)*length+iLen];
                }
            }
            break;
        case BASIS_ACCU_ACCURATE:
            // using exact computation of Gaussian functions in basis_cart_gaussian()
            array_new(arrayGauss, (nmax+1)*length);
            basis_cart_gaussian(arrayCoord, length,
                                arrayGauss,
                                nmax, sigma,
                                accuracy);

            //#pragma omp parallel for private(iOrder,jOrder,fPre,fCval)
            for (iLen = 0; iLen < length; iLen++) {
                fPre  = std::sqrt(1 / (sigma * std::sqrt((T)M_PI)));

                for (iOrder = 0; iOrder <= nmax; iOrder++) {
                    fPre   *= std::sqrt((iOrder == 0) ? 1 : (T) (2*iOrder));
                    arrayBasis[iOrder*length+iLen] = 0;

                    for (jOrder = 0; jOrder <= iOrder/2; jOrder++) {
                        fCval = fPre * std::pow((T)-1,(T)jOrder) /
                            ((T) factorial(jOrder) *
                             (T) factorial(iOrder-2*jOrder) *
                             std::pow((T)2,(T)(2*jOrder)));

                        arrayBasis[iOrder*length+iLen] += fCval *
                            arrayGauss[(iOrder-2*jOrder)*length+iLen];
                    }
                }
            }
            array_delete(arrayGauss);
            break;
        default:
            ERROR("basis_cart_hermite", "unsupported accuracy");
    }
}

void basis_cart_hermite(const floatmp* const arrayCoord, size_t length,
                              floatmp* const arrayBasis,
                        unsigned int nmax, floatmp sigma)
{
    assert(arrayCoord != NULL);
    assert(arrayBasis != NULL);
    assert(length > 0);
    assert(sigma  > 0);

    unsigned int    iOrder, jOrder;
    size_t          iLen;
    floatmp         fPre, fExp1, fExp2, fCval;
    floatmp         x, deltaX = 0;
    floatmp*        arrayGauss = NULL;

    // make sure that deltaX is constant
    if (length > 1) deltaX = arrayCoord[1] - arrayCoord[0];

    // using exact computation of Gaussian functions in basis_cart_gaussian()
    arrayGauss = new floatmp [(nmax+1)*length];
    basis_cart_gaussian(arrayCoord, length,
                        arrayGauss,
                        nmax, sigma);

    #pragma omp parallel for private(iOrder,jOrder,fPre,fCval)
    for (iLen = 0; iLen < length; iLen++) {
        fPre  = bmp::sqrt(1 / (sigma * bmp::sqrt(bmath::constants::pi<floatmp>())));

        for (iOrder = 0; iOrder <= nmax; iOrder++) {
            fPre   *= bmp::sqrt((iOrder == 0) ? 1 : (floatmp) (2*iOrder));
            arrayBasis[iOrder*length+iLen] = 0;

            for (jOrder = 0; jOrder <= iOrder/2; jOrder++) {
                fCval = fPre * std::pow(-1,jOrder) /
                    (bmath::factorial<floatmp>(jOrder) *
                     bmath::factorial<floatmp>(iOrder-2*jOrder) *
                     (floatmp) std::pow(2,2*jOrder));

                arrayBasis[iOrder*length+iLen] += fCval * arrayGauss[(iOrder-2*jOrder)*length+iLen];
            }
        }
    }

    delete [] arrayGauss;
}

// instantiation
template
void basis_cart_hermite<float >(const float*  const arrayCoord, size_t length,
                                      float*  const arrayBasis,
                                unsigned int nmax, float  sigma,
                                eBasisAccuracy accuracy);
template
void basis_cart_hermite<double>(const double* const arrayCoord, size_t length,
                                      double* const arrayBasis,
                                unsigned int nmax, double sigma,
                                eBasisAccuracy accuracy);

} // namespace gem
