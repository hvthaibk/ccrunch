/***********************************************************************
 *  File:       cBasisGaussian.cpp
 *
 *  Purpose:    Implementation of a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cBasisGaussian.hpp"

namespace gem {

template <typename T>
cBasisGaussian<T>::cBasisGaussian(void)
{
}

template <typename T>
cBasisGaussian<T>:: ~cBasisGaussian()
{
    memFree();
}

template <typename T>
void cBasisGaussian<T>::memFree(void)
{
    _basisRow.memFree();
    _basisCol.memFree();
    _basisSec.memFree();
}

template <typename T>
void cBasisGaussian<T>::printData(std::string message, int width, int precision, bool scientific) const
{
    if (_basisRow.getNelement() > 0) {
        _basisRow.printData(message, width, precision, scientific);
    }

    if (_basisCol.getNelement() > 0) {
        _basisCol.printData(message, width, precision, scientific);
    }

    if (_basisSec.getNelement() > 0) {
        _basisSec.printData(message, width, precision, scientific);
    }
}

template <typename T>
void cBasisGaussian<T>::printSize(std::string message, int width, int precision) const
{
    _basisRow.printSize(message, width, precision);
    _basisCol.printSize(message, width, precision);
    _basisSec.printSize(message, width, precision);
}

template <typename T>
void cBasisGaussian<T>::genBasisRow(T sigma, eBasisAccuracy accuracy)
{
    const std::string   funcName("void cBasisGaussian<T>::genBasisRow("
                                    "T sigma, eBasisAccuracy accuracy)");

    cBasis<T>::_coordRow.requireNonEmpty(funcName);

    _basisRow.memReAlloc(cSize3(cBasis<T>::_nmax+1,
                         cBasis<T>::_coordRow.getNrow(),
                         1));

    switch (accuracy) {
        case BASIS_ACCU_APPROX:
        case BASIS_ACCU_ACCURATE:
            basis_cart_gaussian(cBasis<T>::_coordRow.getAddrData(), cBasis<T>::_coordRow.getNrow(),
                                           _basisRow.getAddrData(),
                                cBasis<T>::_nmax, sigma,
                                accuracy);
            break;
        case BASIS_ACCU_MP:
            floatmp     *arrayCoord, *arrayBasis;

            arrayCoord = new floatmp [cBasis<T>::_coordRow.getNrow()];
            arrayBasis = new floatmp [cBasis<T>::_coordRow.getNrow()*(cBasis<T>::_nmax+1)];

            for (size_t i = 0; i < cBasis<T>::_coordRow.getNrow(); i++) {
                arrayCoord[i] = cBasis<T>::_coordRow[i];
            }

            basis_cart_gaussian(arrayCoord, cBasis<T>::_coordRow.getNrow(),
                                arrayBasis,
                                cBasis<T>::_nmax, (floatmp) sigma);

            for (size_t i = 0; i < cBasis<T>::_coordRow.getNrow()*(cBasis<T>::_nmax+1); i++) {
                _basisRow[i] = arrayBasis[i].convert_to<T>();
            }

            delete [] arrayCoord;
            delete [] arrayBasis;

            break;
        default:
            ERROR(funcName, "unsupported accuracy");
    }
}

template <typename T>
void cBasisGaussian<T>::genBasisCol(T sigma, eBasisAccuracy accuracy)
{
    const std::string   funcName("void cBasisGaussian<T>::genBasisCol("
                                    "T sigma, eBasisAccuracy accuracy)");

    cBasis<T>::_coordCol.requireNonEmpty(funcName);

    _basisCol.memReAlloc(cSize3(cBasis<T>::_nmax+1,
                         cBasis<T>::_coordCol.getNrow(),
                         1));

    switch (accuracy) {
        case BASIS_ACCU_APPROX:
        case BASIS_ACCU_ACCURATE:
            basis_cart_gaussian(cBasis<T>::_coordCol.getAddrData(), cBasis<T>::_coordCol.getNrow(),
                                           _basisCol.getAddrData(),
                                cBasis<T>::_nmax, sigma,
                                accuracy);
            break;
        case BASIS_ACCU_MP:
            floatmp     *arrayCoord, *arrayBasis;

            arrayCoord = new floatmp [cBasis<T>::_coordCol.getNrow()];
            arrayBasis = new floatmp [cBasis<T>::_coordCol.getNrow()*(cBasis<T>::_nmax+1)];

            for (size_t i = 0; i < cBasis<T>::_coordCol.getNrow(); i++) {
                arrayCoord[i] = cBasis<T>::_coordCol[i];
            }

            basis_cart_gaussian(arrayCoord, cBasis<T>::_coordCol.getNrow(),
                                arrayBasis,
                                cBasis<T>::_nmax, (floatmp) sigma);

            for (size_t i = 0; i < cBasis<T>::_coordCol.getNrow()*(cBasis<T>::_nmax+1); i++) {
                _basisCol[i] = arrayBasis[i].convert_to<T>();
            }

            delete [] arrayCoord;
            delete [] arrayBasis;

            break;
        default:
            ERROR(funcName, "unsupported accuracy");
    }
}

template <typename T>
void cBasisGaussian<T>::genBasisSec(T sigma, eBasisAccuracy accuracy)
{
    const std::string   funcName("void cBasisGaussian<T>::genBasisSec("
                                    "T sigma, eBasisAccuracy accuracy)");

    cBasis<T>::_coordSec.requireNonEmpty(funcName);

    _basisSec.memReAlloc(cSize3(cBasis<T>::_nmax+1,
                         cBasis<T>::_coordSec.getNrow(),
                         1));

    switch (accuracy) {
        case BASIS_ACCU_APPROX:
        case BASIS_ACCU_ACCURATE:
            basis_cart_gaussian(cBasis<T>::_coordSec.getAddrData(), cBasis<T>::_coordSec.getNrow(),
                                           _basisSec.getAddrData(),
                                cBasis<T>::_nmax, sigma,
                                accuracy);
            break;
        case BASIS_ACCU_MP:
            floatmp     *arrayCoord, *arrayBasis;

            arrayCoord = new floatmp [cBasis<T>::_coordSec.getNrow()];
            arrayBasis = new floatmp [cBasis<T>::_coordSec.getNrow()*(cBasis<T>::_nmax+1)];

            for (size_t i = 0; i < cBasis<T>::_coordSec.getNrow(); i++) {
                arrayCoord[i] = cBasis<T>::_coordSec[i];
            }

            basis_cart_gaussian(arrayCoord, cBasis<T>::_coordSec.getNrow(),
                                arrayBasis,
                                cBasis<T>::_nmax, (floatmp) sigma);

            for (size_t i = 0; i < cBasis<T>::_coordSec.getNrow()*(cBasis<T>::_nmax+1); i++) {
                _basisSec[i] = arrayBasis[i].convert_to<T>();
            }

            delete [] arrayCoord;
            delete [] arrayBasis;

            break;
        default:
            ERROR(funcName, "unsupported accuracy");
    }
}

template <typename T>
void cBasisGaussian<T>::genBasis2D(const cInt3& order, cData3<T>& basis2D)
{
    const std::string   funcName("void cBasisGaussian<T>::genBasis2D("
                                    "const cInt3& order, cData3& basis2D)");

    require(order[0] >= 0, funcName + ": order[0] is negative");
    require(order[1] >= 0, funcName + ": order[1] is negative");
    require((unsigned int) order[0] < _basisRow.getNrow(), funcName + ": order[0] is too high");
    require((unsigned int) order[1] < _basisCol.getNrow(), funcName + ": order[1] is too high");
    require(order[2] == 0, funcName + ": order[2] != 0");

    basis2D.memReAlloc(cSize3(_basisRow.getNcol(),
                              _basisCol.getNcol(),
                              1));

    size_t      offsetRow = order[0] * _basisRow.getNcol();
    size_t      offsetCol = order[1] * _basisCol.getNcol();

    array_linalg_prodvec2mat(  basis2D.getAddrData(),
                             _basisRow.getAddrData()+offsetRow, _basisRow.getNcol(),
                             _basisCol.getAddrData()+offsetCol, _basisCol.getNcol());
}

template <typename T>
void cBasisGaussian<T>::genBasis3D(const cInt3& order, cData3<T>& basis2D)
{
    const std::string   funcName("void cBasisGaussian<T>::genBasis3D("
                                    "const cInt3& order, cData3& basis2D)");

    require(order[0] >= 0, funcName + ": order[0] is negative");
    require(order[1] >= 0, funcName + ": order[1] is negative");
    require(order[2] >= 0, funcName + ": order[2] is negative");
    require((unsigned int) order[0] < _basisRow.getNrow(), funcName + ": order[0] is too high");
    require((unsigned int) order[1] < _basisCol.getNrow(), funcName + ": order[1] is too high");
    require((unsigned int) order[2] < _basisSec.getNrow(), funcName + ": order[2] is too high");

    basis2D.memReAlloc(cSize3(_basisRow.getNcol(),
                              _basisCol.getNcol(),
                              _basisSec.getNcol()));

    size_t      offsetRow = order[0] * _basisRow.getNcol();
    size_t      offsetCol = order[1] * _basisCol.getNcol();
    size_t      offsetSec = order[2] * _basisSec.getNcol();

    array_linalg_prodvec2mat(  basis2D.getAddrData(),
                             _basisRow.getAddrData()+offsetRow, _basisRow.getNcol(),
                             _basisCol.getAddrData()+offsetCol, _basisCol.getNcol(),
                             _basisSec.getAddrData()+offsetSec, _basisSec.getNcol());
}

// instantiation
template class cBasisGaussian<float >;
template class cBasisGaussian<double>;

} // namespace gem
