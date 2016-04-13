/***********************************************************************
 *  File:       cBasis.cpp
 *
 *  Purpose:    Implementation of a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cBasis.hpp"

#include "cNDgrid.hpp"

namespace gem {

template <typename T>
cBasis<T>::cBasis(void)
{
    _nmax = 0;
    _order.clear();

    _dim = 0;
}

template <typename T>
void cBasis<T>::printInfo(void)
{
        std::cout << "Basis common info: "         << std::endl;
        std::cout << "     nmax    = " << _nmax    << std::endl;
        std::cout << "     dim     = " << _dim     << std::endl;
        std::cout << "     origin  = " << _origin  << std::endl;
        std::cout << "     spacing = " << _spacing << std::endl;
}

template <typename T>
void cBasis<T>::printOrder(void)
{
    const std::string   funcName("void cBasis<T>::printOrder(void)");

    require(_order.size() > 0, funcName + ": order set is empty");

    std::cout << "#order = " << _order.size() << std::endl;
    for (size_t i = 0; i < _order.size(); i++) {
        std::cout << _order[i] << std::endl;
    }
    std::cout << std::endl;
}

template <typename T>
void cBasis<T>::printCoordXY(void)
{
    const std::string   funcName("void cBasis<T>::printCoordXY(void)");

    _coordX.requireNonEmpty(funcName);
    _coordY.requireNonEmpty(funcName);

    cData3<T>       X, Y;

    cNDgrid<T>().genXY(_coordX, _coordY, X, Y);

    X.printData("X");
    Y.printData("Y");
}

template <typename T>
void cBasis<T>::printCoordXYZ(void)
{
    const std::string   funcName("void cBasis<T>::printCoordXYZ(void)");

    _coordX.requireNonEmpty(funcName);
    _coordY.requireNonEmpty(funcName);
    _coordZ.requireNonEmpty(funcName);

    cData3<T>       X, Y, Z;

    cNDgrid<T>().genXYZ(_coordX, _coordY, _coordZ, X, Y, Z);

    X.printData("X");
    Y.printData("Y");
    Z.printData("Z");
}

template <typename T>
void cBasis<T>::printCoordCartesian(void)
{
    const std::string   funcName("void cBasis<T>::printCoordCartesian(void)");

    _coordRow.requireNonEmpty(funcName);

                                     _coordRow.printData("coordRow");
    if (_coordCol.getNelement() > 0) _coordCol.printData("coordCol");
    if (_coordSec.getNelement() > 0) _coordSec.printData("coordSec");

    _coordX.requireNonEmpty(funcName);

                                   _coordX.printData("coordX");
    if (_coordY.getNelement() > 0) _coordY.printData("coordY");
    if (_coordZ.getNelement() > 0) _coordZ.printData("coordZ");
}

template <typename T>
void cBasis<T>::printCoordPolar(void)
{
    const std::string   funcName("void cBasis<T>::printCoordPolar(void)");

    _coordRad.requireNonEmpty(funcName);
    _coordThe.requireNonEmpty(funcName);

    _coordRad.printData("coordRad");
    _coordThe.printData("coordThe");
}

template <typename T>
void cBasis<T>::printCoordSpherical(void)
{
    const std::string   funcName("void cBasis<T>::printCoordSpherical(void)");

    _coordRad.requireNonEmpty(funcName);
    _coordThe.requireNonEmpty(funcName);
    _coordPhi.requireNonEmpty(funcName);

    _coordRad.printData("coordRad");
    _coordThe.printData("coordThe");
    _coordPhi.printData("coordPhi");
}

template <typename T>
void cBasis<T>::genOrderFull2DLozenge(unsigned int nmax)
{
    _nmax = nmax;
    _order.clear();
    _order.resize((_nmax+1)*(_nmax+2)/2);

    for (unsigned int n = 0, i = 0; n <= _nmax; n++) {
        for (unsigned int m = 0; m <= _nmax-n; m++, i++) {

            _order[i] = cInt3(n,m,0);
        }
    }
}

template <typename T>
void cBasis<T>::genOrderFull2DSquare(unsigned int nmax)
{
    _nmax = nmax;
    _order.clear();
    _order.resize(pow2(_nmax+1));

    for (unsigned int n = 0, i = 0; n <= _nmax; n++) {
        for (unsigned int m = 0; m <= _nmax; m++, i++) {

            _order[i] = cInt3(n,m,0);
        }
    }
}

template <typename T>
void cBasis<T>::genOrderFull3DCube(unsigned int nmax)
{
    _nmax = nmax;
    _order.clear();
    _order.resize(pow3(_nmax+1));

    for (unsigned int n = 0, i = 0; n <= _nmax; n++) {
        for (unsigned int l = 0; l <= _nmax; l++) {
            for (unsigned int m = 0; m <= _nmax; m++, i++) {

                _order[i] = cInt3(n,l,m);
            }
        }
    }
}

template <typename T>
void cBasis<T>::genOrderFull3DOctahedron(unsigned int nmax)
{
    _nmax = nmax;
    _order.clear();
    _order.resize((_nmax+1)*(_nmax+2)*(_nmax+3)/6);

    for (unsigned int n = 0, i = 0; n <= _nmax; n++) {
        for (unsigned int l = 0; l <= _nmax-n; l++) {
            for (unsigned int m = 0; m <= _nmax-n-l; m++, i++) {

                _order[i] = cInt3(n,l,m);
            }
        }
    }
}

template <typename T>
void cBasis<T>::genOrderInc2DLozenge(unsigned int n)
{
    unsigned int            i(0);

    _order.resize(n+1);

    for (unsigned int m = 0; m <= n; m++) {

        _order[i++] = cInt3(m,n-m,0);
    }
}

template <typename T>
void cBasis<T>::genOrderInc2DSquare(unsigned int n)
{
    unsigned int            i(0);

    _order.resize(pow2(n+1)-pow2(n));

    _order[i++] = cInt3(n,n,0);

    for (unsigned int m = 0; m < n; m++) {

        _order[i++] = cInt3(n,m,0);
        _order[i++] = cInt3(m,n,0);
    }
}

template <typename T>
void cBasis<T>::genOrderInc3DCube(unsigned int n)
{
    unsigned int            i(0);

    _order.resize(pow3(n+1)-pow3(n));

    _order[i++] = cInt3(n,n,n);

    for (unsigned int l = 0; l < n; l++) {

        _order[i++] = cInt3(n,n,l);
        _order[i++] = cInt3(n,l,n);
        _order[i++] = cInt3(l,n,n);

        _order[i++] = cInt3(n,l,l);
        _order[i++] = cInt3(l,n,l);
        _order[i++] = cInt3(l,l,n);

        for (unsigned int m = 0; m < l; m++) {

            _order[i++] = cInt3(n,l,m);
            _order[i++] = cInt3(n,m,l);
            _order[i++] = cInt3(l,n,m);
            _order[i++] = cInt3(m,n,l);
            _order[i++] = cInt3(l,m,n);
            _order[i++] = cInt3(m,l,n);
        }
    }
}

template <typename T>
void cBasis<T>::genOrderInc3DOctahedron(unsigned int n)
{
    unsigned int            i(0);

    _order.resize((n+1)*(n+2)/2);

    for (unsigned int m = 0; m <= n; m++) {
        for (unsigned int l = 0; l <= n-m; l++) {

            _order[i++] = cInt3(m,l,n-m-l);
        }
    }
}

template <typename T>
void cBasis<T>::genCoordCartesian(unsigned int       dim,
                                  const cVector3<T>& origin,
                                  const cVector3<T>& spacing,
                                  const cSize3&      size)
{
    const std::string   funcName("void cBasis<T>::genCoordCartesian("
                                    "unsigned int dim, "
                                    "const cVector3<T>& origin, "
                                    "const cVector3<T>& spacing, "
                                    "const cSize3& size)");

    require(size.getProduct() > 1, funcName + ": invalid size");

    _dim     = dim;
    _origin  = origin;
    _spacing = spacing;
    _size    = size;

    switch (_dim) {
        case 3:
            break;
        case 2:
            require(_size[2] == 1, funcName + ": invalid size[2] for 2D");
            break;
        case 1:
            require(_size[2] == 1, funcName + ": invalid size[2] for 1D");
            require(_size[1] == 1, funcName + ": invalid size[1] for 1D");
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    switch (_dim) {
        case 3:
            require(_spacing[2] > 0, funcName + ": invalid spacing[2]");
            require(   _size[2] > 0, funcName + ": invalid size[2]");

            _coordSec.linspace(_origin[2], _origin[2]+((T)_size[2]-1)*_spacing[2], _size[2]);
        case 2:
            require(_spacing[1] > 0, funcName + ": invalid spacing[1]");
            require(   _size[1] > 0, funcName + ": invalid size[1]");

            _coordCol.linspace(_origin[1], _origin[1]+((T)_size[1]-1)*_spacing[1], _size[1]);
        case 1:
            require(_spacing[0] > 0, funcName + ": invalid spacing[0]");
            require(   _size[0] > 0, funcName + ": invalid size[0]");

            _coordRow.linspace(_origin[0], _origin[0]+((T)_size[0]-1)*_spacing[0], _size[0]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    switch (_dim) {
        case 1:
            _coordX = _coordRow;
            break;
        case 2:
            _coordX = _coordCol;
            _coordY = _coordRow;    _coordY.opMathInv();
            break;
        case 3:
            _coordX = _coordCol;
            _coordY = _coordRow;    _coordY.opMathInv();
            _coordZ = _coordSec;
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cBasis<T>::convertCoordCart2Pol(void)
{
    const std::string   funcName("void cBasis<T>::convertCoordCart2Pol(void)");

    cData3<T>       X, Y;

    cNDgrid<T>().genXY(_coordX, _coordY, X, Y);

    _coordRad.memReAlloc(cSize3(X.getNrow(),X.getNcol(),1));
    _coordThe.memReAlloc(cSize3(X.getNrow(),X.getNcol(),1));

    coord_cart2pol(        X.getAddrData(),         Y.getAddrData(),
                   _coordRad.getAddrData(), _coordThe.getAddrData(),
                   _coordRad.getNelement());
}

template <typename T>
void cBasis<T>::convertCoordCart2Sph(void)
{
    const std::string   funcName("void cBasis<T>::convertCoordCart2Sph(void)");

    cData3<T>       X, Y, Z;

    cNDgrid<T>().genXYZ(_coordX, _coordY, _coordZ, X, Y, Z);

    _coordRad.memReAlloc(cSize3(X.getNrow(),X.getNcol(),X.getNsec()));
    _coordThe.memReAlloc(cSize3(X.getNrow(),X.getNcol(),X.getNsec()));
    _coordPhi.memReAlloc(cSize3(X.getNrow(),X.getNcol(),X.getNsec()));

    coord_cart2sph(        X.getAddrData(),         Y.getAddrData(),         Z.getAddrData(),
                   _coordRad.getAddrData(), _coordThe.getAddrData(), _coordPhi.getAddrData(),
                   _coordRad.getNelement());
}

// instantiation
template class cBasis<float >;
template class cBasis<double>;

} // namespace gem
