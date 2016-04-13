/***********************************************************************
 *  File:       cNDgrid.cpp
 *
 *  Purpose:    Implementation of a ndgrid generation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cNDgrid.hpp"

namespace gem {

template <typename T>
void cNDgrid<T>::genRC(const cData3<T>& rVec, const cData3<T>& cVec,
                             cData3<T>& R,          cData3<T>& C)
{
    const std::string   funcName("void cNDgrid<T>::genRC("
                                    "const cData3<T>& rVec, "
                                    "const cData3<T>& cVec, "
                                    "cData3<T>& R, "
                                    "cData3<T>& C)");

    require(rVec.getNrow() > 0 && rVec.getNcol() == 1 && rVec.getNsec() == 1,
            funcName + ": rVec is not a row vector");
    require(cVec.getNrow() > 0 && cVec.getNcol() == 1 && cVec.getNsec() == 1,
            funcName + ": cVec is not a row vector");

    size_t      nr = rVec.getNrow();
    size_t      nc = cVec.getNrow();

    R.memReAlloc(cSize3(nr,nc,1));
    C.memReAlloc(cSize3(nr,nc,1));

    for (size_t ir = 0; ir < nr; ir ++) {
        R.opReplace(rVec[ir],cSize3(1,nc,1),cSize3(ir,0,0));
    }

    for (size_t ic = 0; ic < nc; ic ++) {
        C.opReplace(cVec[ic],cSize3(nr,1,1),cSize3(0,ic,0));
    }
}

template <typename T>
void cNDgrid<T>::genXY(const cData3<T>& xVec, const cData3<T>& yVec,
                             cData3<T>& X,          cData3<T>& Y)
{
    const std::string   funcName("void cNDgrid<T>::genXY("
                                    "const cData3<T>& xVec, "
                                    "const cData3<T>& yVec, "
                                    "cData3<T>& X, "
                                    "cData3<T>& Y)");

    require(xVec.getNrow() > 0 && xVec.getNcol() == 1 && xVec.getNsec() == 1,
            funcName + ": xVec is not a row vector");
    require(yVec.getNrow() > 0 && yVec.getNcol() == 1 && yVec.getNsec() == 1,
            funcName + ": yVec is not a row vector");

    size_t      nx = xVec.getNrow();
    size_t      ny = yVec.getNrow();

    X.memReAlloc(cSize3(ny,nx,1));
    Y.memReAlloc(cSize3(ny,nx,1));

    for (size_t iy = 0; iy < ny; iy ++) {
        Y.opReplace(yVec[iy],cSize3(1,nx,1),cSize3(iy,0,0));
    }

    for (size_t ix = 0; ix < nx; ix ++) {
        X.opReplace(xVec[ix],cSize3(ny,1,1),cSize3(0,ix,0));
    }
}

template <typename T>
void cNDgrid<T>::genRCS(const cData3<T>& rVec, const cData3<T>& cVec, const cData3<T>& sVec,
                              cData3<T>& R,          cData3<T>& C,          cData3<T>& S)
{
    const std::string   funcName("void cNDgrid<T>::genRCS("
                                    "const cData3<T>& rVec, "
                                    "const cData3<T>& cVec, "
                                    "const cData3<T>& sVec, "
                                    "cData3<T>& R, "
                                    "cData3<T>& C, "
                                    "cData3<T>& S)");

    require(rVec.getNrow() > 0 && rVec.getNcol() == 1 && rVec.getNsec() == 1,
            funcName + ": rVec is not a row vector");
    require(cVec.getNrow() > 0 && cVec.getNcol() == 1 && cVec.getNsec() == 1,
            funcName + ": cVec is not a row vector");
    require(sVec.getNrow() > 0 && sVec.getNcol() == 1 && sVec.getNsec() == 1,
            funcName + ": sVec is not a row vector");

    size_t      nr = rVec.getNrow();
    size_t      nc = cVec.getNrow();
    size_t      ns = sVec.getNrow();

    R.memReAlloc(cSize3(nr,nc,ns));
    C.memReAlloc(cSize3(nr,nc,ns));
    S.memReAlloc(cSize3(nr,nc,ns));

    for (size_t ir = 0; ir < nr; ir ++) {
        R.opReplace(rVec[ir],cSize3(1,nc,ns),cSize3(ir,0,0));
    }

    for (size_t ic = 0; ic < nc; ic ++) {
        C.opReplace(cVec[ic],cSize3(nr,1,ns),cSize3(0,ic,0));
    }

    for (size_t is = 0; is < ns; is ++) {
        S.opReplace(sVec[is],cSize3(nr,nc,1),cSize3(0,0,is));
    }
}

template <typename T>
void cNDgrid<T>::genXYZ(const cData3<T>& xVec, const cData3<T>& yVec, const cData3<T>& zVec,
                              cData3<T>& X,          cData3<T>& Y,          cData3<T>& Z)
{
    const std::string   funcName("void cNDgrid<T>::genXYZ("
                                    "const cData3<T>& xVec, "
                                    "const cData3<T>& yVec, "
                                    "const cData3<T>& zVec, "
                                    "cData3<T>& X, "
                                    "cData3<T>& Y, "
                                    "cData3<T>& Z)");

    require(xVec.getNrow() > 0 && xVec.getNcol() == 1 && xVec.getNsec() == 1,
            funcName + ": xVec is not a row vector");
    require(yVec.getNrow() > 0 && yVec.getNcol() == 1 && yVec.getNsec() == 1,
            funcName + ": yVec is not a row vector");
    require(zVec.getNrow() > 0 && zVec.getNcol() == 1 && zVec.getNsec() == 1,
            funcName + ": zVec is not a row vector");

    size_t      nx = xVec.getNrow();
    size_t      ny = yVec.getNrow();
    size_t      nz = zVec.getNrow();

    X.memReAlloc(cSize3(ny,nx,nz));
    Y.memReAlloc(cSize3(ny,nx,nz));
    Z.memReAlloc(cSize3(ny,nx,nz));

    for (size_t iy = 0; iy < ny; iy ++) {
        Y.opReplace(yVec[iy],cSize3(1,nx,nz),cSize3(iy,0,0));
    }

    for (size_t ix = 0; ix < nx; ix ++) {
        X.opReplace(xVec[ix],cSize3(ny,1,nz),cSize3(0,ix,0));
    }

    for (size_t iz = 0; iz < nz; iz ++) {
        Z.opReplace(zVec[iz],cSize3(ny,nx,1),cSize3(0,0,iz));
    }
}

template <typename T>
void cNDgrid<T>::genRCfreq(size_t nrow,  size_t ncol,
                           cData3<T>& R, cData3<T>& C)
{
    const std::string   funcName("void cNDgrid<T>::genRCfreq("
                                    "size_t nrow, size_t ncol, "
                                    "cData3<T>& R, cData3<T>& C)");

    cData3<T>   rVec, cVec;

    cVec.linspace(0, (T) (ncol/2), ncol/2+1);
    if (ncol % 2 == 0) {
        cVec[ncol/2] = - (T) (ncol/2);
    }

    rVec.linspace(0, (T) (nrow-1), nrow);
    for (size_t i = (nrow+1)/2; i < nrow; i++) {
        rVec[i] = rVec[i] - (T) nrow;
    }

    genRC(rVec, cVec, R, C);
}

// instantiation
template class cNDgrid<float >;
template class cNDgrid<double>;

} // namespace gem
