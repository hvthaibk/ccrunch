/***********************************************************************
 *  File:       cNDgrid.hpp
 *
 *  Purpose:    Header file for a ndgrid generation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CNDGRID_HPP__
#define __GEM_CNDGRID_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cNDgrid
{
private:

public:
     cNDgrid(void) {};
    ~cNDgrid()    {};

    // 2D R-C
    void genRC(const cData3<T>& rVec, const cData3<T>& cVec,
                     cData3<T>& R,          cData3<T>& C);

    // 2D X-Y
    void genXY(const cData3<T>& xVec, const cData3<T>& yVec,
                     cData3<T>& X,          cData3<T>& Y);

    // 3D R-C-S
    void genRCS(const cData3<T>& rVec, const cData3<T>& cVec, const cData3<T>& sVec,
                      cData3<T>& R,          cData3<T>& C,          cData3<T>& S);

    // 3D X-Y-Z
    void genXYZ(const cData3<T>& xVec, const cData3<T>& yVec, const cData3<T>& zVec,
                      cData3<T>& X,          cData3<T>& Y,          cData3<T>& Z);

    // 2D R-C for FFTW r2c_2d
    void genRCfreq(size_t nrow,  size_t ncol,
                   cData3<T>& R, cData3<T>& C);
};

} // namespace gem

#endif
