/***********************************************************************
 *  File:       cCoordinate.hpp
 *
 *  Purpose:    Header file for a coordinate system class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CCOORDINATE_HPP__
#define __GEM_CCOORDINATE_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cCoordinate
{
private:

public:
     cCoordinate(void) {};
    ~cCoordinate()     {};

    // 2D
    void computeAngle          (const cData3<T>& X, const cData3<T>& Y, cData3<T>& theta);
    void computeDistance       (const cData3<T>& X, const cData3<T>& Y, cData3<T>& dist );
    void computeDistanceSquared(const cData3<T>& X, const cData3<T>& Y, cData3<T>& dist2);

    // 3D
    void computeAngle          (const cData3<T>& X, const cData3<T>& Y, const cData3<T>& Z, cData3<T>& angle);
    void computeDistance       (const cData3<T>& X, const cData3<T>& Y, const cData3<T>& Z, cData3<T>& dist );
    void computeDistanceSquared(const cData3<T>& X, const cData3<T>& Y, const cData3<T>& Z, cData3<T>& dist2);
};

} // namespace gem

#endif
