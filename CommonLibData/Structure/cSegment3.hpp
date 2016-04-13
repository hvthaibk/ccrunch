/***********************************************************************
 *  File:       cSegment3.hpp
 *
 *  Purpose:    Header file for a 3D segment class
 *
 *  Author:     Lucia Martin Reixach, Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Lucia Martin Reixach, Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSEGMENT3_HPP__
#define __GEM_CSEGMENT3_HPP__

#include "cVector3.hpp"

namespace gem {

template <typename T>
class cSegment3
{
private:
    cVector3<T>     _p[2];

public:
    cSegment3(void);
    cSegment3(const cVector3<T>& p1, const cVector3<T>& p2);
    cSegment3(const cSegment3& other);
    ~cSegment3() {};

    // operator overloading
    cVector3<T>&         operator[](size_t i)       { assert(i < 2);  return _p[i]; }
    cVector3<T>          operator[](size_t i) const { assert(i < 2);  return _p[i]; }

    cSegment3&           operator= (const cSegment3& other);

    template <typename T1>
    friend std::ostream& operator<<(std::ostream& os, const cSegment3<T1>& segment);

    // property
    T       getAngle2D (void) const;
    T       getLength2D(void) const;
    T       getLength3D(void) const;

    // intersection
    bool    isCoincide (const cSegment3& other) const;
    bool    isIntersect(const cSegment3& other) const;
    bool    isParallel (const cSegment3& other) const;
    bool    isTouched  (const cSegment3& other) const;

    cVector3<T>     getIntersectPoint(const cSegment3& other) const;
};

} // namespace gem

#endif
