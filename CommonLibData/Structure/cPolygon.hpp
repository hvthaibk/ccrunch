/***********************************************************************
 *  File:       cPolygon.hpp
 *
 *  Purpose:    Header file for a polygon class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CPOLYGON_HPP__
#define __GEM_CPOLYGON_HPP__

#include "cSetPoint3.hpp"
#include "cSegment3.hpp"
#include "cData3.hpp"
#include "vector"
#include "algorithm"

namespace gem {

template <typename T>
class cPolygon : public cSetPoint3<T>
{
protected:

public:
    cPolygon(void) {};
    cPolygon(const cPolygon& other);
    ~cPolygon()    {};

    // property
    bool    isPolygon(void) const;
    bool    isSimple (void) const;
    void    opFill(const cVector3<size_t>& imSize,cData3<T>& imData) const;

    // coordenates
    T       getCoordMinX(void) const;
    T       getCoordMaxX(void) const;
    T       getCoordMinY(void) const;
    T       getCoordMaxY(void) const;

    // segments
    cSegment3<T>    *getSegments(void) const;
};

} // namespace gem

#endif
