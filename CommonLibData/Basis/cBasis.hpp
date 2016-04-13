/***********************************************************************
 *  File:       cBasis.hpp
 *
 *  Purpose:    Header file for a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CBASIS_HPP__
#define __GEM_CBASIS_HPP__

#include "basis.hpp"
#include "cData3.hpp"

#include <vector>

namespace gem {

template <typename T>
class cBasis
{
protected:
    unsigned int            _nmax;
    std::vector<cInt3>      _order;

    unsigned int            _dim;
    cVector3<T>             _origin;
    cVector3<T>             _spacing;
    cSize3                  _size;

    // Cartesian
    cData3<T>               _coordRow;
    cData3<T>               _coordCol;
    cData3<T>               _coordSec;

    cData3<T>               _coordX;
    cData3<T>               _coordY;
    cData3<T>               _coordZ;

    // polar / spherical
    cData3<T>               _coordRad;
    cData3<T>               _coordThe;
    cData3<T>               _coordPhi;

public:
     cBasis(void);
    ~cBasis()  {};

    // debug
    void    printInfo          (void);
    void    printOrder         (void);
    void    printCoordXY       (void);
    void    printCoordXYZ      (void);
    void    printCoordCartesian(void);
    void    printCoordPolar    (void);
    void    printCoordSpherical(void);

    // generate order set
    void    genOrderFull2DLozenge   (unsigned int nmax);
    void    genOrderFull2DSquare    (unsigned int nmax);
    void    genOrderFull3DCube      (unsigned int nmax);
    void    genOrderFull3DOctahedron(unsigned int nmax);

    void    genOrderInc2DLozenge    (unsigned int n);
    void    genOrderInc2DSquare     (unsigned int n);
    void    genOrderInc3DCube       (unsigned int n);
    void    genOrderInc3DOctahedron (unsigned int n);

    // generate and covert coordinates
    void    genCoordCartesian(unsigned int       dim,
                              const cVector3<T>& origin,
                              const cVector3<T>& spacing,
                              const cSize3&      size);

    void    convertCoordCart2Pol(void);
    void    convertCoordCart2Sph(void);

    // set/get params
    void                        setNmax (unsigned int nmax) { _nmax = nmax;  }
    const std::vector<cInt3>&   getOrder(void)              { return _order; }
};

} // namespace gem

#endif
