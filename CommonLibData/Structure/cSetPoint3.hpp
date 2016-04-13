/***********************************************************************
 *  File:       cSetPoint3.hpp
 *
 *  Purpose:    Header file for a 3D point set class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CSETPOINT3_HPP__
#define __GEM_CSETPOINT3_HPP__

#include "macro.hpp"
#include "cVector3.hpp"

#include <vector>

namespace gem {

enum eFlipPlane
{
    FLIP_PLANE_YZ = 0,
    FLIP_PLANE_ZX = 1,
    FLIP_PLANE_XY = 2
};

template <typename T>
class cSetPoint3 : public std::vector< cVector3<T> >
{
protected:

public:
    cSetPoint3(void) {};
    cSetPoint3(const cSetPoint3& other);
    ~cSetPoint3()    {};

    // std::vector wrap-ups
    void                push_back(const cVector3<T>& other) {        std::vector< cVector3<T> >::push_back(other);                             }
    void                clear    (void)                     {        std::vector< cVector3<T> >::clear();                                      }
    void                erase    (size_t i)                 {        std::vector< cVector3<T> >::erase(std::vector< cVector3<T> >::begin()+i); }
    void                resize   (size_t n)                 {        std::vector< cVector3<T> >::resize(n);                                    }
    size_t              size     (void)     const           { return std::vector< cVector3<T> >::size();                                       }
    const cVector3<T>   at       (size_t i) const           { return std::vector< cVector3<T> >::at(i);                                        }
          cVector3<T>&  at       (size_t i)                 { return std::vector< cVector3<T> >::at(i);                                        }

    // statistics
    void                requireNonEmpty(void) const { require(size() > 0, "point set is empty"); }
    const cVector3<T>   getMean(void) const;

    // operator overloading
          cVector3<T>&  operator[](size_t i)       { assert(i < size());  return at(i); }
    const cVector3<T>   operator[](size_t i) const { assert(i < size());  return at(i); }

    cSetPoint3&         operator= (const cSetPoint3& other);

    template <typename T1>
    friend std::ostream& operator<<(std::ostream& os, const cSetPoint3<T1>& sp3);

    // transform
    void    opFlip     (eFlipPlane plane = FLIP_PLANE_YZ);

    void    opRotate   (const cVector3<T>& angle );
    void    opScale    (const cVector3<T>& factor);
    void    opTranslate(const cVector3<T>& offset);

    // comparison
    T       computeSD  (const cSetPoint3& other) const;
    T       computeRMSD(const cSetPoint3& other) const;
};

} // namespace gem

#endif
