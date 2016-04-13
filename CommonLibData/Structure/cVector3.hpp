/***********************************************************************
 *  File:       cVector3.hpp
 *
 *  Purpose:    Header file for a 3D vector class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CVECTOR3_HPP__
#define __GEM_CVECTOR3_HPP__

#include <stdint.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace gem {

template <typename T>
class cVector3
{
protected:
    T       _v[3];

public:
    cVector3(T v1 = 0, T v2 = 0, T v3 = 0);
    cVector3(T v[3]);
    cVector3(const cVector3& v);
    ~cVector3() {};

    // operator overloading
    T&                   operator[](size_t i)       { assert(i < 3);  return _v[i]; }
    T                    operator[](size_t i) const { assert(i < 3);  return _v[i]; }

    cVector3&            operator= (T value);
    template <typename U>
    cVector3&            operator= (const cVector3<U>& other);

    const cVector3       operator+ (T value)           const;
          cVector3&      operator+=(T value);
    const cVector3       operator+ (const cVector3& v) const;
          cVector3&      operator+=(const cVector3& v);
    const cVector3       operator- (T value)           const;
          cVector3&      operator-=(T value);
    const cVector3       operator- (const cVector3& v) const;
          cVector3&      operator-=(const cVector3& v);
    const cVector3       operator* (T value)           const;
          cVector3&      operator*=(T value);
    const cVector3       operator* (const cVector3& v) const;
          cVector3&      operator*=(const cVector3& v);
    const cVector3       operator/ (T value)           const;
          cVector3&      operator/=(T value);
    const cVector3       operator/ (const cVector3& v) const;
          cVector3&      operator/=(const cVector3& v);
    const cVector3       operator- (void)              const;
    bool                 operator==(const cVector3& v) const;
    bool                 operator!=(const cVector3& v) const;
    bool                 operator< (const cVector3& v) const;
    bool                 operator<=(const cVector3& v) const;
    bool                 operator> (const cVector3& v) const;
    bool                 operator>=(const cVector3& v) const;
    T                    operator| (const cVector3& v) const;
    const cVector3       operator^ (const cVector3& v) const;

    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const cVector3<U>& v);

    // property
    T       getMax    (void) const;
    T       getMin    (void) const;
    T       getSum    (void) const;
    T       getSum2   (void) const;
    T       getLength (void) const;
    T       getProduct(void) const;

    // type casting
    template <typename U>
    void    opTypeCast(const cVector3<U>& other);

    // indexing
    void    ind2sub(const cVector3<T>& size, T idx);
    T       sub2ind(const cVector3<T>& size) const;
};

template <typename T>
cVector3<T> operator*(T value, const cVector3<T>& v);

typedef     cVector3<unsigned int>  cUInt3;
typedef     cVector3<  signed int>  cInt3;
typedef     cVector3<size_t   >     cSize3;
typedef     cVector3<ptrdiff_t>     cOffset3;
typedef     cVector3<float    >     cPoint3f;
typedef     cVector3<double   >     cPoint3d;

// ---------------------------------------------------------------------
// implementation for templated member function
// ---------------------------------------------------------------------

template <typename T> template <typename U>
cVector3<T>& cVector3<T>::operator=(const cVector3<U>& other)
{
    _v[0] = (T) other[0];
    _v[1] = (T) other[1];
    _v[2] = (T) other[2];

    return *this;
}

template <typename T> template <typename U>
void cVector3<T>::opTypeCast(const cVector3<U>& other)
{
    _v[0] = (T) other[0];
    _v[1] = (T) other[1];
    _v[2] = (T) other[2];
}

} // namespace gem

#endif
