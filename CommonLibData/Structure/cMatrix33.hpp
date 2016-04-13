/***********************************************************************
 *  File:       cMatrix33.hpp
 *
 *  Purpose:    Header file for a 3x3 matrix class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CMATRIX33_HPP__
#define __GEM_CMATRIX33_HPP__

#include "cVector3.hpp"

namespace gem {

template <typename T>
class cMatrix33
{
private:
    cVector3<T>     _m[3];

public:
    cMatrix33(T m11 = 0, T m12 = 0, T m13 = 0,
              T m21 = 0, T m22 = 0, T m23 = 0,
              T m31 = 0, T m32 = 0, T m33 = 0);
    cMatrix33(const T v[9]);
    cMatrix33(const cVector3<T>& v1, const cVector3<T>& v2, const cVector3<T>& v3);
    cMatrix33(const cVector3<T> v[3]);
    cMatrix33(const cMatrix33& m);
    ~cMatrix33() {};

    T    getLength(void) const;
    void transpose(void);

    cVector3<T>&        operator[](size_t i)       { assert(i < 3);  return _m[i]; }
    cVector3<T>         operator[](size_t i) const { assert(i < 3);  return _m[i]; }

    cMatrix33&          operator= (T value);
    cMatrix33&          operator= (const cMatrix33& m);

    const cMatrix33     operator+ (T value)            const;
          cMatrix33&    operator+=(T value);
    const cMatrix33     operator+ (const cMatrix33& m) const;
          cMatrix33&    operator+=(const cMatrix33& m);
    const cMatrix33     operator- (T value)            const;
          cMatrix33&    operator-=(T value);
    const cMatrix33     operator- (const cMatrix33& m) const;
          cMatrix33&    operator-=(const cMatrix33& m);
    const cMatrix33     operator* (T value)            const;
          cMatrix33&    operator*=(T value);
    const cMatrix33     operator* (const cMatrix33& m) const;
          cMatrix33&    operator*=(const cMatrix33& m);
    const cMatrix33     operator/ (T value)            const;
          cMatrix33&    operator/=(T value);
    const cMatrix33     operator/ (const cMatrix33& m) const;
          cMatrix33&    operator/=(const cMatrix33& m);
    const cMatrix33     operator- (void)               const;
    bool                operator==(const cMatrix33& m) const;
    bool                operator!=(const cMatrix33& m) const;
    bool                operator< (const cMatrix33& m) const;
    bool                operator<=(const cMatrix33& m) const;
    bool                operator> (const cMatrix33& m) const;
    bool                operator>=(const cMatrix33& m) const;
    T                   operator| (const cMatrix33& m) const;

    template <typename T1>
    friend std::ostream& operator<<(std::ostream& os, const cMatrix33<T1>& m);
};

template <typename T>
cMatrix33<T> operator*(T d, const cMatrix33<T>& m);
template <typename T>
cVector3<T>  operator*(const cMatrix33<T>& m, const cVector3<T>& v);
template <typename T>
cVector3<T>  operator*(const cVector3<T>& v, const cMatrix33<T>& m);

} // namespace gem

#endif
