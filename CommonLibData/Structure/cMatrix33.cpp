/***********************************************************************
 *  File:       cMatrix33.cpp
 *
 *  Purpose:    Implementation of a 3x3 matrix class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "macro.hpp"
#include "cMatrix33.hpp"

#include <iomanip>

namespace gem {

template <typename T>
cMatrix33<T>::cMatrix33(T m11, T m12, T m13,
                        T m21, T m22, T m23,
                        T m31, T m32, T m33)
{
    _m[0][0] = m11;     _m[0][1] = m12;     _m[0][2] = m13;
    _m[1][0] = m21;     _m[1][1] = m22;     _m[1][2] = m23;
    _m[2][0] = m31;     _m[2][1] = m32;     _m[2][2] = m33;
}

template <typename T>
cMatrix33<T>::cMatrix33(const T m[9])
{
    _m[0][0] = m[0];    _m[0][1] = m[1];    _m[0][2] = m[2];
    _m[1][0] = m[3];    _m[1][1] = m[4];    _m[1][2] = m[5];
    _m[2][0] = m[6];    _m[2][1] = m[7];    _m[2][2] = m[8];
}

template <typename T>
cMatrix33<T>::cMatrix33(const cMatrix33& m)
{
    _m[0][0] = m[0][0]; _m[0][1] = m[0][1]; _m[0][2] = m[0][2];
    _m[1][0] = m[1][0]; _m[1][1] = m[1][1]; _m[1][2] = m[1][2];
    _m[2][0] = m[2][0]; _m[2][1] = m[2][1]; _m[2][2] = m[2][2];
}

template <typename T>
cMatrix33<T>::cMatrix33(const cVector3<T>& v1, const cVector3<T>& v2, const cVector3<T>& v3)
{
    _m[0] = v1;     _m[1] = v2;     _m[2] = v3;
}

template <typename T>
cMatrix33<T>::cMatrix33(const cVector3<T> v[3])
{
    _m[0] = v[0];   _m[1] = v[1];   _m[2] = v[2];
}

template <typename T>
T cMatrix33<T>::getLength(void) const
{
    return (T) std::sqrt(pow2(_m[0][0])+pow2(_m[0][1])+pow2(_m[0][2])+
                         pow2(_m[1][0])+pow2(_m[1][1])+pow2(_m[1][2])+
                         pow2(_m[2][0])+pow2(_m[2][1])+pow2(_m[2][2]));
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator=(T value)
{
    _m[0] = value;
    _m[1] = value;
    _m[2] = value;

    return *this;
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator=(const cMatrix33& m)
{
    _m[0] = m[0];
    _m[1] = m[1];
    _m[2] = m[2];

    return *this;
}

template <typename T>
void cMatrix33<T>::transpose(void)
{
    std::swap(_m[0][1], _m[1][0]);
    std::swap(_m[0][2], _m[2][0]);
    std::swap(_m[1][0], _m[0][1]);
    std::swap(_m[1][2], _m[2][1]);
    std::swap(_m[2][0], _m[0][2]);
    std::swap(_m[2][1], _m[1][2]);
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator+(T value) const
{
    return cMatrix33(_m[0] + value,
                     _m[1] + value,
                     _m[2] + value);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator+=(T value)
{
    _m[0] += value;
    _m[1] += value;
    _m[2] += value;

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator+(const cMatrix33& m) const
{
    return cMatrix33(_m[0] + m[0],
                     _m[1] + m[1],
                     _m[2] + m[2]);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator+=(const cMatrix33& m)
{
    _m[0] += m[0];
    _m[1] += m[1];
    _m[2] += m[2];

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator-(T value) const
{
    return cMatrix33(_m[0] - value,
                     _m[1] - value,
                     _m[2] - value);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator-=(T value)
{
    _m[0] -= value;
    _m[1] -= value;
    _m[2] -= value;

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator-(const cMatrix33& m) const
{
    return cMatrix33(_m[0] - m[0],
                     _m[1] - m[1],
                     _m[2] - m[2]);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator-=(const cMatrix33& m)
{
    _m[0] -= m[0];
    _m[1] -= m[1];
    _m[2] -= m[2];

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator*(T value) const
{
    return cMatrix33(_m[0] * value,
                     _m[1] * value,
                     _m[2] * value);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator*=(T value)
{
    _m[0] *= value;
    _m[1] *= value;
    _m[2] *= value;

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator*(const cMatrix33& m) const
{
    return cMatrix33(_m[0] * m[0],
                     _m[1] * m[1],
                     _m[2] * m[2]);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator*=(const cMatrix33& m)
{
    _m[0] *= m[0];
    _m[1] *= m[1];
    _m[2] *= m[2];

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator/(T value) const
{
    assert(value != 0);

    return cMatrix33(_m[0] / value,
                     _m[1] / value,
                     _m[2] / value);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator/=(T value)
{
    assert(value != 0);

    _m[0] /= value;
    _m[1] /= value;
    _m[2] /= value;

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator/(const cMatrix33& m) const
{
    return cMatrix33(_m[0] / m[0],
                     _m[1] / m[1],
                     _m[2] / m[2]);
}

template <typename T>
cMatrix33<T>& cMatrix33<T>::operator/=(const cMatrix33& m)
{
    _m[0] /= m[0];
    _m[1] /= m[1];
    _m[2] /= m[2];

    return *this;
}

template <typename T>
const cMatrix33<T> cMatrix33<T>::operator-(void) const
{
    return cMatrix33(-_m[0],
                     -_m[1],
                     -_m[2]);
}

template <typename T>
bool cMatrix33<T>::operator==(const cMatrix33& m) const
{
    return (_m[0] == m[0]) && (_m[1] == m[1]) && (_m[2] == m[2]);
}

template <typename T>
bool cMatrix33<T>::operator!=(const cMatrix33& m) const
{
    return (_m[0] != m[0]) || (_m[1] != m[1]) || (_m[2] != m[2]);
}

template <typename T>
bool cMatrix33<T>::operator<(const cMatrix33& m) const
{
    return (_m[0] < m[0]) && (_m[1] < m[1]) && (_m[2] < m[2]);
}

template <typename T>
bool cMatrix33<T>::operator<=(const cMatrix33& m) const
{
    return (_m[0] <= m[0]) && (_m[1] <= m[1]) && (_m[2] <= m[2]);
}

template <typename T>
bool cMatrix33<T>::operator>(const cMatrix33& m) const
{
    return (_m[0] > m[0]) && (_m[1] > m[1]) && (_m[2] > m[2]);
}

template <typename T>
bool cMatrix33<T>::operator>=(const cMatrix33& m) const
{
    return (_m[0] >= m[0]) && (_m[1] >= m[1]) && (_m[2] >= m[2]);
}

template <typename T>
T cMatrix33<T>::operator|(const cMatrix33& m) const
{
    return (_m[0] | m[0]) + (_m[1] | m[1]) + (_m[2] | m[2]);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cMatrix33<T>& m)
{
    std::cout << m[0] << std::endl
              << m[1] << std::endl
              << m[2];

    return os;
}

template <typename T>
cMatrix33<T> operator*(T value, const cMatrix33<T>& m)
{
    return cMatrix33<T>(value*m[0], value*m[1], value*m[2]);
}

template <typename T>
cVector3<T> operator*(const cMatrix33<T>& m, const cVector3<T>& v)
{
    return cVector3<T>(m[0] | v,
                       m[1] | v,
                       m[2] | v);
}

template <typename T>
cVector3<T> operator*(const cVector3<T>& v, const cMatrix33<T>& m)
{
    return cVector3<T>(v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0],
                       v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1],
                       v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2]);
}

// instantiation
template class cMatrix33<float >;
template class cMatrix33<double>;

template std::ostream& operator<<(std::ostream& os, const cMatrix33<float >& m);
template std::ostream& operator<<(std::ostream& os, const cMatrix33<double>& m);

template cMatrix33<float > operator*(float  value, const cMatrix33<float >& m);
template cMatrix33<double> operator*(double value, const cMatrix33<double>& m);

template cVector3<float >  operator*(const cMatrix33<float >& m, const cVector3<float >& v);
template cVector3<double>  operator*(const cMatrix33<double>& m, const cVector3<double>& v);

template cVector3<float >  operator*(const cVector3<float >& v, const cMatrix33<float >& m);
template cVector3<double>  operator*(const cVector3<double>& v, const cMatrix33<double>& m);

} // namespace gem
