/***********************************************************************
 *  File:       cVector3.cpp
 *
 *  Purpose:    Implementation of a 3D vector class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "macro.hpp"
#include "cVector3.hpp"

#include <iomanip>

namespace gem {

template <typename T>
cVector3<T>::cVector3(T v1, T v2, T v3)
{
    _v[0] = v1;
    _v[1] = v2;
    _v[2] = v3;
}

template <typename T>
cVector3<T>::cVector3(T v[3])
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
}

template <typename T>
cVector3<T>::cVector3(const cVector3& v)
{
    _v[0] = v[0];
    _v[1] = v[1];
    _v[2] = v[2];
}

template <typename T>
cVector3<T>& cVector3<T>::operator=(T value)
{
    _v[0] = value;
    _v[1] = value;
    _v[2] = value;

    return *this;
}
template <typename T>
const cVector3<T> cVector3<T>::operator+(T value) const
{
    return cVector3(_v[0] + value,
                    _v[1] + value,
                    _v[2] + value);
}

template <typename T>
cVector3<T>& cVector3<T>::operator+=(T value)
{
    _v[0] += value;
    _v[1] += value;
    _v[2] += value;

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator+(const cVector3& v) const
{
    return cVector3(_v[0] + v[0],
                    _v[1] + v[1],
                    _v[2] + v[2]);
}

template <typename T>
cVector3<T>& cVector3<T>::operator+=(const cVector3& v)
{
    _v[0] += v[0];
    _v[1] += v[1];
    _v[2] += v[2];

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator-(T value) const
{
    return cVector3(_v[0] - value,
                    _v[1] - value,
                    _v[2] - value);
}

template <typename T>
cVector3<T>& cVector3<T>::operator-=(T value)
{
    _v[0] -= value;
    _v[1] -= value;
    _v[2] -= value;

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator-(const cVector3& v) const
{
    return cVector3(_v[0] - v[0],
                    _v[1] - v[1],
                    _v[2] - v[2]);
}

template <typename T>
cVector3<T>& cVector3<T>::operator-=(const cVector3& v)
{
    _v[0] -= v[0];
    _v[1] -= v[1];
    _v[2] -= v[2];

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator*(T value) const
{
    return cVector3(_v[0] * value,
                    _v[1] * value,
                    _v[2] * value);
}

template <typename T>
cVector3<T>& cVector3<T>::operator*=(T value)
{
    _v[0] *= value;
    _v[1] *= value;
    _v[2] *= value;

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator*(const cVector3& v) const
{
    return cVector3(_v[0] * v[0],
                    _v[1] * v[1],
                    _v[2] * v[2]);
}

template <typename T>
cVector3<T>& cVector3<T>::operator*=(const cVector3& v)
{
    _v[0] *= v[0];
    _v[1] *= v[1];
    _v[2] *= v[2];

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator/(T value) const
{
    assert(value != 0);

    return cVector3(_v[0] / value,
                    _v[1] / value,
                    _v[2] / value);
}

template <typename T>
cVector3<T>& cVector3<T>::operator/=(T value)
{
    assert(value != 0);

    _v[0] /= value;
    _v[1] /= value;
    _v[2] /= value;

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator/(const cVector3& v) const
{
    assert((v[0] != 0) && (v[1] != 0) && (v[2] != 0));

    return cVector3(_v[0] / v[0],
                    _v[1] / v[1],
                    _v[2] / v[2]);
}

template <typename T>
cVector3<T>& cVector3<T>::operator/=(const cVector3& v)
{
    assert((v[0] != 0) && (v[1] != 0) && (v[2] != 0));

    _v[0] /= v[0];
    _v[1] /= v[1];
    _v[2] /= v[2];

    return *this;
}

template <typename T>
const cVector3<T> cVector3<T>::operator-(void) const
{
    return cVector3(-_v[0],
                    -_v[1],
                    -_v[2]);
}

template <typename T>
bool cVector3<T>::operator==(const cVector3& v) const
{
    return (_v[0] == v[0]) && (_v[1] == v[1]) && (_v[2] == v[2]);
}

template <typename T>
bool cVector3<T>::operator!=(const cVector3& v) const
{
    return (_v[0] != v[0]) || (_v[1] != v[1]) || (_v[2] != v[2]);
}

template <typename T>
bool cVector3<T>::operator<(const cVector3& v) const
{
    return (_v[0] < v[0]) && (_v[1] < v[1]) && (_v[2] < v[2]);
}

template <typename T>
bool cVector3<T>::operator<=(const cVector3& v) const
{
    return (_v[0] <= v[0]) && (_v[1] <= v[1]) && (_v[2] <= v[2]);
}

template <typename T>
bool cVector3<T>::operator>(const cVector3& v) const
{
    return (_v[0] > v[0]) && (_v[1] > v[1]) && (_v[2] > v[2]);
}

template <typename T>
bool cVector3<T>::operator>=(const cVector3& v) const
{
    return (_v[0] >= v[0]) && (_v[1] >= v[1]) && (_v[2] >= v[2]);
}

template <typename T>
T cVector3<T>::operator|(const cVector3& v) const
{
    return _v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2];
}

template <typename T>
const cVector3<T> cVector3<T>::operator^(const cVector3& v) const
{
    return cVector3(_v[1]*v[2] - _v[2]*v[1],
                    _v[2]*v[0] - _v[0]*v[2],
                    _v[0]*v[1] - _v[1]*v[0]);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cVector3<T>& v)
{
    std::cout << std::setprecision(4) << std::fixed
              << "[ "
              << std::setw(10) << v[0]
              << std::setw(10) << v[1]
              << std::setw(10) << v[2]
              << " ]";

    return os;
}

template <typename T>
T cVector3<T>::getMax(void) const
{
    return std::max(std::max(_v[0],_v[1]),_v[2]);
}

template <typename T>
T cVector3<T>::getMin(void) const
{
    return std::min(std::min(_v[0],_v[1]),_v[2]);
}

template <typename T>
T cVector3<T>::getSum(void) const
{
    return _v[0]+_v[1]+_v[2];
}

template <typename T>
T cVector3<T>::getSum2(void) const
{
    return pow2(_v[0])+pow2(_v[1])+pow2(_v[2]);
}

template <typename T>
T cVector3<T>::getLength(void) const
{
    return (T) std::sqrt(pow2(_v[0])+pow2(_v[1])+pow2(_v[2]));
}

template <typename T>
T cVector3<T>::getProduct(void) const
{
    return _v[0]*_v[1]*_v[2];
}

template <typename T>
void cVector3<T>::ind2sub(const cVector3<T>& size, T idx)
{
    assert(idx < size.getProduct());

    _v[0] = (T) std::floor(idx / (size[1]*size[2]));
    _v[1] = (T) std::floor(idx / size[2] - _v[0]*size[1]);
    _v[2] = idx - _v[0]*size[1]*size[2] - _v[1]*size[2];
}

template <typename T>
T cVector3<T>::sub2ind(const cVector3<T>& size) const
{
    assert(*this < size);

    return _v[0]*size[1]*size[2] + _v[1]*size[2] + _v[2];
}

template <typename T>
cVector3<T> operator*(T value, const cVector3<T>& v)
{
    return v * value;
}

// instantiation
template class cVector3<int32_t >;
template class cVector3<uint32_t>;
template class cVector3<int64_t >;
template class cVector3<uint64_t>;
template class cVector3<float   >;
template class cVector3<double  >;

template std::ostream& operator<<(std::ostream& os, const cVector3<int32_t >& v);
template std::ostream& operator<<(std::ostream& os, const cVector3<uint32_t>& v);
template std::ostream& operator<<(std::ostream& os, const cVector3<int64_t >& v);
template std::ostream& operator<<(std::ostream& os, const cVector3<uint64_t>& v);
template std::ostream& operator<<(std::ostream& os, const cVector3<float   >& v);
template std::ostream& operator<<(std::ostream& os, const cVector3<double  >& v);

template cVector3<int32_t > operator*(int32_t  value, const cVector3<int32_t >& v);
template cVector3<uint32_t> operator*(uint32_t value, const cVector3<uint32_t>& v);
template cVector3<int64_t > operator*(int64_t  value, const cVector3<int64_t >& v);
template cVector3<uint64_t> operator*(uint64_t value, const cVector3<uint64_t>& v);
template cVector3<float   > operator*(float    value, const cVector3<float   >& v);
template cVector3<double  > operator*(double   value, const cVector3<double  >& v);

} // namespace gem
