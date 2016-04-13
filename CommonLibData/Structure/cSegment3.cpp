/***********************************************************************
 *  File:       cSegment3.cpp
 *
 *  Purpose:    Implementation of a 3D segment class
 *
 *  Author:     Lucia Martin Reixach, Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Lucia Martin Reixach, Thai V. Hoang, INRIA
 **********************************************************************/

#include "cSegment3.hpp"
#include "macro.hpp"

namespace gem {

template <typename T>
cSegment3<T>::cSegment3(void)
{
    _p[0] = 0;
    _p[1] = 0;
}

template <typename T>
cSegment3<T>::cSegment3(const cVector3<T>& p1, const cVector3<T>& p2)
{
    _p[0] = p1;
    _p[1] = p2;
}

template <typename T>
cSegment3<T>::cSegment3(const cSegment3& other)
{
    _p[0] = other._p[0];
    _p[1] = other._p[1];
}

template <typename T>
cSegment3<T>& cSegment3<T>::operator=(const cSegment3& other)
{
    _p[0] = other._p[0];
    _p[1] = other._p[1];

    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cSegment3<T>& segment)
{
    std::cout << segment._p[0] << std::endl;
    std::cout << segment._p[1] << std::endl;

    return os;
}

template <typename T>
T cSegment3<T>::getAngle2D(void) const
{
    return std::atan2(_p[1][1]-_p[0][1], _p[1][0]-_p[0][0]);
}

template <typename T>
T cSegment3<T>::getLength2D(void) const
{
    return dist(_p[0][0], _p[1][0], _p[0][1], _p[1][1]);
}

template <typename T>
T cSegment3<T>::getLength3D(void) const
{
    return dist(_p[0][0], _p[1][0], _p[0][1], _p[1][1], _p[0][2], _p[1][2]);
}

template <typename T>
bool cSegment3<T>::isCoincide (const cSegment3& other) const
{
    T   m1, h1, m2, h2;   // y = m*x + h

    m1 = (_p[1][1] - _p[0][1]) / (_p[1][0] - _p[0][0]);
    h1 = _p[0][1] - m1 * _p[0][0];
    m2 = (other._p[1][1] - other._p[0][1]) / (other._p[1][0] - other._p[0][0]);
    h2 = other._p[0][1] - m2 * other._p[0][0];

    if (std::abs(m2 - m1) < GEM_PRECISON && std::abs(h2 - h1) < GEM_PRECISON)
        return true;
    else
        return false;
}

template <typename T>
bool cSegment3<T>::isIntersect(const cSegment3& other) const
{
    cVector3<T>     r, s, dif, t, u;
    T               prod;

    r = _p[1] - _p[0];
    s = other._p[1] - other._p[0];

    prod = (r ^ s)[2];
    dif = other._p[0] - _p[0];

    if (prod == 0)
        return false;
    else
    {
        t = (dif ^ s) / prod;
        u = (dif ^ r) / prod;

        if ((t[2] > 0 && t[2] < 1) && (u[2] > 0 && u[2] < 1))
            return true;
        else
            return false;
    }
}

template <typename T>
bool cSegment3<T>::isParallel (const cSegment3& other) const
{
    return ((_p[1] - _p[0]) ^ (other._p[1] - other._p[0])).getLength() <= GEM_PRECISON;
}

template <typename T>
cVector3<T> cSegment3<T>::getIntersectPoint(const cSegment3& other) const
{
    cVector3<T>     r, s;
    T               t;

    r = _p[1] - _p[0];
    s = other._p[1] - other._p[0];

    t = ((other._p[0] - _p[0]) ^ s).getLength() / (r ^ s).getLength();

    return _p[0] + (r * t);
}

template <typename T>
bool cSegment3<T>::isTouched (const cSegment3& other) const
{
    T   m, h;   // y = m*x + h

    for (size_t i = 0; i < 2; i++)
    {
       if ((_p[0][0] <= other._p[i][0] && other._p[i][0] <= _p[1][0]) || (_p[1][0] <= other._p[i][0] && other._p[i][0] <= _p[0][0]))
        {
            m = (_p[1][1] - _p[0][1]) / (_p[1][0] - _p[0][0]);
            h = _p[0][1] - m * _p[0][0];

            if (m * other._p[i][0] + h == other._p[i][1])
                return true;
        }
    }
    return false;
}

// instantiation
template class cSegment3<float >;
template class cSegment3<double>;

template std::ostream& operator<<(std::ostream& os, const cSegment3<float >& atom);
template std::ostream& operator<<(std::ostream& os, const cSegment3<double>& atom);

} // namespace gem
