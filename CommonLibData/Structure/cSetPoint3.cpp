/***********************************************************************
 *  File:       cSetPoint3.cpp
 *
 *  Purpose:    Implementation of a 3D point set class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cSetPoint3.hpp"

#include "transform.hpp"
#include "cMatrix33.hpp"

namespace gem {

template <typename T>
cSetPoint3<T>::cSetPoint3(const cSetPoint3& other)
    : std::vector< cVector3<T> >(other)
{
/*    resize(other.size());

    #pragma omp parallel for
    for (size_t i = 0; i < other.size(); i++) {
        at[i] = other[i];
    }*/
}

template <typename T>
const cVector3<T> cSetPoint3<T>::getMean(void) const
{
    requireNonEmpty();

    cVector3<T>     pointSum;

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
        pointSum += at(i);
    }

    return pointSum / (T) size();
}

template <typename T>
cSetPoint3<T>& cSetPoint3<T>::operator=(const cSetPoint3& other)
{
    other.requireNonEmpty();

    resize(other.size());

    #pragma omp parallel for
    for (size_t i = 0; i < other.size(); i++) {
       at(i) = other[i];
    }

    return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const cSetPoint3<T>& sp3)
{
    sp3.requireNonEmpty();

    for (size_t i = 0; i < sp3.size(); i++) {
        std::cout << sp3[i] << std::endl;
    }

    return os;
}

template <typename T>
void cSetPoint3<T>::opFlip(eFlipPlane plane)
{
    requireNonEmpty();

    switch (plane) {
        case FLIP_PLANE_YZ:
            #pragma omp parallel for
            for (size_t i = 0; i < size(); i++) {
               at(i)[0] = -at(i)[0];
            }
            break;
        case FLIP_PLANE_ZX:
            #pragma omp parallel for
            for (size_t i = 0; i < size(); i++) {
               at(i)[1] = -at(i)[1];
            }
            break;
        case FLIP_PLANE_XY:
            #pragma omp parallel for
            for (size_t i = 0; i < size(); i++) {
               at(i)[2] = -at(i)[2];
            }
            break;
        default:
            ERROR("cSetPoint3<T>::opFlip", "unsupported flip plane mode");
    }
}

template <typename T>
void cSetPoint3<T>::opRotate(const cVector3<T>& angle)
{
    requireNonEmpty();
    require(angle != cVector3<T>(0,0,0),
            "cSetPoint3<T>::opRotate(): rotate with zero angle");

    cVector3<T>     pointTmp;
    cMatrix33<T>    matRot;

    transform_rotmat(angle[0], angle[1], angle[2],
                     matRot[0][0], matRot[0][1], matRot[0][2],
                     matRot[1][0], matRot[1][1], matRot[1][2],
                     matRot[2][0], matRot[2][1], matRot[2][2],
                     ROT3D_RIGHT_ZYZ);

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
        pointTmp = matRot *at(i);;

       at(i) = pointTmp;
    }
}

template <typename T>
void cSetPoint3<T>::opScale(const cVector3<T>& factor)
{
    requireNonEmpty();
    require((factor[0] != 0) && (factor[1] != 0) && (factor[2] != 0),
            "cSetPoint3<T>::opScale(): scale with zero factor");

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
       at(i) /= factor;
    }
}

template <typename T>
void cSetPoint3<T>::opTranslate(const cVector3<T>& offset)
{
    requireNonEmpty();
    require(offset != cVector3<T>(0,0,0),
            "cSetPoint3<T>::opTranslate(): translate with zero offset");

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
       at(i) += offset;
    }
}

template <typename T>
T cSetPoint3<T>::computeSD(const cSetPoint3& other) const
{
    requireNonEmpty();
    require(size() == other.size(),
            "cSetPoint3<T>::computeSD(): two sets do not have the same number of points");

    T     sqrDistSum = 0;

    #pragma omp parallel for reduction(+:sqrDistSum)
    for (size_t i = 0; i < size(); i++) {
        sqrDistSum += pow2(at(i)[0] - other[i][0] ) +
                      pow2(at(i)[1] - other[i][1] ) +
                      pow2(at(i)[2] - other[i][2] );
    }

    return sqrDistSum;
}

template <typename T>
T cSetPoint3<T>::computeRMSD(const cSetPoint3& other) const
{
    return std::sqrt(computeSD(other) / (T) size());
}

// instantiation
/*template class cSetPoint3<int32_t >;
template class cSetPoint3<uint32_t>;
template class cSetPoint3<int64_t >;
template class cSetPoint3<uint64_t>;*/
template class cSetPoint3<float   >;
template class cSetPoint3<double  >;

//template std::ostream& operator<<(std::ostream& os, const cSetPoint3<float >& sp3);
//template std::ostream& operator<<(std::ostream& os, const cSetPoint3<double>& sp3);

} // namespace gem
