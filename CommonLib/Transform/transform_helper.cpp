/***********************************************************************
 *  File:       transform_helper.cpp
 *
 *  Purpose:    Implementation of transformation functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "transform.hpp"

namespace gem {

/**************************
 * Helper functions
 *************************/

// CONVENTIONS
//    - For a rotation on a certain plane (XY), the rotation axis is
//      the remaining axis (Z).
//    - RIGHT is the rotating direction of this axis (looking from
//      origin to its positive end-point) as screwing a screwdriver.
//    - Opposite for LEFT.
//
// So for images:
//    - RIGHT: rotate counter-clockwise
//    - LEFT : rotate         clockwise

// 2D rotation matrix
// SOURCE: http://en.wikipedia.org/wiki/Rotation_matrix
template <typename T>
void transform_rotmat(T theta,
                      T& m11, T& m12,
                      T& m21, T& m22,
                      eRot2D order = ROT2D_RIGHT)
{
    T    sTheta = std::sin(theta);
    T    cTheta = std::cos(theta);

    switch (order) {
        case ROT2D_RIGHT:
            m11 =  cTheta;  m12 = -sTheta;
            m21 =  sTheta;  m22 =  cTheta;
            break;
        default:
            ERROR("transform_rotmat", "unsupported rotation mode");
    }
}

// instantiation
template
void transform_rotmat<float >(float  theta,
                              float&  m11, float&  m12,
                              float&  m21, float&  m22,
                              eRot2D order);
template
void transform_rotmat<double>(double theta,
                              double& m11, double& m12,
                              double& m21, double& m22,
                              eRot2D order);

// 3D rotation matrix
// SOURCE: http://en.wikipedia.org/wiki/Euler_angles
template <typename T>
void transform_rotmat(T alpha, T beta, T gamma,
                      T& m11,  T& m12, T& m13,
                      T& m21,  T& m22, T& m23,
                      T& m31,  T& m32, T& m33,
                      eRot3D order)
{
    T    sAl = std::sin(alpha);
    T    cAl = std::cos(alpha);
    T    sBe = std::sin(beta);
    T    cBe = std::cos(beta);
    T    sGa = std::sin(gamma);
    T    cGa = std::cos(gamma);

    switch (order) {
        case ROT3D_RIGHT_ZYZ:
            m11 = cAl*cBe*cGa - sAl*sGa;
            m12 = -cGa*sAl - cAl*cBe*sGa;
            m13 = cAl*sBe;
            m21 = cAl*sGa + cBe*cGa*sAl;
            m22 = cAl*cGa - cBe*sAl*sGa;
            m23 = sAl*sBe;
            m31 = -cGa*sBe;
            m32 = sBe*sGa;
            m33 = cBe;
            break;
        default:
            ERROR("transform_rotmat", "unsupported rotation mode");
    }
}

// instantiation
template
void transform_rotmat<float >(float  alpha, float  beta, float  gamma,
                              float&  m11,  float&  m12, float&  m13,
                              float&  m21,  float&  m22, float&  m23,
                              float&  m31,  float&  m32, float&  m33,
                              eRot3D order);
template
void transform_rotmat<double>(double alpha, double beta, double gamma,
                              double& m11,  double& m12, double& m13,
                              double& m21,  double& m22, double& m23,
                              double& m31,  double& m32, double& m33,
                              eRot3D order);

// center of coordinate
template <typename T>
T transform_centerCoord(size_t nRow)
{
    return (nRow%2) ? ((T) (nRow>>1)) : ((T) (nRow>>1) - (T) 0.5);
}

template <typename T>
T transform_centerCoord(T cSrc, T factor)
{
    return factor * cSrc;
}

// instantiation
template
float  transform_centerCoord<float >(size_t nRow);
template
double transform_centerCoord<double>(size_t nRow);
template
float  transform_centerCoord<float >(float  cSrc, float  factor);
template
double transform_centerCoord<double>(double cSrc, double factor);

// scaled size
template <typename T>
size_t transform_scaleSize(size_t nPixSrc, T factor)
{
    return nPixSrc > 1 ? (size_t) std::ceil(factor * (T) nPixSrc) : 1;
}

template <typename T>
size_t transform_scaleSize(size_t nPixSrc, T cSrc, T factor)
{
    ptrdiff_t   nPixDst1 = (ptrdiff_t) std::ceil(factor * (              cSrc));
    ptrdiff_t   nPixDst2 = (ptrdiff_t) std::ceil(factor * ((T) nPixSrc - cSrc));

    return (size_t) (nPixDst1 + nPixDst2);
}

// instantiation
template
size_t transform_scaleSize<float >(size_t nPixSrc, float  factor);
template
size_t transform_scaleSize<double>(size_t nPixSrc, double factor);
template
size_t transform_scaleSize<float >(size_t nPixSrc, float  cSrc, float  factor);
template
size_t transform_scaleSize<double>(size_t nPixSrc, double cSrc, double factor);

// binned size
size_t transform_binSize(size_t nPixSrc, size_t factor)
{
    return nPixSrc > 1 ? nPixSrc / factor : 1;
}

} // namespace gem
