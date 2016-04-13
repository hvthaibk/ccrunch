/***********************************************************************
 *  File:       cCoordinate.cpp
 *
 *  Purpose:    Implementation of a coordinate system class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cCoordinate.hpp"

namespace gem {

template <typename T>
void cCoordinate<T>::computeAngle(const cData3<T>& X,
                                  const cData3<T>& Y,
                                  cData3<T>& theta)
{
    const std::string   funcName("void cCoordinate<T>::computeAngle("
                                    "const cData3<T>& X, "
                                    "const cData3<T>& Y, "
                                    "cData3<T>& theta)");

    require(X.getSize() == Y.getSize(),
            funcName + ": input data are not compatible");
    require(X.getNrow() > 1 && X.getNcol() > 1 && X.getNsec() == 1,
            funcName + ": input data is not 2D");

    theta.memReAlloc(X.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < X.getNelement(); i++) {
        theta[i] = std::atan2(Y[i],X[i]);
    }
}

template <typename T>
void cCoordinate<T>::computeDistance(const cData3<T>& X,
                                     const cData3<T>& Y,
                                     cData3<T>& dist)
{
    const std::string   funcName("void cCoordinate<T>::computeDistance("
                                    "const cData3<T>& X, "
                                    "const cData3<T>& Y, "
                                    "cData3<T>& dist)");

    require(X.getSize() == Y.getSize(),
            funcName + ": input data are not compatible");
    require(X.getNrow() > 1 && X.getNcol() > 1 && X.getNsec() == 1,
            funcName + ": input data is not 2D");

    dist.memReAlloc(X.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < X.getNelement(); i++) {
        dist[i] = std::sqrt(pow2(X[i]) + pow2(Y[i]));
    }
}

template <typename T>
void cCoordinate<T>::computeDistanceSquared(const cData3<T>& X,
                                            const cData3<T>& Y,
                                            cData3<T>& dist2)
{
    const std::string   funcName("void cCoordinate<T>::computeDistanceSquared("
                                    "const cData3<T>& X, "
                                    "const cData3<T>& Y, "
                                    "cData3<T>& dist2)");

    require(X.getSize() == Y.getSize(),
            funcName + ": input data are not compatible");
    require(X.getNrow() > 1 && X.getNcol() > 1 && X.getNsec() == 1,
            funcName + ": input data is not 2D");

    dist2.memReAlloc(X.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < X.getNelement(); i++) {
        dist2[i] = pow2(X[i]) + pow2(Y[i]);
    }
}

template <typename T>
void cCoordinate<T>::computeDistance(const cData3<T>& X,
                                     const cData3<T>& Y,
                                     const cData3<T>& Z,
                                     cData3<T>& dist)
{
    const std::string   funcName("void cCoordinate<T>::computeDistance("
                                    "const cData3<T>& X, "
                                    "const cData3<T>& Y, "
                                    "const cData3<T>& Z, "
                                    "cData3<T>& dist)");

    require(X.getSize() == Y.getSize() &&
            X.getSize() == Z.getSize(),
            funcName + ": input data are not compatible");
    require(X.getNrow() > 1 && X.getNcol() > 1 && X.getNsec() == 1,
            funcName + ": input data is not 2D");

    dist.memReAlloc(X.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < X.getNelement(); i++) {
        dist[i] = std::sqrt(pow2(X[i]) + pow2(Y[i]) + pow2(Z[i]));
    }
}

template <typename T>
void cCoordinate<T>::computeDistanceSquared(const cData3<T>& X,
                                            const cData3<T>& Y,
                                            const cData3<T>& Z,
                                            cData3<T>& dist2)
{
    const std::string   funcName("void cCoordinate<T>::computeDistanceSquared("
                                    "const cData3<T>& X, "
                                    "const cData3<T>& Y, "
                                    "const cData3<T>& Z, "
                                    "cData3<T>& dist2)");

    require(X.getSize() == Y.getSize() &&
            X.getSize() == Z.getSize(),
            funcName + ": input data are not compatible");
    require(X.getNrow() > 1 && X.getNcol() > 1 && X.getNsec() == 1,
            funcName + ": input data is not 2D");

    dist2.memReAlloc(X.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < X.getNelement(); i++) {
        dist2[i] = pow2(X[i]) + pow2(Y[i]) + pow2(Z[i]);
    }
}

// instantiation
template class cCoordinate<float >;
template class cCoordinate<double>;

} // namespace gem
