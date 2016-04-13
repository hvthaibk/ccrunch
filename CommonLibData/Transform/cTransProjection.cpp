/***********************************************************************
 *  File:       cTransProjection.cpp
 *
 *  Purpose:    Implementation of a 2D/3D projection class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTransProjection.hpp"

namespace gem {

template <typename T>
void cTransProjection<T>::project2D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& theta, eInter inter)
{
    const std::string   funcName("void cTransProjection::project2D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& theta, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(dataSrc.getDimension() == 2, funcName + ": input data is not 2D");
    require(theta.size() > 0, funcName + ": no angle to project");

    dataDst.memReAlloc(cVector3<size_t>(theta.size(),dataSrc.getNcol(),1));

    transform_project(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                      dataDst.getAddrData(),
                      theta,
                      inter);
}

template <typename T>
void cTransProjection<T>::project3D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& alpha, std::vector<T>& beta, eInter inter)
{
    const std::string   funcName("void cTransProjection<T>::project3D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& alpha, std::vector<T>& beta, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(dataSrc.getDimension() == 3, funcName + ": input data is not 3D");
    require(alpha.size() > 0, funcName + ": no angle to project");
    require(alpha.size() == beta.size(), funcName + ": alpha and beta do not have the same size");

    dataDst.memReAlloc(cSize3(alpha.size(),dataSrc.getNcol(),dataSrc.getNsec()));

    transform_project(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),dataSrc.getNsec(),
                      dataDst.getAddrData(),
                      alpha, beta,
                      inter);
}

template <typename T>
void cTransProjection<T>::trace2D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& theta, eTrace trace, eInter inter)
{
    const std::string   funcName("void cTransProjection<T>::trace2D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& theta, eTrace trace, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(dataSrc.getDimension() == 2, funcName + ": input data is not 2D");
    require(theta.size() > 0, funcName + ": no angle to project");

    dataDst.memReAlloc(cSize3(theta.size(),dataSrc.getNcol(),1));

    transform_trace(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                    dataDst.getAddrData(),
                    theta,
                    trace,
                    inter);
}

template <typename T>
void cTransProjection<T>::trace3D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& alpha, std::vector<T>& beta, std::vector<T>& gamma, eTrace trace, eInter inter)
{
    const std::string   funcName("void cTransProjection<T>::trace3D(const cData3<T>& dataSrc, cData3<T>& dataDst, std::vector<T>& alpha, std::vector<T>& beta, std::vector<T>& gamma, eTrace trace, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(dataSrc.getDimension() == 3, funcName + ": input data is not 3D");
    require(alpha.size() > 0, funcName + ": no angle to project");
    require(alpha.size() == beta.size(), funcName + ": alpha and beta do not have the same size");

    dataDst.memReAlloc(cSize3(alpha.size(),dataSrc.getNrow(),dataSrc.getNcol()));

    transform_trace(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),dataSrc.getNsec(),
                    dataDst.getAddrData(),
                    alpha, beta, gamma,
                    trace,
                    inter);
}

// instantiation
template class cTransProjection<float >;
template class cTransProjection<double>;

} // namespace gem
