/***********************************************************************
 *  File:       cTransMorpho.cpp
 *
 *  Purpose:    Implementation of a transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTransMorpho.hpp"

namespace gem {

template <typename T>
void cTransMorpho<T>::dilate(cData3<T>& dataSrc, T radius)
{
    const std::string   funcName("void cTransMorpho<T>::dilate(cData3<T>& dataSrc, T radius)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>    other(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(),
                               other.getAddrData(),
                             radius);
            break;
        case 2:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                               other.getAddrData(),
                             radius);
            break;
        case 3:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                               other.getAddrData(),
                             radius);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(other);
}

template <typename T>
void cTransMorpho<T>::dilate(const cData3<T>& dataSrc, cData3<T>& dataDst, T radius)
{
    const std::string   funcName("void cTransMorpho<T>::dilate(const cData3<T>& dataSrc, cData3<T>& dataDst, T radius)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(),
                             dataDst.getAddrData(),
                             radius);
            break;
        case 2:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                             dataDst.getAddrData(),
                             radius);
            break;
        case 3:
            transform_dilate(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                             dataDst.getAddrData(),
                             radius);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

// instantiation
template class cTransMorpho<float >;
template class cTransMorpho<double>;

} // namespace gem
