/***********************************************************************
 *  File:       cTransDistance.cpp
 *
 *  Purpose:    Implementation of a distance transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTransDistance.hpp"

namespace gem {

template <typename T>
void cTransDistance<T>::distanceGeneric(const cData3<T>& dataSrc, cData3<T>& dataDst)
{
    const std::string   funcName("void cTransDistance<T>::distanceGeneric(const cData3<T>& dataSrc, cData3<T>& dataDst)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            transform_distance(dataSrc.getAddrData(), dataSrc.getNrow(),
                               dataDst.getAddrData());
            break;
        case 2:
            transform_distance(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                               dataDst.getAddrData());
            break;
        case 3:
            transform_distance(dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                               dataDst.getAddrData());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataDst.opMathSqrt();
}

template <typename T>
void cTransDistance<T>::distanceNormal(const cData3<T>& dataSrc, cData3<T>& dataDst)
{
    const std::string   funcName("void cTransDistance<T>::distanceNormal(const cData3<T>& dataSrc, cData3<T>& dataDst)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    cData3<T>       dataTmp(dataSrc.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < dataTmp.getNelement(); i++) {
        if (dataSrc[i] > 0) {
            dataTmp[i] = 0;
        }
        else if (dataSrc[i] == 0) {
            dataTmp[i] = (T) GEM_DISTANCE_MAX;
        }
        else {
            ERROR(funcName, "dataSrc elements must be >= 0");
        }
    }

    dataDst.memReAlloc(dataTmp.getSize());

    switch (dataTmp.getDimension()) {
        case 1:
            transform_distance(dataTmp.getAddrData(), dataTmp.getNrow(),
                               dataDst.getAddrData());
            break;
        case 2:
            transform_distance(dataTmp.getAddrData(), dataTmp.getNrow(), dataTmp.getNcol(),
                               dataDst.getAddrData());
            break;
        case 3:
            transform_distance(dataTmp.getAddrData(), dataTmp.getNrow(), dataTmp.getNcol(), dataTmp.getNsec(),
                               dataDst.getAddrData());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataDst.opMathSqrt();
}

template <typename T>
void cTransDistance<T>::voronoi(const cData3<T>& dataSrc, const cArray3<size_t>& labelSrc,
                                      cData3<T>& dataDst,       cArray3<size_t>& labelDst)
{
    const std::string   funcName("void cTransDistance<T>::voronoi(const cData3<T>& dataSrc, const cArray3<size_t>& labelSrc, cData3<T>& dataDst, cArray3<size_t>& labelDst)");

    dataSrc.requireNonEmpty(funcName);
    labelSrc.requireNonEmpty(funcName);
    require(&dataSrc  != &dataDst,  funcName + ": two objects are the same");
    require(&labelSrc != &labelDst, funcName + ": two objects are the same");
    require(dataSrc.getSize() == labelSrc.getSize(), funcName + ": data and label do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());
    labelDst.memReAlloc(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            transform_distance( dataSrc.getAddrData(), dataSrc.getNrow(),
                                dataDst.getAddrData(),
                               labelSrc.getAddrData(),
                               labelDst.getAddrData());
            break;
        case 2:
            transform_distance( dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                                dataDst.getAddrData(),
                               labelSrc.getAddrData(),
                               labelDst.getAddrData());
            break;
        case 3:
            transform_distance( dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                dataDst.getAddrData(),
                               labelSrc.getAddrData(),
                               labelDst.getAddrData());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataDst.opMathSqrt();
}

template <typename T>
void cTransDistance<T>::watershed(const cData3<T>& dataSrc, cArray3<size_t>& labelDst,
                                  T cutoff, T seedDistMin, eWaterDir water)
{
    const std::string   funcName("void cTransDistance<T>::watershed(const cData3<T>& dataSrc, cArray3<size_t>& labelDst, T cutoff, T seedDistMin, eWaterDir water)");

    dataSrc.requireNonEmpty(funcName);
    require(seedDistMin > 0, funcName + ": seedDistMin is not positive");

    labelDst.memReAllocZero(dataSrc.getSize());

    std::vector< std::set<size_t> > neighborList;

    transform_watershed( dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                        labelDst.getAddrData(),
                         dataSrc.getDimension(),
                        cutoff,
                        seedDistMin,
                        water,
                        neighborList);

    neighborList.clear();
}

template <typename T>
void cTransDistance<T>::watershedReverseLabel(cArray3<size_t>& labelSrc)
{
    const std::string   funcName("void cTransDistance<T>::watershedReverseLabel(cArray3<size_t>& labelSrc)");

    labelSrc.requireNonEmpty(funcName);

    // reverse the label values so that
    // max label corresponds to the peak of interest
    size_t    maxLabel = labelSrc.getMax();

    #pragma omp parallel for
    for (size_t i = 0; i < labelSrc.getNelement(); i++) {

        if (labelSrc[i] > 0) {
            labelSrc[i] = maxLabel - labelSrc[i] + 1;
        }
    }
}

template <typename T>
void cTransDistance<T>::label(const cData3<T>& dataSrc, cArray3<size_t>& labelDst)
{
    const std::string   funcName("void cTransDistance<T>::label(const cData3<T>& dataSrc, cArray3<size_t>& labelDst)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>   dataTmp1(dataSrc), dataTmp2;

    dataTmp1 *= (T) GEM_DISTANCE_MAX;

    distanceGeneric(dataTmp1, dataTmp2);

    watershed(dataTmp2, labelDst, (T) 0.5, (T) 2, WATER_DIR_FALL);
/*
    const std::string   funcName("void cTransDistance<T>::label(const cData3<T>& dataSrc, cArray3<size_t>& labelDst, T cutoff)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>   dataTmp;

    dataTmp.opThreshold(dataSrc, cutoff);

    labelDst.memReAlloc(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            transform_label( dataTmp.getAddrData(), dataTmp.getNrow(),
                            labelDst.getAddrData());
            break;
        case 2:
            transform_label( dataTmp.getAddrData(), dataTmp.getNrow(), dataTmp.getNcol(),
                            labelDst.getAddrData());
            break;
        case 3:
            transform_label( dataTmp.getAddrData(), dataTmp.getNrow(), dataTmp.getNcol(), dataTmp.getNsec(),
                            labelDst.getAddrData());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
 */
}

// instantiation
template class cTransDistance<float >;
template class cTransDistance<double>;

} // namespace gem
