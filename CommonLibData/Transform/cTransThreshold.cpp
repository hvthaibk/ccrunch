/***********************************************************************
 *  File:       cTransThreshold.cpp
 *
 *  Purpose:    Implementation of a thresholding transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTransThreshold.hpp"

#include <algorithm>

namespace gem {

template <typename T>
T cTransThreshold<T>::getCutOffValue(const cData3<T>& dataSrc, T voxelVolume, T mapVolume, eWaterDir water)
{
    const std::string   funcName("T cTransThreshold<T>::getCutOffValue("
                                    "const cData3<T>& dataSrc, "
                                    "T voxelVolume, "
                                    "T mapVolume, "
                                    "eWaterDir water)");

    dataSrc.requireNonEmpty(funcName);

    std::vector<T>  dataTmp(dataSrc.getAddrData(),
                            dataSrc.getAddrData() + dataSrc.getNelement());

    switch (water) {
        case WATER_DIR_RISE:
            std::sort(dataTmp.begin(), dataTmp.end(), compLesser<T>);
            break;
        case WATER_DIR_FALL:
            std::sort(dataTmp.begin(), dataTmp.end(), compGreater<T>);
            break;
        default:
            ERROR(funcName, "unsupported water mode");
    }

    size_t          numeVoxel = (size_t) round(mapVolume / voxelVolume);

    return dataTmp[numeVoxel-1];
}

template <typename T>
void cTransThreshold<T>::setValueBelowCutOff(cData3<T>& dataSrc, T cutoff, T value)
{
    const std::string   funcName("void cTransThreshold<T>::setValueBelowCutOff("
                                    "cData3<T>& dataSrc, T cutoff, T value)");

    dataSrc.requireNonEmpty(funcName);

    #pragma omp parallel for
    for (size_t i = 0; i < dataSrc.getNelement(); i++) {

        dataSrc[i] = (dataSrc[i] < cutoff) ? value : dataSrc[i];
    }
}

template <typename T>
void cTransThreshold<T>::setValueAboveCutOff(cData3<T>& dataSrc, T cutoff, T value)
{
    const std::string   funcName("void cTransThreshold<T>::setValueAboveCutOff("
                                    "cData3<T>& dataSrc, T cutoff, T value)");

    dataSrc.requireNonEmpty(funcName);

    #pragma omp parallel for
    for (size_t i = 0; i < dataSrc.getNelement(); i++) {

        dataSrc[i] = (dataSrc[i] > cutoff) ? value : dataSrc[i];
    }
}

template <typename T>
void cTransThreshold<T>::threshold(cData3<T>& dataSrc, T cutoff)
{
    const std::string   funcName("void cTransThreshold<T>::threshold("
                                    "cData3<T>& dataSrc, T cutoff)");

    dataSrc.requireNonEmpty(funcName);

    array_threshold(dataSrc.getAddrData(), dataSrc.getNelement(), cutoff);
}

template <typename T>
void cTransThreshold<T>::threshold(const cData3<T>& dataSrc, cData3<T>& dataDst, T cutoff)
{
    const std::string   funcName("void cTransThreshold<T>::threshold("
                                    "const cData3<T>& dataSrc, "
                                    "cData3<T>& dataDst, "
                                    "T cutoff)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    array_threshold(dataSrc.getAddrData(), dataSrc.getNelement(), cutoff,
                    dataDst.getAddrData());
}

// instantiation
template class cTransThreshold<float >;
template class cTransThreshold<double>;

} // namespace gem
