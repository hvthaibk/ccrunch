/***********************************************************************
 *  File:       cuTransSimilarity.cpp
 *
 *  Purpose:    Implementation of a similarity transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuTransSimilarity.hpp"

namespace gem {

template <typename T>
void cuTransSimilarity<T>::extendRotate(cuData3<T>& dataSrc, eExpRot exprot)
{
    const std::string   funcName("void cuTransSimilarity<T>::extendRotate(cuData3<T>& dataSrc, eExpRot exprot)");

    dataSrc.requireNonEmpty(funcName);

    T*          dataTmp = NULL;
    size_t      nRowTmp, nColTmp, nSecTmp = 1;

    switch (dataSrc.getDimension()) {
        case 2:
            switch (exprot) {
                case EXP_ROT_VALID:
                    ERROR(funcName, "EXP_ROT_VALID is not supported");
                    break;
                case EXP_ROT_FULL:
                    cuda_transform_expand_rotate_full(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                        dataTmp,               nRowTmp,           nColTmp );
                    break;
            }
            break;
        case 3:
            switch (exprot) {
                case EXP_ROT_VALID:
                    ERROR(funcName, "EXP_ROT_VALID is not supported");
                    break;
                case EXP_ROT_FULL:
                    cuda_transform_expand_rotate_full(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                        dataTmp,               nRowTmp,           nColTmp,           nSecTmp );
                    break;
            }
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memFree();
    dataSrc.arrayHandle(dataTmp, nRowTmp, nColTmp, nSecTmp);
}

template <typename T>
void cuTransSimilarity<T>::translate(cuData3<T>& dataSrc, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::translate(cuData3<T>& dataSrc, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cuData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataSrc.getNrow(),
                                     offset[0],
                                     inter);
            break;
        case 2:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataSrc.getNrow(), dataSrc.getNcol(),
                                     offset[0],         offset[1],
                                     inter);
            break;
        case 3:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                     offset[0],         offset[1],         offset[2],
                                     inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(dataDst);
}

template <typename T>
void cuTransSimilarity<T>::translate(const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::translate(const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 1:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataDst.getNrow(),
                                     offset[0],
                                     inter);
            break;
        case 2:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataDst.getNrow(), dataDst.getNcol(),
                                     offset[0],         offset[1],
                                     inter);
            break;
        case 3:
            cuda_transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                     dataDst.getNrow(), dataDst.getNcol(), dataDst.getNsec(),
                                     offset[0],         offset[1],         offset[2],
                                     inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuTransSimilarity<T>::rotate(cuData3<T>& dataSrc, const cVector3<T>& angle, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::rotate(cuData3<T>& dataSrc, const cVector3<T>& angle, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cuData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 2:
            require(angle[1] == 0, funcName + ": unused angle != 0");
            require(angle[2] == 0, funcName + ": unused angle != 0");

            cuda_transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                  dataSrc.getNrow(), dataSrc.getNcol(),
                                  angle[0],
                                  inter);
            break;
        case 3:
            cuda_transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                  dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                  angle[0], angle[1], angle[2],
                                  inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(dataDst);
}

template <typename T>
void cuTransSimilarity<T>::rotate(const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& angle, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::rotate(const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& angle, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 2:
            require(angle[1] == 0, funcName + ": unused angle != 0");
            require(angle[2] == 0, funcName + ": unused angle != 0");

            cuda_transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                  dataDst.getNrow(), dataDst.getNcol(),
                                  angle[0],
                                  inter);
            break;
        case 3:
            cuda_transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                  dataDst.getNrow(), dataDst.getNcol(), dataDst.getNsec(),
                                  angle[0], angle[1], angle[2],
                                  inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuTransSimilarity<T>::scale(cuData3<T>& dataSrc, T factor, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::scale(cuData3<T>& dataSrc, T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(factor > 0, funcName + ": negative scaling factor");

    cuData3<T>    dataDst(cSize3(transform_scaleSize(dataSrc.getNrow(), factor),
                                 transform_scaleSize(dataSrc.getNcol(), factor),
                                 transform_scaleSize(dataSrc.getNsec(), factor)));

    switch (dataSrc.getDimension()) {
        case 1:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(),
                                 factor,
                                 inter);
            break;
        case 2:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(), dataSrc.getNcol(),
                                 factor,
                                 inter);
            break;
        case 3:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                 factor,
                                 inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(dataDst);
}

template <typename T>
void cuTransSimilarity<T>::scale(const cuData3<T>& dataSrc, cuData3<T>& dataDst, T factor, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::scale(const cuData3<T>& dataSrc, cuData3<T>& dataDst, T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(factor > 0, funcName + ": negative scaling factor");

    dataDst.memReAlloc(cSize3(transform_scaleSize(dataSrc.getNrow(), factor),
                              transform_scaleSize(dataSrc.getNcol(), factor),
                              transform_scaleSize(dataSrc.getNsec(), factor)));

    switch (dataDst.getDimension()) {
        case 1:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(),
                                 factor,
                                 inter);
            break;
        case 2:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(), dataSrc.getNcol(),
                                 factor,
                                 inter);
            break;
        case 3:
            cuda_transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                                 dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                 factor,
                                 inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuTransSimilarity<T>::rotTrans(cuData3<T>& dataSrc,
        const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::rotTrans(cuData3<T>& dataSrc, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cuData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 2:
            require(angle[1]  == 0, funcName + ": unused angle != 0");
            require(angle[2]  == 0, funcName + ": unused angle != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                            dataSrc.getNrow(), dataSrc.getNcol(),
                                             angle[0],
                                            offset[0],         offset[1],
                                            inter);
            break;
        case 3:
            cuda_transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                            dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                                             angle[0],          angle[1],          angle[2],
                                            offset[0],         offset[1],         offset[2],
                                            inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(dataDst);
}

template <typename T>
void cuTransSimilarity<T>::rotTrans(const cuData3<T>& dataSrc, cuData3<T>& dataDst,
        const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::rotTrans(const cuData3<T>& dataSrc, cuData3<T>& dataDst, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 2:
            require(angle[1]  == 0, funcName + ": unused angle != 0");
            require(angle[2]  == 0, funcName + ": unused angle != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                            dataDst.getNrow(), dataDst.getNcol(),
                                             angle[0],
                                            offset[0],         offset[1],
                                            inter);
            break;
        case 3:
            cuda_transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                            dataDst.getNrow(), dataDst.getNcol(), dataDst.getNsec(),
                                             angle[0],          angle[1],          angle[2],
                                            offset[0],         offset[1],         offset[2],
                                            inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuTransSimilarity<T>::interp(const cuData3<T>& R,
                                  const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                                  eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::interp("
                                    "const cuData3<T>& R, "
                                    "const cuData3<T>& dataSrc, cuData3<T>& dataDst, "
                                    "eInter inter)");

    require(dataSrc.getDimension() == 1,
            funcName + ": dataSrc is not 1D");
    require(dataSrc.getSize() == R.getSize(),
            funcName + ": dataSrc and R do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());

    cuda_transform_interp(      R.getAddrData(),
                          dataSrc.getAddrData(),
                          dataDst.getAddrData(),
                          R.getNrow(),
                          inter);
}

template <typename T>
void cuTransSimilarity<T>::interp(const cuData3<T>& R, const cuData3<T>& C,
                                  const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                                  eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::interp("
                                    "const cuData3<T>& R, const cuData3<T>& C, "
                                    "const cuData3<T>& dataSrc, cuData3<T>& dataDst, "
                                    "eInter inter)");

    require(dataSrc.getDimension() == 2,
            funcName + ": dataSrc is not 2D");
    require(dataSrc.getSize() == R.getSize(),
            funcName + ": dataSrc and R do not have the same size");
    require(dataSrc.getSize() == C.getSize(),
            funcName + ": dataSrc and C do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());

    cuda_transform_interp(      R.getAddrData(),
                                C.getAddrData(),
                          dataSrc.getAddrData(),
                          dataDst.getAddrData(),
                          R.getNrow(), R.getNcol(),
                          inter);
}

template <typename T>
void cuTransSimilarity<T>::interp(const cuData3<T>& R, const cuData3<T>& C, const cuData3<T>& S,
                                  const cuData3<T>& dataSrc, cuData3<T>& dataDst,
                                  eInter inter)
{
    const std::string   funcName("void cuTransSimilarity<T>::interp("
                                    "const cuData3<T>& R, const cuData3<T>& C, const cuData3<T>& S, "
                                    "const cuData3<T>& dataSrc, cuData3<T>& dataDst, "
                                    "eInter inter)");

    require(dataSrc.getDimension() == 3,
            funcName + ": dataSrc is not 3D");
    require(dataSrc.getSize() == R.getSize(),
            funcName + ": dataSrc and R do not have the same size");
    require(dataSrc.getSize() == C.getSize(),
            funcName + ": dataSrc and C do not have the same size");
    require(dataSrc.getSize() == S.getSize(),
            funcName + ": dataSrc and S do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());

    cuda_transform_interp(      R.getAddrData(),
                                C.getAddrData(),
                                S.getAddrData(),
                          dataSrc.getAddrData(),
                          dataDst.getAddrData(),
                          R.getNrow(), R.getNcol(), R.getNsec(),
                          inter);
}

// instantiation
template class cuTransSimilarity<float >;
template class cuTransSimilarity<double>;

} // namespace gem
