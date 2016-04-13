/***********************************************************************
 *  File:       cTransSimilarity.cpp
 *
 *  Purpose:    Implementation of a similarity transformation class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTransSimilarity.hpp"

namespace gem {

template <typename T>
void cTransSimilarity<T>::bin(cData3<T>& dataSrc, size_t factor)
{
    const std::string   funcName("void cTransSimilarity<T>::bin("
                                    "cData3<T>& dataSrc, size_t factor)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>    other(cVector3<size_t>(transform_binSize(dataSrc.getNrow(), factor),
                                        transform_binSize(dataSrc.getNcol(), factor),
                                        transform_binSize(dataSrc.getNsec(), factor)));

    switch (dataSrc.getDimension()) {
        case 1:
            transform_bin(dataSrc.getAddrData(), other.getAddrData(),
                          dataSrc.getNrow(),
                          factor);
            break;
        case 2:
            transform_bin(dataSrc.getAddrData(), other.getAddrData(),
                          dataSrc.getNrow(), dataSrc.getNcol(),
                          factor);
            break;
        case 3:
            transform_bin(dataSrc.getAddrData(), other.getAddrData(),
                          dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                          factor);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(other);
}

template <typename T>
void cTransSimilarity<T>::bin(const cData3<T>& dataSrc, cData3<T>& dataDst, size_t factor)
{
    const std::string   funcName("void cTransSimilarity<T>::bin("
                                    "const cData3<T>& dataSrc, "
                                    "cData3<T>& dataDst, size_t factor)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(cVector3<size_t>(transform_binSize(dataSrc.getNrow(), factor),
                                        transform_binSize(dataSrc.getNcol(), factor),
                                        transform_binSize(dataSrc.getNsec(), factor)));

    switch (dataSrc.getDimension()) {
        case 1:
            transform_bin(dataSrc.getAddrData(), dataDst.getAddrData(),
                          dataSrc.getNrow(),
                          factor);
            break;
        case 2:
            transform_bin(dataSrc.getAddrData(), dataDst.getAddrData(),
                          dataSrc.getNrow(), dataSrc.getNcol(),
                          factor);
            break;
        case 3:
            transform_bin(dataSrc.getAddrData(), dataDst.getAddrData(),
                          dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                          factor);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cTransSimilarity<T>::extendRotate(cData3<T>& dataSrc, eExpRot exprot)
{
    const std::string   funcName("void cTransSimilarity<T>::extendRotate("
                                    "cData3<T>& dataSrc, eExpRot exprot)");

    dataSrc.requireNonEmpty(funcName);

    T*          dataTmp = NULL;
    size_t      nRowTmp, nColTmp, nSecTmp = 1;

    switch (dataSrc.getDimension()) {
        case 2:
            switch (exprot) {
                case EXP_ROT_VALID:
                    transform_expand_rotate_valid(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                        dataTmp,               nRowTmp,           nColTmp );
                    break;
                case EXP_ROT_FULL:
                    transform_expand_rotate_full(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                        dataTmp,               nRowTmp,           nColTmp );
                    break;
            }
            break;
        case 3:
            switch (exprot) {
                case EXP_ROT_VALID:
                    transform_expand_rotate_valid(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                        dataTmp,               nRowTmp,           nColTmp,           nSecTmp );
                    break;
                case EXP_ROT_FULL:
                    transform_expand_rotate_full(
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
void cTransSimilarity<T>::extendRotate(const cData3<T>& dataSrc, cData3<T>& dataDst, eExpRot exprot)
{
    const std::string   funcName("void cTransSimilarity<T>::extendRotate("
                                    "const cData3<T>& dataSrc, "
                                    "cData3<T>& dataDst, eExpRot exprot)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    T*          dataTmp = NULL;
    size_t      nRowTmp, nColTmp, nSecTmp = 1;

    switch (dataSrc.getDimension()) {
        case 2:
            switch (exprot) {
                case EXP_ROT_VALID:
                    transform_expand_rotate_valid(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                        dataTmp,               nRowTmp,           nColTmp );
                    break;
                case EXP_ROT_FULL:
                    transform_expand_rotate_full(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(),
                        dataTmp,               nRowTmp,           nColTmp );
                    break;
            }
            break;
        case 3:
            switch (exprot) {
                case EXP_ROT_VALID:
                    transform_expand_rotate_valid(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                        dataTmp,               nRowTmp,           nColTmp,           nSecTmp );
                    break;
                case EXP_ROT_FULL:
                    transform_expand_rotate_full(
                        dataSrc.getAddrData(), dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                        dataTmp,               nRowTmp,           nColTmp,           nSecTmp );
                    break;
            }
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataDst.memFree();
    dataDst.arrayHandle(dataTmp, nRowTmp, nColTmp, nSecTmp);
}

template <typename T>
void cTransSimilarity<T>::translate(cData3<T>& dataSrc, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::translate(cData3<T>& dataSrc, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                dataSrc.getNrow(),
                                offset[0],
                                inter);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                dataSrc.getNrow(), dataSrc.getNcol(),
                                offset[0],         offset[1],
                                inter);
            break;
        case 3:
            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
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
void cTransSimilarity<T>::translate(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::translate(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                dataDst.getNrow(),
                                offset[0],
                                inter);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                dataDst.getNrow(), dataDst.getNcol(),
                                offset[0],         offset[1],
                                inter);
            break;
        case 3:
            transform_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                dataDst.getNrow(), dataDst.getNcol(), dataDst.getNsec(),
                                offset[0],         offset[1],         offset[2],
                                inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cTransSimilarity<T>::rotate(cData3<T>& dataSrc, const cVector3<T>& angle, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::rotate(cData3<T>& dataSrc, const cVector3<T>& angle, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 2:
            require(angle[1] == 0, funcName + ": unused angle != 0");
            require(angle[2] == 0, funcName + ": unused angle != 0");

            transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                             dataSrc.getNrow(), dataSrc.getNcol(),
                             angle[0],
                             inter);
            break;
        case 3:
            transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
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
void cTransSimilarity<T>::rotate(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& angle, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::rotate(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& angle, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 2:
            require(angle[1] == 0, funcName + ": unused angle != 0");
            require(angle[2] == 0, funcName + ": unused angle != 0");

            transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                             dataDst.getNrow(), dataDst.getNcol(),
                             angle[0],
                             inter);
            break;
        case 3:
            transform_rotate(dataSrc.getAddrData(), dataDst.getAddrData(),
                             dataDst.getNrow(), dataDst.getNcol(), dataDst.getNsec(),
                             angle[0], angle[1], angle[2],
                             inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cTransSimilarity<T>::scale(cData3<T>& dataSrc, const T factor, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::scale(cData3<T>& dataSrc, const T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(factor > 0, funcName + ": negative scaling factor");

    cData3<T>    dataDst(cSize3(transform_scaleSize(dataSrc.getNrow(), factor),
                                transform_scaleSize(dataSrc.getNcol(), factor),
                                transform_scaleSize(dataSrc.getNsec(), factor)));

    switch (dataSrc.getDimension()) {
        case 1:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(),
                            factor,
                            inter);
            break;
        case 2:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(),
                            factor,
                            inter);
            break;
        case 3:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
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
void cTransSimilarity<T>::scale(const cData3<T>& dataSrc, cData3<T>& dataDst, T factor, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::scale(const cData3<T>& dataSrc, cData3<T>& dataDst, T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(factor > 0, funcName + ": negative scaling factor");

    dataDst.memReAlloc(cSize3(transform_scaleSize(dataSrc.getNrow(), factor),
                              transform_scaleSize(dataSrc.getNcol(), factor),
                              transform_scaleSize(dataSrc.getNsec(), factor)));

    switch (dataDst.getDimension()) {
        case 1:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(),
                            factor,
                            inter);
            break;
        case 2:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(),
                            factor,
                            inter);
            break;
        case 3:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                            factor,
                            inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cTransSimilarity<T>::scale(cData3<T>& dataSrc, const cVector3<T>& origin, T factor, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::scale(cData3<T>& dataSrc, const cVector3<T>& origin, T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(factor > 0, funcName + ": negative scaling factor");

    cData3<T>    dataDst(cSize3(transform_scaleSize(dataSrc.getNrow(), origin[0], factor),
                                transform_scaleSize(dataSrc.getNcol(), origin[1], factor),
                                transform_scaleSize(dataSrc.getNsec(), origin[2], factor)));

    switch (dataSrc.getDimension()) {
        case 1:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(),
                            origin[0],
                            factor,
                            inter);
            break;
        case 2:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(),
                            origin[0],         origin[1],
                            factor,
                            inter);
            break;
        case 3:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                            origin[0],         origin[1],         origin[2],
                            factor,
                            inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    dataSrc.memSwap(dataDst);
}

template <typename T>
void cTransSimilarity<T>::scale(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& origin, T factor, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::scale(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& origin, T factor, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");
    require(factor > 0, funcName + ": negative scaling factor");

    dataDst.memReAlloc(cSize3(transform_scaleSize(dataSrc.getNrow(), origin[0], factor),
                              transform_scaleSize(dataSrc.getNcol(), origin[1], factor),
                              transform_scaleSize(dataSrc.getNsec(), origin[2], factor)));

    switch (dataDst.getDimension()) {
        case 1:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(),
                            origin[0],
                            factor,
                            inter);
            break;
        case 2:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(),
                            origin[0],         origin[1],
                            factor,
                            inter);
            break;
        case 3:
            transform_scale(dataSrc.getAddrData(), dataDst.getAddrData(),
                            dataSrc.getNrow(), dataSrc.getNcol(), dataSrc.getNsec(),
                            origin[0],         origin[1],         origin[2],
                            factor,
                            inter);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cTransSimilarity<T>::rotTrans(cData3<T>& dataSrc, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::rotTrans(cData3<T>& dataSrc, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);

    cData3<T>    dataDst(dataSrc.getSize());

    switch (dataSrc.getDimension()) {
        case 2:
            require(angle[1]  == 0, funcName + ": unused angle != 0");
            require(angle[2]  == 0, funcName + ": unused angle != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                       dataSrc.getNrow(), dataSrc.getNcol(),
                                        angle[0],
                                       offset[0],         offset[1],
                                       inter);
            break;
        case 3:
            transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
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
void cTransSimilarity<T>::rotTrans(const cData3<T>& dataSrc, cData3<T>& dataDst,
        const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::rotTrans(const cData3<T>& dataSrc, cData3<T>& dataDst, const cVector3<T>& angle, const cVector3<T>& offset, eInter inter)");

    dataSrc.requireNonEmpty(funcName);
    require(&dataSrc != &dataDst, funcName + ": two objects are the same");

    dataDst.memReAlloc(dataSrc.getSize());

    switch (dataDst.getDimension()) {
        case 2:
            require(angle[1]  == 0, funcName + ": unused angle != 0");
            require(angle[2]  == 0, funcName + ": unused angle != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
                                       dataDst.getNrow(), dataDst.getNcol(),
                                        angle[0],
                                       offset[0],         offset[1],
                                       inter);
            break;
        case 3:
            transform_rotate_translate(dataSrc.getAddrData(), dataDst.getAddrData(),
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
void cTransSimilarity<T>::interp(const cData3<T>& R,
                                 const cData3<T>& dataSrc, cData3<T>& dataDst,
                                 eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::interp("
                                    "const cData3<T>& R, "
                                    "const cData3<T>& dataSrc, cData3<T>& dataDst, "
                                    "eInter inter)");

    require(dataSrc.getDimension() == 1,
            funcName + ": dataSrc is not 1D");
    require(dataSrc.getSize() == R.getSize(),
            funcName + ": dataSrc and R do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());

    transform_interp(      R.getAddrData(),
                     dataSrc.getAddrData(),
                     dataDst.getAddrData(),
                     R.getNrow(),
                     inter);
}

template <typename T>
void cTransSimilarity<T>::interp(const cData3<T>& R, const cData3<T>& C,
                                 const cData3<T>& dataSrc, cData3<T>& dataDst,
                                 eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::interp("
                                    "const cData3<T>& R, const cData3<T>& C, "
                                    "const cData3<T>& dataSrc, cData3<T>& dataDst, "
                                    "eInter inter)");

    require(dataSrc.getDimension() == 2,
            funcName + ": dataSrc is not 2D");
    require(dataSrc.getSize() == R.getSize(),
            funcName + ": dataSrc and R do not have the same size");
    require(dataSrc.getSize() == C.getSize(),
            funcName + ": dataSrc and C do not have the same size");

    dataDst.memReAlloc(dataSrc.getSize());

    transform_interp(      R.getAddrData(),
                           C.getAddrData(),
                     dataSrc.getAddrData(),
                     dataDst.getAddrData(),
                     R.getNrow(), R.getNcol(),
                     inter);
}

template <typename T>
void cTransSimilarity<T>::interp(const cData3<T>& R, const cData3<T>& C, const cData3<T>& S,
                                 const cData3<T>& dataSrc, cData3<T>& dataDst,
                                 eInter inter)
{
    const std::string   funcName("void cTransSimilarity<T>::interp("
                                    "const cData3<T>& R, const cData3<T>& C, const cData3<T>& S, "
                                    "const cData3<T>& dataSrc, cData3<T>& dataDst, "
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

    transform_interp(      R.getAddrData(),
                           C.getAddrData(),
                           S.getAddrData(),
                     dataSrc.getAddrData(),
                     dataDst.getAddrData(),
                     R.getNrow(), R.getNcol(), R.getNsec(),
                     inter);
}

// instantiation
template class cTransSimilarity<float >;
template class cTransSimilarity<double>;

} // namespace gem
