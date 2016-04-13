/***********************************************************************
 *  File:       cData3XCorr.cpp
 *
 *  Purpose:    Implementation of a correlation-data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cData3XCorr.hpp"

#include "fftSize.hpp"
#include "xcorr.hpp"
#include "cTransDistance.hpp"
#include "cTransSimilarity.hpp"
#include "cTransThreshold.hpp"

namespace gem {

template <typename T>
cData3XCorr<T>::cData3XCorr(const cSize3& size)
    : cData3<T>(size), _id(0), _fileName(""), _mskArea(0)
{}

template <typename T>
cData3XCorr<T>::cData3XCorr(const cData3XCorr& other)
    : cData3<T>(other), _id(0), _fileName(""), _mskArea(0)
{}

template <typename T>
void cData3XCorr<T>::opFFTSizeChange(cSize3& sizeOri, const cSize3& sizeMin)
{
    const std::string   funcName("void cData3XCorr<T>::opFFTSizeChange(cSize3& sizeOri, const cSize3& sizeMin)");

    cArray3<T>::requireNonEmpty(funcName);

    cSize3    sizeFFT;

    switch (cArray3<T>::getDimension()) {
        case 2:
            sizeFFT = cSize3(fftGoodSize[std::max(cArray3<T>::_size[0],sizeMin[0])],
                                       fftGoodSize[std::max(cArray3<T>::_size[1],sizeMin[1])],
                                       1);

            if (cArray3<T>::_size != sizeFFT) {
/*#ifndef NDEBUG
                std::cout << "FFT size: change  size from "
                          << cArray3<T>::_size[0] << "x"
                          << cArray3<T>::_size[1] << " to "
                          << sizeFFT[0] << "x"
                          << sizeFFT[1] << std::endl;
#endif*/

                sizeOri = cSize3(cArray3<T>::_size[0],
                                           cArray3<T>::_size[1],
                                           1);

                cArray3<T>::opPad(sizeFFT, cSize3(0,0,0));
            }
            else {
                sizeOri = cSize3(cArray3<T>::_size[0],
                                           cArray3<T>::_size[1],
                                           1);
            }

            break;

        case 3:
            sizeFFT = cSize3(fftGoodSize3D[std::max(cArray3<T>::_size[0],sizeMin[0])],
                                       fftGoodSize3D[std::max(cArray3<T>::_size[1],sizeMin[1])],
                                       fftGoodSize3D[std::max(cArray3<T>::_size[2],sizeMin[2])]);

            if (cArray3<T>::_size != sizeFFT) {
/*#ifndef NDEBUG
                std::cout << "FFT size: change  size from "
                          << cArray3<T>::_size[0] << "x"
                          << cArray3<T>::_size[1] << "x"
                          << cArray3<T>::_size[2] << " to "
                          << sizeFFT[0] << "x"
                          << sizeFFT[1] << "x"
                          << sizeFFT[2] << std::endl;
#endif*/

                sizeOri = cSize3(cArray3<T>::_size[0],
                                           cArray3<T>::_size[1],
                                           cArray3<T>::_size[2]);

                cArray3<T>::opPad(sizeFFT, cSize3(0,0,0));
            }
            else {
                sizeOri = cSize3(cArray3<T>::_size[0],
                                           cArray3<T>::_size[1],
                                           cArray3<T>::_size[2]);
            }

            break;

        default:
            ERROR(funcName, "unsupported dimension");
    }

    require(sizeOri > (size_t) 0, funcName + "sizeOri is zero");
}

template <typename T>
void cData3XCorr<T>::opFFTSizeRestore(const cSize3& sizeOri)
{
    const std::string   funcName("void cData3XCorr<T>::opFFTSizeRestore(const cSize3& sizeOri)");

    cArray3<T>::requireNonEmpty(funcName);

    require(sizeOri > (size_t) 0, funcName + "sizeOri is zero");

    switch (cArray3<T>::getDimension()) {
        case 2:
            if (sizeOri[0]*sizeOri[1] > 0) {
/*#ifndef NDEBUG
                std::cout << "FFT size: restore size from "
                          << cArray3<T>::_size[0] << "x"
                          << cArray3<T>::_size[1] << " to "
                          << sizeOri[0] << "x"
                          << sizeOri[1] << std::endl << std::endl;
#endif*/

                if (cArray3<T>::_size != sizeOri) {
                    cArray3<T>::opCrop(sizeOri, cSize3(0,0,0));
                }
                else {
                    //WARNING(funcName, "no need to restore");
                }
            }
            break;
        case 3:
            if (sizeOri[0]*sizeOri[1]*sizeOri[2] > 0) {
/*#ifndef NDEBUG
                std::cout << "FFT size: restore size from "
                          << cArray3<T>::_size[0] << "x"
                          << cArray3<T>::_size[1] << "x"
                          << cArray3<T>::_size[2] << " to "
                          << sizeOri[0] << "x"
                          << sizeOri[1] << "x"
                          << sizeOri[2] << std::endl << std::endl;
#endif*/

                if (cArray3<T>::_size != sizeOri) {
                    cArray3<T>::opCrop(sizeOri, cSize3(0,0,0));
                }
                else {
                    //WARNING(funcName, "no need to restore");
                }
            }
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cData3XCorr<T>::opNormalizeTplECC(bool useMask)
{
    const std::string   funcName("void cData3XCorr<T>::opNormalizeTplECC(bool useMask)");

    cArray3<T>::requireNonEmpty(funcName);

    if (useMask) {
        require(cArray3<T>::getSize() == _mask.getSize(),
                funcName + ": data and mask do not have the same size");

        xcorrNormTemplateECC(cArray3<T>::_data,
                             cArray3<T>::getNelement(),
                             _mask.getAddrData());
    }
    else {
        xcorrNormTemplateECC(cArray3<T>::_data,
                             cArray3<T>::getNelement());
    }
}

template <typename T>
void cData3XCorr<T>::opNormalizeTplNCC(bool useMask)
{
    const std::string   funcName("void cData3XCorr<T>::opNormalizeTplNCC(bool useMask)");

    cArray3<T>::requireNonEmpty(funcName);

    if (useMask) {
        require(cArray3<T>::getSize() == _mask.getSize(),
                funcName + ": data and mask do not have the same size");

        _mskArea = xcorrNormTemplateNCC(cArray3<T>::_data,
                                        cArray3<T>::getNelement(),
                                        _mask.getAddrData());
    }
    else {
        _mskArea = xcorrNormTemplateNCC(cArray3<T>::_data,
                                        cArray3<T>::getNelement());

    }
}

template <typename T>
void cData3XCorr<T>::opNormalizeTplWCC1(bool useMask)
{
    const std::string   funcName("void cData3XCorr<T>::opNormalizeTplWCC1(bool useMask)");

    cArray3<T>::requireNonEmpty(funcName);
    require(cArray3<T>::getSize() == _wghtTpl.getSize(),
            funcName + ": data and weight do not have the same size");

    if (useMask) {
        require(cArray3<T>::getSize() == _mask.getSize(),
                funcName + ": data and mask do not have the same size");

        xcorrNormTemplateWCC1(cArray3<T>::_data,
                              cArray3<T>::getNelement(),
                              _mask.getAddrData(),
                              _wghtTpl.getAddrData());
    }
    else {
        xcorrNormTemplateWCC1(cArray3<T>::_data,
                              cArray3<T>::getNelement(),
                              _wghtTpl.getAddrData());
    }
}

template <typename T>
void cData3XCorr<T>::opNormalizeTplWCC2(bool useMask)
{
    const std::string   funcName("void cData3XCorr<T>::opNormalizeTplWCC2(bool useMask)");

    cArray3<T>::requireNonEmpty(funcName);
    require(cArray3<T>::getSize() == _wghtTpl.getSize() &&
            cArray3<T>::getSize() == _wghtRef.getSize(),
            funcName + ": data and weights do not have the same size");

    if (useMask) {
        require(cArray3<T>::getSize() == _mask.getSize(),
                funcName + ": data and mask do not have the same size");

        xcorrNormTemplateWCC2(cArray3<T>::_data,
                              cArray3<T>::getNelement(),
                              _mask.getAddrData(),
                              _wghtTpl.getAddrData(),
                              _wghtRef.getAddrData());
    }
    else {
        xcorrNormTemplateWCC2(cArray3<T>::_data,
                              cArray3<T>::getNelement(),
                              _wghtTpl.getAddrData(),
                              _wghtRef.getAddrData());
    }
}

template <typename T>
void cData3XCorr<T>::opLimitCorrValue(T valMin, T valMax)
{
    const std::string   funcName("void cData3XCorr<T>::opLimitCorrValue(T valMin, T valMax)");

    cArray3<T>::requireNonEmpty(funcName);

    #pragma omp parallel for
    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        if (cArray3<T>::_data[i] < valMin || cArray3<T>::_data[i] > valMax) {
            cArray3<T>::_data[i] = 0;
        }
    }
}

template <typename T>
void cData3XCorr<T>::opClearBorderMargin(size_t marginSize)
{
    const std::string   funcName("void cData3XCorr<T>::opClearBorderMargin(size_t marginSize)");

    cArray3<T>::requireNonEmpty(funcName);

    size_t  nrow = cArray3<T>::getNrow();
    size_t  ncol = cArray3<T>::getNrow();
    size_t  irow, icol;

    require(marginSize < nrow/2 && marginSize < ncol/2,
            funcName + ": margin size is too large");

    #pragma omp parallel for private(irow, icol)
    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        irow = i / ncol;
        icol = i % ncol;

        if (irow < marginSize || irow >= nrow-marginSize ||
            icol < marginSize || icol >= ncol-marginSize) {
                cArray3<T>::_data[i] = 0;
        }

    }
}

template <typename T>
void cData3XCorr<T>::opReplaceUsingMask(const cData3XCorr& mask, T value, const cVector3<ptrdiff_t>& offset, bool maskInside)
{
    const std::string   funcName("void cData3XCorr<T>::opReplaceUsingMask(const cData3XCorr& mask, T value, const cVector3<ptrdiff_t>& offset, bool maskInside)");

    cArray3<T>::requireNonEmpty(funcName);
           mask.requireNonEmpty(funcName);
    require(cArray3<T>::getDimension() == mask.getDimension(),
            funcName + ": data and mask do not have the same dimension");

    switch (cArray3<T>::getDimension()) {
        case 1:
            array_mask_replace(cArray3<T>::_data, cArray3<T>::_size[0],
                                      mask._data,        mask._size[0],
                                           value,            offset[0],
                                      maskInside);
            break;
        case 2:
            array_mask_replace(cArray3<T>::_data, cArray3<T>::_size[0], cArray3<T>::_size[1],
                                      mask._data,        mask._size[0],        mask._size[1],
                                           value,            offset[0],            offset[1],
                                      maskInside);
            break;
        case 3:
            array_mask_replace(cArray3<T>::_data, cArray3<T>::_size[0], cArray3<T>::_size[1], cArray3<T>::_size[2],
                                      mask._data,        mask._size[0],        mask._size[1],        mask._size[2],
                                           value,            offset[0],            offset[1],            offset[2],
                                      maskInside);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cData3XCorr<T>::opCopyUsingMask(const cData3XCorr& mask, cData3XCorr& output, const cVector3<ptrdiff_t>& offset)
{
    const std::string   funcName("void cData3XCorr<T>::opCopyUsingMask(const cData3XCorr& mask, cData3XCorr& output, const cVector3<ptrdiff_t>& offset)");

    cArray3<T>::requireNonEmpty(funcName);
           mask.requireNonEmpty(funcName);
         output.requireNonEmpty(funcName);
    require(mask.getSize() == output.getSize(),
            funcName + ": mask and output do not have the same size");

    switch (cArray3<T>::getDimension()) {
        case 1:
            array_mask_copy(cArray3<T>::_data, cArray3<T>::_size[0],
                                   mask._data,        mask._size[0],
                                 output._data,            offset[0]);
            break;
        case 2:
            array_mask_copy(cArray3<T>::_data, cArray3<T>::_size[0], cArray3<T>::_size[1],
                                   mask._data,        mask._size[0],        mask._size[1],
                                 output._data,            offset[0],            offset[1]);
            break;
        case 3:
            array_mask_copy(cArray3<T>::_data, cArray3<T>::_size[0], cArray3<T>::_size[1], cArray3<T>::_size[2],
                                   mask._data,        mask._size[0],        mask._size[1],        mask._size[2],
                                 output._data,            offset[0],            offset[1],            offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cData3XCorr<T>::maskCreateByCopy(const cData3XCorr& object)
{
    const std::string   funcName("void cData3XCorr<T>::maskCreateByCopy(const cData3XCorr& object)");

    object.requireNonEmpty(funcName);

    _mask = object;
    cTransThreshold<T>().threshold(_mask, 0.5f);
}

template <typename T>
void cData3XCorr<T>::maskCreateByRotation(const cData3XCorr& object, const cVector3<T>& angle, eInter inter)
{
    const std::string   funcName("void cData3XCorr<T>::maskCreateByRotation(const cData3XCorr& object, const cVector3<T>& angle, eInter inter)");

    object.requireNonEmpty(funcName);

    cTransSimilarity<T>().rotate(object, _mask, angle, inter);
    cTransThreshold<T>().threshold(_mask, 0.5f);
}

template <typename T>
void cData3XCorr<T>::maskCreateByOnes(void)
{
    const std::string   funcName("void cData3XCorr<T>::maskCreateByOnes(void)");

    cArray3<T>::requireNonEmpty(funcName);

    // cannot use the cArray3<T>::operator=(1) here because
    // the compiler may misunderstand it as cData3::operator=(const cuData3<T>& other)
    _mask.memReAlloc(cArray3<T>::getSize());
    _mask.memSetVal(1);
}

template <typename T>
void cData3XCorr<T>::maskCreateByThresholding(T value)
{
    const std::string   funcName("void cData3XCorr<T>::maskCreateByThresholding(T value)");

    cArray3<T>::requireNonEmpty(funcName);

    _mask = *this;
    cTransThreshold<T>().threshold(_mask, value);
}

template <typename T>
void cData3XCorr<T>::maskCheckValue(void)
{
    const std::string   funcName("void cData3XCorr<T>::maskCheckvalue(void)");

    cArray3<T>::requireNonEmpty(funcName);
    _mask.requireNonEmpty(funcName);

    for (size_t i = 0; i < cArray3<T>::getNelement(); i++) {
        require(_mask[i] == 0 || _mask[i] == 1, funcName + ": mask data is not 0 / 1");
    }
}

template <typename T>
void cData3XCorr<T>::maskGetIndexLimit(size_t &indRowMin, size_t &indRowMax,
                                       size_t &indColMin, size_t &indColMax)
{
    const std::string   funcName("void cData3XCorr<T>::maskGetIndexLimit(size_t &indRowMin, size_t &indRowMax, size_t &indColMin, size_t &indColMax)");

    cArray3<T>::requireNonEmpty(funcName);
    require(cArray3<T>::getDimension() == 2, funcName + ": only support 2D mask");

    size_t    *arrayInd = NULL, nFound;
    size_t    *arrayIndRow = NULL;
    size_t    *arrayIndCol = NULL;

    array_new(arrayInd, cArray3<T>::getNelement());

    array_find(cArray3<T>::_data, cArray3<T>::getNelement(), arrayInd, nFound);

    array_new(arrayIndRow, nFound);
    array_new(arrayIndCol, nFound);

    array_index_ind2sub(arrayInd, arrayIndRow, arrayIndCol,
                        nFound, cArray3<T>::_size[0], cArray3<T>::_size[1]);

    indRowMin = array_reduce_min(arrayIndRow, nFound);
    indRowMax = array_reduce_max(arrayIndRow, nFound);
    indColMin = array_reduce_min(arrayIndCol, nFound);
    indColMax = array_reduce_max(arrayIndCol, nFound);

    array_delete(arrayInd);
    array_delete(arrayIndRow);
    array_delete(arrayIndCol);
}

template <typename T>
void cData3XCorr<T>::maskShrinkToLimit(size_t indRowMin, size_t indRowMax,
                                       size_t indColMin, size_t indColMax)
{
    const std::string   funcName("void cData3XCorr<T>::maskShrinkToLimit(size_t indRowMin, size_t indRowMax, size_t indColMin, size_t indColMax)");

    cArray3<T>::requireNonEmpty(funcName);
    require(cArray3<T>::getDimension() == 2, funcName + ": only support 2D mask");

    size_t    iRow, iCol, indx;

    for (iRow = 0, indx = 0; iRow < cArray3<T>::_size[0]; iRow++) {
        for (iCol = 0; iCol < cArray3<T>::_size[1]; iCol++, indx++) {
            if (iRow < indRowMin || iRow > indRowMax ||
                iCol < indColMin || iCol > indColMax) {
                cArray3<T>::_data[indx] = 0;
            }
        }
    }
}

template <typename T>
void cData3XCorr<T>::wghtTplCreate(unsigned int mode)
{
    const std::string   funcName("void cData3XCorr<T>::wghtTplCreate(unsigned int mode)");

    cArray3<T>::requireNonEmpty(funcName);

    cData3<T>   maskInv;

    switch (mode) {
        case 0:     // constant value, WCC->NCC
            _wghtTpl.memReAlloc(cArray3<T>::getSize());
            _wghtTpl.memSetVal(1);

            break;

        case 1:     // distance transform
            maskInv = _mask;
            maskInv.opMathInv();
            maskInv += 1;

            cTransDistance<T>().distanceNormal(maskInv, _wghtTpl);

            break;

        case 2:     // template data
            _wghtTpl = *this;
            _wghtTpl.opMathAbs();

            break;

        case 3:     // Gaussian

            break;

        default:
            ERROR(funcName, "unsupported mode for weight definition");

    }

    // Gaussian

    // mode 3: proportional to the inversed absolute value of template
    /*for (size_t i = 0; i < _wghtTpl.getNelement(); i++) {
        if (_mask[i] > 0) {
            _wghtTpl[i] = std::abs(1/_wghtTpl[i]);
        }
    }
    _wghtTpl /= _wghtTpl.getMax();*/

    // mode 4: simulate the Laplacian behaviour
    /*_wghtTpl.filterCreateLaplacian();
    _wghtTpl.opMathAbs();
    _wghtTpl /= _wghtTpl.getMax();*/

    // mode 5:
    /*_wghtTpl -= _wghtTpl.getMin();
    _wghtTpl.opMathSqr();
    _wghtTpl /= _wghtTpl.getMax();*/

    // mode 6:
    /*_wghtTpl.filterLaplacian(3.0/16, 0, 3);
    _wghtTpl.opMathAbs();
    _wghtTpl.opMathSqr();
    _wghtTpl /= _wghtTpl.getMax();*/
}

template <typename T>
void cData3XCorr<T>::wghtRefCreateUsingWghtTpl(void)
{
    const std::string   funcName("void cData3XCorr<T>::wghtRefCreateUsingWghtTpl(void)");

    _wghtTpl.requireNonEmpty(funcName);

    _wghtRef = _wghtTpl;
}

// instantiation
template class cData3XCorr<float >;
template class cData3XCorr<double>;

} // namespace gem
