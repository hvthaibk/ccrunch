/***********************************************************************
 *  File:       cArray3.hpp
 *
 *  Purpose:    Header file for a 3D array class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CARRAY3_HPP__
#define __GEM_CARRAY3_HPP__

#include "config.hpp"
#include "macro.hpp"
#include "transform.hpp"
#include "cVector3.hpp"

#include <cstddef>

namespace gem {

#ifdef __GEM_USE_CUDA__
template <typename T> class cuArray3;
#endif

template <typename T>
class cArray3
{
protected:
    T*          _data;
    cSize3      _size;

public:
     cArray3(const cSize3&  size = cSize3(0,0,0));
     cArray3(const cArray3& other);
#ifdef __GEM_USE_CUDA__
     cArray3(const cuArray3<T>& other);
#endif
    ~cArray3();
#ifdef __GEM_USE_STD11__
     cArray3(cArray3&& other);
#endif

    // empty checking
    void    requireEmpty   (const std::string& funcName) const { require(_data == NULL, funcName + ": cArray3 object is not empty"); }
    void    requireNonEmpty(const std::string& funcName) const { require(_data != NULL, funcName + ": cArray3 object is empty");     }

    // dimension and length
    size_t           getNrow        (void) const { return _size[0];                   }
    size_t           getNcol        (void) const { return _size[1];                   }
    size_t           getNsec        (void) const { return _size[2];                   }
    size_t           getNelement    (void) const { return _size[0]*_size[1]*_size[2]; }
    T*               getAddrData    (void) const { return _data;                      }
    cSize3           getSize        (void) const { return _size;                      }
    unsigned int     getDimension   (void) const;

    // operator overloading
    T&              operator[](size_t i)       { assert(i < getNelement());  return _data[i]; }
    T               operator[](size_t i) const { assert(i < getNelement());  return _data[i]; }

    cArray3&        operator= (T value);
    cArray3&        operator= (const cArray3&     other);
#ifdef __GEM_USE_CUDA__
    cArray3&        operator= (const cuArray3<T>& other);
#endif
#ifdef __GEM_USE_STD11__
    cArray3&        operator= (cArray3&& other);
#endif

    const cArray3   operator+ (T value)              const;
          cArray3&  operator+=(T value);
    const cArray3   operator+ (const cArray3& other) const;
          cArray3&  operator+=(const cArray3& other);
    const cArray3   operator- (T value)              const;
          cArray3&  operator-=(T value);
    const cArray3   operator- (const cArray3& other) const;
          cArray3&  operator-=(const cArray3& other);
    const cArray3   operator* (T value)              const;
          cArray3&  operator*=(T value);
    const cArray3   operator* (const cArray3& other) const;
          cArray3&  operator*=(const cArray3& other);
    const cArray3   operator/ (T value)              const;
          cArray3&  operator/=(T value);
    const cArray3   operator/ (const cArray3& other) const;
          cArray3&  operator/=(const cArray3& other);
    const cArray3   operator- (void)                 const;
    bool            operator==(const cArray3& other) const;
    bool            operator!=(const cArray3& other) const;
    bool            operator< (const cArray3& other) const;
    bool            operator<=(const cArray3& other) const;
    bool            operator> (const cArray3& other) const;
    bool            operator>=(const cArray3& other) const;

    // element access
    T       getElement(      size_t  index) const;
    T       getElement(const cSize3& index) const;
    void    setElement(      size_t  index, T value);
    void    setElement(const cSize3& index, T value);

    // reduction
    T       getMax    (void) const;
    T       getMin    (void) const;
    T       getSum    (void) const;
    T       getSum2   (void) const;
    T       getMean   (void) const;
    T       getStd    (void) const;
    size_t  getMaxIndx(void) const;
    size_t  getMinIndx(void) const;

    void    getMax    (const cArray3& other, eReduceNDim dim);
    void    getMin    (const cArray3& other, eReduceNDim dim);
    void    getSum    (const cArray3& other, eReduceNDim dim);
    void    getSum2   (const cArray3& other, eReduceNDim dim);
    void    getMean   (const cArray3& other, eReduceNDim dim);
    void    getStd    (const cArray3& other, eReduceNDim dim);

    // alloc / free / swap
    void    memAlloc            (const cSize3& size);
    void    memAllocZero        (const cSize3& size);
    void    memReAlloc          (const cSize3& size);
    void    memReAllocZero      (const cSize3& size);
    void    memReAllocReduceNDim(const cSize3& size, unsigned int ndim, eReduceNDim dim);
    void    memFree             (void);
    void    memSetZero          (void);
    void    memSetVal           (T value);
    void    memSwap             (cArray3& other);
    void    memReShape          (const cSize3& size);

    // debugging
    void    printData(std::string message = "", int width = 8, int precision = 4, bool scientific = 0) const;
    void    printSize(std::string message = "", int width = 8, int precision = 4) const;

    // visualisation
    void    autocontrast8bit(int cutoff);
    void    autocontrast8bit(T min, T max);
    void    normalize8bit(void);
    void    normalize8bit(const cArray3& other);

    // math
    void    opMathAbs (void);
    void    opMathInv (void);
    void    opMathNorm(void);
    void    opMathSqr (void);
    void    opMathSqrt(void);

    // 1D ray and 2D slide
    void    getRayRow(const cArray3& other, size_t iCol, size_t iSec);
    void    getRayCol(const cArray3& other, size_t iRow, size_t iSec);
    void    getRaySec(const cArray3& other, size_t iRow, size_t iCol);

    void    getSlideRow(const cArray3& other, size_t iRow);
    void    getSlideCol(const cArray3& other, size_t iCol);
    void    getSlideSec(const cArray3& other, size_t iSec);

    // array operations
    void    opCircShift(                      const cOffset3& offset);
    void    opCircShift(const cArray3& other, const cOffset3& offset);

    void    opShift    (                      const cOffset3& offset);
    void    opShift    (const cArray3& other, const cOffset3& offset);

    void    opPad      (                      const cSize3& size, const cSize3& offset);
    void    opPad      (const cArray3& other, const cSize3& size, const cSize3& offset);

    void    opCrop     (                      const cSize3& size, const cSize3& offset);
    void    opCrop     (const cArray3& other, const cSize3& size, const cSize3& offset);

    void    opReplace  (T value, const cSize3& size, const cSize3& offset);
    void    opReplace  (const cArray3& other,        const cSize3& offset);

    void    opPermute  (                      ePermute permute = PERMUTE2D);
    void    opPermute  (const cArray3& other, ePermute permute = PERMUTE2D);

    void    opCopy     (const cArray3& other, const cSize3& size, const cOffset3& offset);

    // indexing
    void    opIndexCol2Row(void);
    void    opIndexRow2Col(void);

    template <typename U>   void    opIndexCol2Row(const cArray3<U>& other);
    template <typename U>   void    opIndexRow2Col(const cArray3<U>& other);

    // type casting
    template <typename U>   void    opTypeCast    (const cArray3<U>& other);

    // conversion from C-type array
    void    arrayImport(const T* const array, size_t nRow, size_t nCol = 1, size_t nSec = 1);
    void    arrayHandle(      T* const array, size_t nRow, size_t nCol = 1, size_t nSec = 1);
    void    arrayExport(      T*       array);
};

template <typename T>
cArray3<T> operator-(T value, const cArray3<T>& array);

// ---------------------------------------------------------------------
// implementation for templated member function
// ---------------------------------------------------------------------

template <typename T> template <typename U>
void cArray3<T>::opIndexCol2Row(const cArray3<U>& other)
{
    const std::string   funcName("void cArray3<T>::opIndexCol2Row(const cArray3<U>& other)");

    other.requireNonEmpty(funcName);

    memReAlloc(other.getSize());

    switch (getDimension()) {
        case 2:
            array_index_col2row(other.getAddrData(), _data,
                                _size[0], _size[1]);
            break;
        case 3:
            array_index_col2row(other.getAddrData(), _data,
                                _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T> template <typename U>
void cArray3<T>::opIndexRow2Col(const cArray3<U>& other)
{
    const std::string   funcName("void cArray3<T>::opIndexRow2Col(const cArray3<U>& other)");

    other.requireNonEmpty(funcName);

    memReAlloc(other.getSize());

    switch (getDimension()) {
        case 2:
            array_index_row2col(other.getAddrData(), _data,
                                _size[0], _size[1]);
            break;
        case 3:
            array_index_row2col(other.getAddrData(), _data,
                                _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T> template <typename U>
void cArray3<T>::opTypeCast(const cArray3<U>& other)
{
    other.requireNonEmpty("void cArray3<T>::opTypeCast(const cArray3<U>& other)");

    memReAlloc(other.getSize());

    if (getNelement() > 0) {
        array_typecast(other.getAddrData(), getNelement(), getAddrData());
    }
}

} // namespace gem

#endif
