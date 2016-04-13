/***********************************************************************
 *  File:       cuArray3.hpp
 *
 *  Purpose:    Header file for a 3D array class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUARRAY3_HPP__
#define __GEM_CUARRAY3_HPP__

#include "cArray3.hpp"

namespace gem {

template <typename T>
class cuArray3
{
protected:
    T*          _data;
    cSize3      _size;

public:
     cuArray3(const cSize3& size = cSize3(0,0,0));
     cuArray3(const cuArray3&   other);
     cuArray3(const cArray3<T>& other);
    ~cuArray3();

    // empty checking
    void    requireEmpty   (const std::string funcName) const { require(_data == NULL, funcName + ": cuArray3 object is not empty"); }
    void    requireNonEmpty(const std::string funcName) const { require(_data != NULL, funcName + ": cuArray3 object is empty");     }

    // dimension and length
    size_t           getNrow        (void) const { return _size[0];                   }
    size_t           getNcol        (void) const { return _size[1];                   }
    size_t           getNsec        (void) const { return _size[2];                   }
    size_t           getNelement    (void) const { return _size[0]*_size[1]*_size[2]; }
    T*               getAddrData    (void) const { return _data;                      }
    cSize3           getSize        (void) const { return _size;                      }
    unsigned int     getDimension   (void) const;

    // operator overloading
    //T&              operator[](size_t i)       { assert(i < getNelement());  return _data[i]; }
    //T               operator[](size_t i) const { assert(i < getNelement());  return _data[i]; }

    cuArray3&        operator= (T value);
    cuArray3&        operator= (const cuArray3&   other);
    cuArray3&        operator= (const cArray3<T>& other);

    const cuArray3   operator+ (T value)               const;
          cuArray3&  operator+=(T value);
    const cuArray3   operator+ (const cuArray3& other) const;
          cuArray3&  operator+=(const cuArray3& other);
    const cuArray3   operator- (T value)               const;
          cuArray3&  operator-=(T value);
    const cuArray3   operator- (const cuArray3& other) const;
          cuArray3&  operator-=(const cuArray3& other);
    const cuArray3   operator* (T value)               const;
          cuArray3&  operator*=(T value);
    const cuArray3   operator* (const cuArray3& other) const;
          cuArray3&  operator*=(const cuArray3& other);
    const cuArray3   operator/ (T value)               const;
          cuArray3&  operator/=(T value);
    const cuArray3   operator/ (const cuArray3& other) const;
          cuArray3&  operator/=(const cuArray3& other);
    const cuArray3   operator- (void)                  const;
    //bool            operator==(const cuArray3& other) const;
    //bool            operator!=(const cuArray3& other) const;
    //bool            operator< (const cuArray3& other) const;
    //bool            operator<=(const cuArray3& other) const;
    //bool            operator> (const cuArray3& other) const;
    //bool            operator>=(const cuArray3& other) const;

    // element access
    //T       getElement(      size_t  index) const;
    //T       getElement(const cSize3& index) const;
    //void    setElement(      size_t  index, T value);
    //void    setElement(const cSize3& index, T value);

    // reduction
    T       getMax    (void) const;
    T       getMin    (void) const;
    T       getSum    (void) const;
    T       getSum2   (void) const;
    T       getMean   (void) const;
    T       getStd    (void) const;
    size_t  getMaxIndx(void) const;
    size_t  getMinIndx(void) const;

    // alloc / free / swap
    void    memAlloc      (const cSize3& size);
    void    memAllocZero  (const cSize3& size);
    void    memReAlloc    (const cSize3& size);
    void    memReAllocZero(const cSize3& size);
    void    memFree       (void);
    void    memSetZero    (void);
    void    memSetVal     (T value);
    void    memSwap       (cuArray3& other);

    // debugging
    void    printSize    (std::string message = "", int width = 8, int precision = 4) const;

    // visualisation
    void    normalize8bit(void);

    // math
    void    opMathAbs (void);
    void    opMathInv (void);
    void    opMathNorm(void);
    void    opMathSqr (void);
    void    opMathSqrt(void);

    // array operations
    void    opCircShift(                       const cOffset3& offset);
    void    opCircShift(const cuArray3& other, const cOffset3& offset);

    void    opShift    (                       const cOffset3& offset);
    void    opShift    (const cuArray3& other, const cOffset3& offset);

    void    opPad      (                       const cSize3& size, const cSize3& offset);
    void    opPad      (const cuArray3& other, const cSize3& size, const cSize3& offset);

    void    opCrop     (                       const cSize3& size, const cSize3& offset);
    void    opCrop     (const cuArray3& other, const cSize3& size, const cSize3& offset);

    void    opReplace  (T value, const cSize3& size, const cSize3& offset);
    void    opReplace  (const cuArray3& other,       const cSize3& offset);

    void    opPermute  (                       ePermute permute = PERMUTE2D);
    void    opPermute  (const cuArray3& other, ePermute permute = PERMUTE2D);

    //void    opCopy     (const cuArray3& other, const cSize3& size, const cOffset3& offset);

    // indexing
    void    opIndexCol2Row(void);
    void    opIndexRow2Col(void);

    template <typename U>   void    opIndexCol2Row(const cuArray3<U>& other);
    template <typename U>   void    opIndexRow2Col(const cuArray3<U>& other);

    // type casting
    template <typename U>   void    opTypeCast    (const cuArray3<U>& other);

    // conversion from C-type array
    void    arrayImport(const T* const array, size_t nRow, size_t nCol = 1, size_t nSec = 1);
    void    arrayHandle(      T* const array, size_t nRow, size_t nCol = 1, size_t nSec = 1);
    void    arrayExport(      T*       array);
};

// ---------------------------------------------------------------------
// implementation for templated member function
// ---------------------------------------------------------------------

template <typename T> template <typename U>
void cuArray3<T>::opIndexCol2Row(const cuArray3<U>& other)
{
    const std::string   funcName("void cuArray3<T>::opIndexCol2Row(const cuArray3<U>& other)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 2:
            cuda_array_index_col2row(other._data, _data,
                                     _size[0], _size[1]);
            break;
        case 3:
            cuda_array_index_col2row(other._data, _data,
                                     _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T> template <typename U>
void cuArray3<T>::opIndexRow2Col(const cuArray3<U>& other)
{
    const std::string   funcName("void cuArray3<T>::opIndexRow2Col(const cuArray3<U>& other)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 2:
            cuda_array_index_row2col(other._data, _data,
                                     _size[0], _size[1]);
            break;
        case 3:
            cuda_array_index_row2col(other._data, _data,
                                     _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T> template <typename U>
void cuArray3<T>::opTypeCast(const cuArray3<U>& other)
{
    other.requireNonEmpty("void cuArray3<T>::opTypeCast(const cuArray3<U>& other)");

    memReAlloc(other.getSize());

    if (getNelement() > 0) {
        cuda_array_typecast(other.getAddrData(), getNelement(), getAddrData());
    }
}

} // namespace gem

#endif
