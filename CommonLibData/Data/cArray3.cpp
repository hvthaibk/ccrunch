/***********************************************************************
 *  File:       cArray3.cpp
 *
 *  Purpose:    Implementation of a 3D array class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cArray3.hpp"

#include "config.hpp"
#include "macro.hpp"
#include "array.hpp"

#ifdef __GEM_USE_CUDA__
#include "cuArray3.hpp"
#include "array.cuh"
#endif  // __GEM_USE_CUDA__

namespace gem {

template <typename T>
cArray3<T>::cArray3(const cSize3& size)
{
    _size = size;
    _data = NULL;

    if (getNelement() > 0) {
        array_new_zero(_data, getNelement());
    }
    else {
#ifndef NDEBUG
        //WARNING("cArray3::cArray3", "an empty array is requested");
#endif
    }
}

template <typename T>
cArray3<T>::cArray3(const cArray3& other)
{
    _size = other._size;
    _data = NULL;

    if (getNelement() > 0) {
        array_new(_data, getNelement());
        array_memcpy(_data, other._data, getNelement());
    }
    else {
#ifndef NDEBUG
        WARNING("cArray3::cArray3", "an empty array is requested");
#endif
    }
}

#ifdef __GEM_USE_CUDA__
template <typename T>
cArray3<T>::cArray3(const cuArray3<T>& other)
{
    _size = other.getSize();
    _data = NULL;

    if (getNelement() > 0) {
        array_new(_data, getNelement());
        cuda_array_memcpy_d2h(_data, other.getAddrData(), getNelement());
    }
    else {
#ifndef NDEBUG
        WARNING("cArray3::cArray3", "an empty array is requested");
#endif
    }
}
#endif  // __GEM_USE_CUDA__

template <typename T>
cArray3<T>::~cArray3()
{
    memFree();
}

#ifdef __GEM_USE_STD11__
template <typename T>
cArray3<T>::cArray3(cArray3&& other)
{
    _size = other._size;
    _data = NULL;

    if (getNelement() > 0) {
        _data = other._data;

        other._size = 0;
        other._data = NULL;
    }
    else {
#ifndef NDEBUG
        WARNING("cArray3::cArray3", "an empty array is requested");
#endif
    }
}
#endif  // __GEM_USE_STD11__

template <typename T>
unsigned int cArray3<T>::getDimension(void) const
{
    unsigned int    ndim = (_size[0]>1)+(_size[1]>1)+(_size[2]>1);

    if (_size[1] > 1 && ndim < 2)  ndim = 2;
    if (_size[2] > 1 && ndim < 3)  ndim = 3;

    return ndim;
}

template <typename T>
cArray3<T>& cArray3<T>::operator=(T value)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator=(T value)");

    requireNonEmpty(funcName);

    array_setval(_data, value, getNelement());

    return *this;
}

template <typename T>
cArray3<T>& cArray3<T>::operator=(const cArray3& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator=(const cArray3& other)");

    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    if (getNelement() > 0) {
        array_memcpy(_data, other._data, getNelement());
    }

    return *this;
}

#ifdef __GEM_USE_CUDA__
template <typename T>
cArray3<T>& cArray3<T>::operator=(const cuArray3<T>& other)
{
    memReAlloc(other.getSize());

    if (getNelement() > 0) {
        cuda_array_memcpy_d2h(_data, other.getAddrData(), getNelement());
    }

    return *this;
}
#endif  // __GEM_USE_CUDA__

#ifdef __GEM_USE_STD11__
template <typename T>
cArray3<T>& cArray3<T>::operator=(cArray3&& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator=(cArray3&& other)");

    require(this != &other, funcName + ": two objects are the same");

    _size = other._size;
    _data = other._data;

    other._size = 0;
    other._data = NULL;

    return *this;
}
#endif  // __GEM_USE_STD11__

template <typename T>
const cArray3<T> cArray3<T>::operator+(T value) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator+(T value) const");

    requireNonEmpty(funcName);

    cArray3<T>    object(_size);

    array_math_add(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator+=(T value)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator+=(T value)");

    requireNonEmpty(funcName);

    array_math_add(_data, value, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator+(const cArray3& other) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator+(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    cArray3<T>    object(_size);

    array_math_add(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator+=(const cArray3& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator+=(const cArray3& other)");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    array_math_add(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator-(T value) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator-(T value) const");

    requireNonEmpty(funcName);

    cArray3<T>    object(_size);

    array_math_sub(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator-=(T value)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator-=(T value)");

    requireNonEmpty(funcName);

    array_math_sub(_data, value, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator-(const cArray3& other) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator-(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    cArray3<T>    object(_size);

    array_math_sub(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator-=(const cArray3& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator-=(const cArray3& other)");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    array_math_sub(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator*(T value) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator*(T value) const");

    requireNonEmpty(funcName);

    cArray3<T>    object(_size);

    array_math_mul(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator*=(T value)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator*=(T value)");

    requireNonEmpty(funcName);

    array_math_mul(_data, value, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator*(const cArray3& other) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator*(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    cArray3<T>    object(_size);

    array_math_mul(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator*=(const cArray3& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator*=(const cArray3& other)");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    array_math_mul(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator/(T value) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator/(T value) const");

    requireNonEmpty(funcName);

    cArray3<T>    object(_size);

    array_math_div(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator/=(T value)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator/=(T value)");

    requireNonEmpty(funcName);

    array_math_div(_data, value, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator/(const cArray3& other) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator/(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    cArray3<T>    object(_size);

    array_math_div(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cArray3<T>& cArray3<T>::operator/=(const cArray3& other)
{
    const std::string   funcName("cArray3<T>& cArray3<T>::operator/=(const cArray3& other)");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    array_math_div(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cArray3<T> cArray3<T>::operator-(void) const
{
    const std::string   funcName("const cArray3<T> cArray3<T>::operator-(void) const");

    requireNonEmpty(funcName);

    cArray3<T>    object(_size);

    array_math_inv(_data, object._data, getNelement());

    return object;
}

template <typename T>
bool cArray3<T>::operator==(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator==(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    float   threshF = 1e-3f;
    T       threshT = (T) threshF;
    if ((float) threshT != threshF)     threshT = 1;

    return array_reduce_compare(_data, other._data, getNelement(), threshT);
}

template <typename T>
bool cArray3<T>::operator!=(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator!=(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    return !(*this == other);
}

template <typename T>
bool cArray3<T>::operator<(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator<(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    bool    returnVal = true;

    for (size_t i = 0; i < getNelement(); i++) {
        if (_data[i] >= other[i]) {
            returnVal = false;
            break;
        }
    }

    return returnVal;
}

template <typename T>
bool cArray3<T>::operator<=(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator<=(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    bool    returnVal = true;

    for (size_t i = 0; i < getNelement(); i++) {
        if (_data[i] > other[i]) {
            returnVal = false;
            break;
        }
    }

    return returnVal;
}

template <typename T>
bool cArray3<T>::operator>(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator>(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    bool    returnVal = true;

    for (size_t i = 0; i < getNelement(); i++) {
        if (_data[i] <= other[i]) {
            returnVal = false;
            break;
        }
    }

    return returnVal;
}

template <typename T>
bool cArray3<T>::operator>=(const cArray3& other) const
{
    const std::string   funcName("bool cArray3<T>::operator>=(const cArray3& other) const");

    requireNonEmpty(funcName);
    require(getNelement() == other.getNelement(),
            funcName + ": two objects do not have the same #elements");

    bool    returnVal = true;

    for (size_t i = 0; i < getNelement(); i++) {
        if (_data[i] < other[i]) {
            returnVal = false;
            break;
        }
    }

    return returnVal;
}

template <typename T>
T cArray3<T>::getElement(size_t index) const
{
    require(index < getNelement(), "T cArray3<T>::getElement(size_t index) const: index out-of-bound");

    return _data[index];
}

template <typename T>
T cArray3<T>::getElement(const cSize3& index) const
{
    require(index < _size, "T cArray3<T>::getElement(const cSize3& index) const: index out-of-bound");

    return _data[(index[0]*_size[1] + index[1])*_size[2] + index[2]];
}

template <typename T>
void cArray3<T>::setElement(size_t index, T value)
{
    require(index < getNelement(), "void cArray3<T>::setElement(size_t index, T value): index out-of-bound");

    _data[index] = value;
}

template <typename T>
void cArray3<T>::setElement(const cSize3& index, T value)
{
    require(index < _size, "void cArray3<T>::setElement(const cSize3& index, T value): index out-of-bound");

    _data[(index[0]*_size[1] + index[1])*_size[2] + index[2]] = value;
}

template <typename T>
T cArray3<T>::getMax(void) const
{
    requireNonEmpty("T cArray3<T>::getMax(void) const");

    return array_reduce_max(_data, getNelement());
}

template <typename T>
T cArray3<T>::getMin(void) const
{
    requireNonEmpty("T cArray3<T>::getMin(void) const");

    return array_reduce_min(_data, getNelement());
}

template <typename T>
T cArray3<T>::getSum(void) const
{
    requireNonEmpty("T cArray3<T>::getSum(void) const");

    return array_reduce_sum(_data, getNelement());
}

template <typename T>
T cArray3<T>::getSum2(void) const
{
    requireNonEmpty("T cArray3<T>::getSum2(void) const");

    return array_reduce_sum2(_data, getNelement());
}

template <typename T>
T cArray3<T>::getMean(void) const
{
    requireNonEmpty("T cArray3<T>::getMean(void) const");

    return array_reduce_mean(_data, getNelement());
}

template <typename T>
T cArray3<T>::getStd (void) const
{
    requireNonEmpty("T cArray3<T>::getStd (void) const");

    return array_reduce_std(_data, getNelement());
}

template <typename T>
size_t cArray3<T>::getMaxIndx(void) const
{
    requireNonEmpty("size_t cArray3<T>::getMaxIndx(void) const");

    return array_reduce_maxindx(_data, getNelement());
}

template <typename T>
size_t cArray3<T>::getMinIndx(void) const
{
    requireNonEmpty("size_t cArray3<T>::getMinIndx(void) const");

    return array_reduce_minindx(_data, getNelement());
}

template <typename T>
void cArray3<T>::getMax(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getMax(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_max(other._data, other._size[0], other._size[1],
                                       _data,
                                 dim);
            break;
        case 3:
            array_reducendim_max(other._data, other._size[0], other._size[1], other._size[2],
                                       _data,
                                 dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::getMin(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getMin(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_min(other._data, other._size[0], other._size[1],
                                       _data,
                                 dim);
            break;
        case 3:
            array_reducendim_min(other._data, other._size[0], other._size[1], other._size[2],
                                       _data,
                                 dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::getSum(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getSum(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_sum(other._data, other._size[0], other._size[1],
                                       _data,
                                 dim);
            break;
        case 3:
            array_reducendim_sum(other._data, other._size[0], other._size[1], other._size[2],
                                       _data,
                                 dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::getSum2(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getSum2(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_sum2(other._data, other._size[0], other._size[1],
                                        _data,
                                  dim);
            break;
        case 3:
            array_reducendim_sum2(other._data, other._size[0], other._size[1], other._size[2],
                                        _data,
                                  dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::getMean(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getMean(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_mean(other._data, other._size[0], other._size[1],
                                        _data,
                                  dim);
            break;
        case 3:
            array_reducendim_mean(other._data, other._size[0], other._size[1], other._size[2],
                                        _data,
                                  dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::getStd(const cArray3& other, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::getStd(const cArray3& other, eReduceNDim dim)");

    other.requireNonEmpty(funcName);

    memReAllocReduceNDim(other._size, other.getDimension(), dim);

    switch (other.getDimension()) {
        case 2:
            array_reducendim_std(other._data, other._size[0], other._size[1],
                                       _data,
                                 dim);
            break;
        case 3:
            array_reducendim_std(other._data, other._size[0], other._size[1], other._size[2],
                                       _data,
                                 dim);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::memAlloc(const cSize3& size)
{
    requireEmpty("void cArray3<T>::memAlloc(const cSize3& size)");

    _size = size;

    if (getNelement() > 0) {
        array_new(_data, getNelement());
    }
}

template <typename T>
void cArray3<T>::memAllocZero(const cSize3& size)
{
    requireEmpty("void cArray3<T>::memAllocZero(const cSize3& size)");

    _size = size;

    if (getNelement() > 0) {
        array_new_zero(_data, getNelement());
    }
}

template <typename T>
void cArray3<T>::memReAlloc(const cSize3& size)
{
    if (_size != size) {
        if (getNelement() != size[0]*size[1]*size[2]) {
            memFree();
            memAlloc(size);
        }

        _size = size;
    }
    /*if (getNelement() != size[0]*size[1]*size[2]) {
        memFree();
        memAlloc(size);
    }
    else {
        _size = size;
    }*/
}

template <typename T>
void cArray3<T>::memReAllocZero(const cSize3& size)
{
    if (_size != size) {
        if (getNelement() != size[0]*size[1]*size[2]) {
            memFree();
            memAlloc(size);
        }

        _size = size;
    }

    memSetZero();
    /*if (getNelement() != size[0]*size[1]*size[2]) {
        memFree();
        memAllocZero(size);
    }
    else {
        _size = size;
        memSetZero();
    }*/
}

template <typename T>
void cArray3<T>::memReAllocReduceNDim(const cSize3& size, unsigned int ndim, eReduceNDim dim)
{
    const std::string   funcName("void cArray3<T>::memReAllocReduceNDim(const cSize3& size, unsigned int ndim, eReduceNDim dim)");

    switch (ndim) {
        case 2:
            switch (dim) {
                case REDUCE_NDIM_ROW:
                    memReAlloc(cSize3(size[1],1,1));
                    break;
                case REDUCE_NDIM_COL:
                    memReAlloc(cSize3(size[0],1,1));
                    break;
                default:
                    ERROR(funcName, "unsupported dimension mode");
            }
            break;
        case 3:
            switch (dim) {
                case REDUCE_NDIM_ROW:
                    memReAlloc(cSize3(size[1],size[2],1));
                    break;
                case REDUCE_NDIM_COL:
                    memReAlloc(cSize3(size[0],size[2],1));
                    break;
                case REDUCE_NDIM_SEC:
                    memReAlloc(cSize3(size[0],size[1],1));
                    break;
                default:
                    ERROR(funcName, "unsupported dimension mode");
            }
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::memFree(void)
{
    require((getNelement() >  0 && _data != NULL) ||
            (getNelement() == 0 && _data == NULL),
            "void cArray3<T>::memFree(void): data and size mismatch");

    if (getNelement() > 0) {
        array_delete(_data);
    }

    _size = 0;
}

template <typename T>
void cArray3<T>::memSetZero(void)
{
    requireNonEmpty("void cArray3<T>::memSetZero(void)");

    array_memset(_data, 0, getNelement());
}

template <typename T>
void cArray3<T>::memSetVal(T value)
{
    requireNonEmpty("void cArray3<T>::memSetVal(T value)");

    array_setval(_data, value, getNelement());
}

template <typename T>
void cArray3<T>::memSwap(cArray3& other)
{
    require(this != &other, "void cArray3<T>::memSwap(cArray3& other): two objects are the same");

    std::swap(_size[0], other._size[0]);
    std::swap(_size[1], other._size[1]);
    std::swap(_size[2], other._size[2]);
    std::swap(_data, other._data);
}

template <typename T>
void cArray3<T>::memReShape(const cSize3& size)
{
    require(getNelement() == size.getProduct(), "void cArray3<T>::memReShape(const cSize3& size): invalid new size");

    _size = size;
}

template <typename T>
void cArray3<T>::printData(std::string message, int width, int precision, bool scientific) const
{
    const std::string   funcName("void cArray3<T>::printData(std::string message, int width, int precision, bool scientific) const");

    requireNonEmpty(funcName);

    if (message.size())
        std::cout << message << ":" << std::endl;

    if (_size[2] > 1) {
        array_print(_data, _size[0], _size[1], _size[2], width, precision, scientific);
    }
    else if(_size[1] > 1) {
        array_print(_data, _size[0], _size[1], width, precision, scientific);
    }
    else {
        array_print(_data, _size[0], width, precision, scientific);
    }

    /* This code is incorrect for size [1 3 4]
    switch (getDimension()) {
        case 1:     array_print(_data, _size[0]                    , width, precision); break;
        case 2:     array_print(_data, _size[0], _size[1]          , width, precision); break;
        case 3:     array_print(_data, _size[0], _size[1], _size[2], width, precision); break;
        default:    ERROR(funcName, "unsupported dimension");
    }*/
}

template <typename T>
void cArray3<T>::printSize(std::string message, int width, int precision) const
{
    if (message.size())
        std::cout << message << ":" << std::endl;

    std::cout.precision(precision);
    std::setw(width);

    std::cout << std::setw(width) << "size = " << _size
              << std::endl << std::endl;
}

template <typename T>
void cArray3<T>::autocontrast8bit(T min, T max)
{
    size_t n = getNelement();

    requireNonEmpty("void cArray3<T>::autocontrast8bit(T min, T max)");

    *this -= min;
    *this /= (max-min) / 255;

    /* clip values that are no more in range */
    #pragma omp parallel for
    for (size_t i = 0; i < n; i++) {
        if (_data[i] <= 0)
            _data[i] = 0;
        if (_data[i] > 255)
            _data[i] = 255;
    }
}

template <typename T>
void cArray3<T>::autocontrast8bit(int cutoff)
{
    size_t histo[256] = {0};
    size_t n       = getNelement();
    size_t cut     = n * cutoff / 100;
    T min = 0, max = 255;
    ssize_t remains;

    requireNonEmpty("void cArray3<T>::autocontrast8bit(int cutoff)");

    normalize8bit();

    /* build histogram */
    for (size_t i=0; i<n; i++)
    {
      histo[(int)_data[i]]++;
    }

    /* find min */
    remains = cut;
    for (size_t i=0; i<256; i++)
    {
      remains -= histo[i];
      if (remains <= 0)
      {
        min = (T)i;
        break;
      }
    }

    /* find max */
    remains = cut;
    for (size_t i=255; i>=1; i--)
    {
      remains -= histo[i];
      if (remains <= 0)
      {
        max = (T)i;
        break;
      }
    }

    /* stretch to that min/max */
    autocontrast8bit(min, max);
}

template <typename T>
void cArray3<T>::normalize8bit(void)
{
    requireNonEmpty("void cArray3<T>::normalize8bit(void)");

    autocontrast8bit(getMin(), getMax());
}

template <typename T>
void cArray3<T>::normalize8bit(const cArray3& other)
{
    other.requireNonEmpty("void cArray3<T>::normalize8bit(const cArray3& other)");

    *this = other;

    normalize8bit();
}

template <typename T>
void cArray3<T>::opMathAbs(void)
{
    requireNonEmpty("void cArray3<T>::opMathAbs(void)");

    array_math_abs(_data, getNelement());
}

template <typename T>
void cArray3<T>::opMathInv(void)
{
    requireNonEmpty("void cArray3<T>::opMathInv(void)");

    array_math_inv(_data, getNelement());
}

template <typename T>
void cArray3<T>::opMathNorm(void)
{
    requireNonEmpty("void cArray3<T>::opMathNorm(void)");

    array_math_norm(_data, getNelement());
}

template <typename T>
void cArray3<T>::opMathSqr(void)
{
    requireNonEmpty("void cArray3<T>::opMathSqr(void)");

    array_math_sqr(_data, getNelement());
}

template <typename T>
void cArray3<T>::opMathSqrt(void)
{
    requireNonEmpty("void cArray3<T>::opMathSqrt(void)");

    array_math_sqrt(_data, getNelement());
}

template <typename T>
void cArray3<T>::getRayRow(const cArray3& other, size_t iCol, size_t iSec)
{
    const std::string   funcName("void cArray3<T>::getRayRow(const cArray3& other, size_t iCol, size_t iSec)");

    other.requireNonEmpty(funcName);
    require(iCol < other._size[1], funcName + ": iCol out of range");
    require(iSec < other._size[2], funcName + ": iSec out of range");

    memReAlloc(cSize3(other._size[0],1,1));

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < other._size[0]; iRow++) {
        _data[iRow] = other._data[sub2ind(iRow,iCol,iSec,other._size[1],other._size[2])];
    }
}

template <typename T>
void cArray3<T>::getRayCol(const cArray3& other, size_t iRow, size_t iSec)
{
    const std::string   funcName("void cArray3<T>::getRayCol(const cArray3& other, size_t iRow, size_t iSec)");

    other.requireNonEmpty(funcName);
    require(iRow < other._size[0], funcName + ": iRow out of range");
    require(iSec < other._size[2], funcName + ": iSec out of range");

    memReAlloc(cSize3(other._size[1],1,1));

    #pragma omp parallel for
    for (size_t iCol = 0; iCol < other._size[1]; iCol++) {
        _data[iCol] = other._data[sub2ind(iRow,iCol,iSec,other._size[1],other._size[2])];
    }
}

template <typename T>
void cArray3<T>::getRaySec(const cArray3& other, size_t iRow, size_t iCol)
{
    const std::string   funcName("void cArray3<T>::getRaySec(const cArray3& other, size_t iRow, size_t iCol)");

    other.requireNonEmpty(funcName);
    require(iRow < other._size[0], funcName + ": iRow out of range");
    require(iCol < other._size[1], funcName + ": iCol out of range");

    memReAlloc(cSize3(other._size[2],1,1));

    array_memcpy(_data, other._data + sub2ind(iRow,iCol,0,other._size[1],other._size[2]), other._size[2]);
}

template <typename T>
void cArray3<T>::getSlideRow(const cArray3& other, size_t iRow)
{
    const std::string   funcName("void cArray3<T>::getSlideRow(const cArray3& other, size_t iRow)");

    other.requireNonEmpty(funcName);
    require(iRow < other._size[0], funcName + ": iRow out of range");

    memReAlloc(cSize3(other._size[1],other._size[2],1));

    array_memcpy(_data, other._data + iRow*getNelement(), getNelement());
}

template <typename T>
void cArray3<T>::getSlideCol(const cArray3& other, size_t iCol)
{
    const std::string   funcName("void cArray3<T>::getSlideCol(const cArray3& other, size_t iCol)");

    other.requireNonEmpty(funcName);
    require(iCol < other._size[1], funcName + ": iCol out of range");

    memReAlloc(cSize3(other._size[0],other._size[2],1));

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < _size[0]; iRow++) {
        for (size_t iSec = 0; iSec < _size[1]; iSec++) {
            _data[sub2ind(iRow,iSec,_size[1])] = other._data[sub2ind(iRow,iCol,iSec,other._size[1],other._size[2])];
        }
    }
}

template <typename T>
void cArray3<T>::getSlideSec(const cArray3& other, size_t iSec)
{
    const std::string   funcName("void cArray3<T>::getSlideSec(const cArray3& other, size_t iSec)");

    other.requireNonEmpty(funcName);
    require(iSec < other._size[2], funcName + ": iSec out of range");

    memReAlloc(cSize3(other._size[0],other._size[1],1));

    #pragma omp parallel for
    for (size_t iRow = 0; iRow < _size[0]; iRow++) {
        for (size_t iCol = 0; iCol < _size[1]; iCol++) {
            _data[sub2ind(iRow,iCol,_size[1])] = other._data[sub2ind(iRow,iCol,iSec,_size[1],_size[2])];
        }
    }
}

template <typename T>
void cArray3<T>::opCircShift(const cOffset3& offset)
{
    const std::string   funcName("void cArray3<T>::opCircShift(const cOffset3& offset)");

    requireNonEmpty(funcName);

    cArray3<T>    other(_size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_circshift(_data, other._data,
                            _size[0],
                            offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_circshift(_data, other._data,
                             _size[0],  _size[1],
                            offset[0], offset[1]);
            break;
        case 3:
            array_circshift(_data, other._data,
                             _size[0],  _size[1],  _size[2],
                            offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opCircShift(const cArray3& other, const cOffset3& offset)
{
    const std::string   funcName("void cArray3<T>::opCircShift(const cArray3& other, const cOffset3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_circshift(other._data, _data,
                            _size[0],
                            offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_circshift(other._data, _data,
                             _size[0],  _size[1],
                            offset[0], offset[1]);
            break;
        case 3:
            array_circshift(other._data, _data,
                             _size[0],  _size[1],  _size[2],
                            offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opShift(const cOffset3& offset)
{
    const std::string   funcName("void cArray3<T>::opShift(const cOffset3& offset)");

    requireNonEmpty(funcName);

    cArray3<T>    other(_size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_shift(_data, other._data,
                         _size[0],
                        offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_shift(_data, other._data,
                         _size[0],  _size[1],
                        offset[0], offset[1]);
            break;
        case 3:
            array_shift(_data, other._data,
                         _size[0],  _size[1],  _size[2],
                        offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opShift(const cArray3& other, const cOffset3& offset)
{
    const std::string   funcName("void cArray3<T>::opShift(const cArray3& other, const cOffset3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_shift(other._data, _data,
                         _size[0],
                        offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_shift(other._data, _data,
                         _size[0],  _size[1],
                        offset[0], offset[1]);
            break;
        case 3:
            array_shift(other._data, _data,
                         _size[0],  _size[1],  _size[2],
                        offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opPad(const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opPad(const cSize3& size, const cSize3& offset)");

    requireNonEmpty(funcName);
    require(size >= _size,          funcName + ": requested size is not larger than current size");
    require(size >= _size + offset, funcName + ": offset value is too large");

    cArray3<T>    other(size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_pad(      _data,       _size[0],
                      other._data, other._size[0],
                                        offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_pad(      _data,       _size[0],       _size[1],
                      other._data, other._size[0], other._size[1],
                                        offset[0],      offset[1]);
            break;
        case 3:
            array_pad(      _data,       _size[0],       _size[1],       _size[2],
                      other._data, other._size[0], other._size[1], other._size[2],
                                        offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opPad(const cArray3& other, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opPad(const cArray3& other, const cSize3& size, const cSize3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other,               funcName + ": two objects are the same");
    require(size >= other._size,          funcName + ": other object is too large");
    require(size >= other._size + offset, funcName + ": offset value is too large");

    memReAllocZero(size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_pad(other._data, other._size[0],
                            _data,       _size[0],
                                        offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_pad(other._data, other._size[0], other._size[1],
                            _data,       _size[0],       _size[1],
                                        offset[0],      offset[1]);
            break;
        case 3:
            array_pad(other._data, other._size[0], other._size[1], other._size[2],
                            _data,       _size[0],       _size[1],       _size[2],
                                        offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opCrop(const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opCrop(const cSize3& size, const cSize3& offset)");

    requireNonEmpty(funcName);
    require(_size >= size,          funcName + ": requested size is not smaller than current size");
    require(_size >= size + offset, funcName + ": offset value is too large");

    cArray3<T>    other(size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_crop(      _data,       _size[0],
                       other._data, other._size[0],
                                         offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_crop(      _data,       _size[0],       _size[1],
                       other._data, other._size[0], other._size[1],
                                         offset[0],      offset[1]);
            break;
        case 3:
            array_crop(      _data,       _size[0],       _size[1],       _size[2],
                       other._data, other._size[0], other._size[1], other._size[2],
                                         offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opCrop(const cArray3& other, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opCrop(const cArray3& other, const cSize3& size, const cSize3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other,               funcName + ": two objects are the same");
    require(other._size >= size,          funcName + ": other object is too small");
    require(other._size >= size + offset, funcName + ": offset value is too large");

    memReAlloc(size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_crop(other._data, other._size[0],
                             _data,       _size[0],
                                         offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_crop(other._data, other._size[0], other._size[1],
                             _data,       _size[0],       _size[1],
                                         offset[0],      offset[1]);
            break;
        case 3:
            array_crop(other._data, other._size[0], other._size[1], other._size[2],
                             _data,       _size[0],       _size[1],       _size[2],
                                         offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opReplace(T value, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opReplace(T value, const cSize3& size, const cSize3& offset)");

    requireNonEmpty("");
    require(size <= _size,        funcName + ": replacing region is too large");
    require(size+offset <= _size, funcName + ": offset value is too large");

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_replace(_data,  _size[0],
                          value,   size[0],
                                 offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_replace(_data,  _size[0],  _size[1],
                          value,   size[0],   size[1],
                                 offset[0], offset[1]);
            break;
        case 3:
            array_replace(_data,  _size[0],  _size[1],  _size[2],
                          value,   size[0],   size[1],   size[2],
                                 offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opReplace(const cArray3& other, const cSize3& offset)
{
    const std::string   funcName("void cArray3<T>::opReplace(const cArray3& other, const cSize3& offset)");

    requireNonEmpty(funcName);
    other.requireNonEmpty(funcName);
    require(this != &other,              funcName + ": two objects are the same");
    require(other._size <= _size,        funcName + ": other object is too large");
    require(other._size+offset <= _size, funcName + ": offset value is too large");

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_replace<T,T>(      _data,       _size[0],
                               other._data, other._size[0],
                                                 offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            array_replace<T,T>(      _data,       _size[0],       _size[1],
                               other._data, other._size[0], other._size[1],
                                                 offset[0],      offset[1]);
            break;
        case 3:
            array_replace<T,T>(      _data,       _size[0],       _size[1],       _size[2],
                               other._data, other._size[0], other._size[1], other._size[2],
                                                 offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cArray3<T>::opPermute(ePermute permute)
{
    const std::string   funcName("void cArray3<T>::opPermute(ePermute permute)");

    requireNonEmpty(funcName);

    cArray3<T>    other;

    switch (permute) {
        case PERMUTE2D:
            require(_size[2] == 1, funcName + ": the 3rd dimension is not 1");
            other.memReAlloc(cSize3(_size[1],_size[0],1));
            break;
        case PERMUTE3D_123:
            ERROR(funcName, "no need to permute");
            break;
        case PERMUTE3D_132:
            other.memReAlloc(cSize3(_size[0],_size[2],_size[1]));
            break;
        case PERMUTE3D_213:
            other.memReAlloc(cSize3(_size[1],_size[0],_size[2]));
            break;
        case PERMUTE3D_231:
            other.memReAlloc(cSize3(_size[1],_size[2],_size[0]));
            break;
        case PERMUTE3D_312:
            other.memReAlloc(cSize3(_size[2],_size[0],_size[1]));
            break;
        case PERMUTE3D_321:
            other.memReAlloc(cSize3(_size[2],_size[1],_size[0]));
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }

    switch (permute) {
        case PERMUTE2D:
            array_permute(_data, other._data,
                          _size[0], _size[1]);
            break;
        case PERMUTE3D_123:
        case PERMUTE3D_132:
        case PERMUTE3D_213:
        case PERMUTE3D_231:
        case PERMUTE3D_312:
        case PERMUTE3D_321:
            array_permute(_data, other._data,
                          _size[0], _size[1], _size[2],
                          permute);
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opPermute(const cArray3& other, ePermute permute)
{
    const std::string   funcName("void cArray3<T>::opPermute(const cArray3& other, ePermute permute)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    switch (permute) {
        case PERMUTE2D:
            memReAlloc(cSize3(other._size[1],other._size[0],1));
            break;
        case PERMUTE3D_123:
            ERROR(funcName, "no need to permute, memcpy is faster");
            break;
        case PERMUTE3D_132:
            memReAlloc(cSize3(other._size[0],other._size[2],other._size[1]));
            break;
        case PERMUTE3D_213:
            memReAlloc(cSize3(other._size[1],other._size[0],other._size[2]));
            break;
        case PERMUTE3D_231:
            memReAlloc(cSize3(other._size[1],other._size[2],other._size[0]));
            break;
        case PERMUTE3D_312:
            memReAlloc(cSize3(other._size[2],other._size[0],other._size[1]));
            break;
        case PERMUTE3D_321:
            memReAlloc(cSize3(other._size[2],other._size[1],other._size[0]));
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }

    switch (permute) {
        case PERMUTE2D:
            array_permute(other._data, _data,
                          other._size[0], other._size[1]);
            break;
        case PERMUTE3D_123:
        case PERMUTE3D_132:
        case PERMUTE3D_213:
        case PERMUTE3D_231:
        case PERMUTE3D_312:
        case PERMUTE3D_321:
            array_permute(other._data, _data,
                          other._size[0], other._size[1], other._size[2],
                          permute);
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }
}

template <typename T>
void cArray3<T>::opCopy(const cArray3& other, const cSize3& size, const cOffset3& offset)
{
    const std::string   funcName("void cArray3<T>::opCopy(const cArray3& other, const cSize3& size, const cOffset3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAllocZero(size);

    cSize3      sizeTmp, offsetTmp;

    if (offset[0] >= 0) { require(offset[0] <= (ptrdiff_t) other._size[0], "out of boundary 10 in " + funcName); }
    else                { require(offset[0]  + (ptrdiff_t)  _size[0] >= 0, "out of boundary 11 in " + funcName); }
    if (offset[1] >= 0) { require(offset[1] <= (ptrdiff_t) other._size[1], "out of boundary 20 in " + funcName); }
    else                { require(offset[1]  + (ptrdiff_t)  _size[1] >= 0, "out of boundary 21 in " + funcName); }
    if (offset[2] >= 0) { require(offset[2] <= (ptrdiff_t) other._size[2], "out of boundary 30 in " + funcName); }
    else                { require(offset[2]  + (ptrdiff_t)  _size[2] >= 0, "out of boundary 31 in " + funcName); }

    offsetTmp[0] = (offset[0] >= 0) ? offset[0] : 0;
    offsetTmp[1] = (offset[1] >= 0) ? offset[1] : 0;
    offsetTmp[2] = (offset[2] >= 0) ? offset[2] : 0;

    sizeTmp[0] = (offsetTmp[0]+_size[0] <= other._size[0]) ? (_size[0] + offset[0] - offsetTmp[0]) : (other._size[0] - offsetTmp[0]);
    sizeTmp[1] = (offsetTmp[1]+_size[1] <= other._size[1]) ? (_size[1] + offset[1] - offsetTmp[1]) : (other._size[1] - offsetTmp[1]);
    sizeTmp[2] = (offsetTmp[2]+_size[2] <= other._size[2]) ? (_size[2] + offset[2] - offsetTmp[2]) : (other._size[2] - offsetTmp[2]);

    cArray3<T>              object(sizeTmp);

    object.opCrop(other, sizeTmp, offsetTmp);

    opReplace(object, cSize3(offsetTmp[0] - offset[0],
                             offsetTmp[1] - offset[1],
                             offsetTmp[2] - offset[2]));
}

template <typename T>
void cArray3<T>::opIndexCol2Row(void)
{
    const std::string   funcName("void cArray3<T>::opIndexCol2Row(void)");

    requireNonEmpty(funcName);

    cArray3<T>    other(_size);

    switch (getDimension()) {
        case 2:
            array_index_col2row(_data, other._data,
                                _size[0], _size[1]);
            break;
        case 3:
            array_index_col2row(_data, other._data,
                                _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::opIndexRow2Col(void)
{
    const std::string   funcName("void cArray3<T>::opIndexRow2Col(void)");

    requireNonEmpty(funcName);

    cArray3<T>    other(_size);

    switch (getDimension()) {
        case 2:
            array_index_row2col(_data, other._data,
                                _size[0], _size[1]);
            break;
        case 3:
            array_index_row2col(_data, other._data,
                                _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cArray3<T>::arrayImport(const T* const array, size_t nRow, size_t nCol, size_t nSec)
{
    const std::string   funcName("void cArray3<T>::arrayImport(const T* const array, size_t nRow, size_t nCol, size_t nSec)");

    if (nRow*nCol*nSec != 0) {
        require(array != NULL, funcName + ": input array cannot be NULL");
    }
    else {
        require(array == NULL, funcName + ": input array must be NULL");
    }

    memReAlloc(cSize3(nRow, nCol, nSec));
    array_memcpy(_data, array, getNelement());
}

template <typename T>
void cArray3<T>::arrayHandle(T* const array, size_t nRow, size_t nCol, size_t nSec)
{
    const std::string   funcName("void cArray3<T>::arrayHandle(T* const array, size_t nRow, size_t nCol, size_t nSec)");

    if (nRow*nCol*nSec != 0) {
        require(array != NULL, funcName + ": input array cannot be NULL");
    }
    else {
        require(array == NULL, funcName + ": input array must be NULL");
    }

    memFree();
    _size = cSize3(nRow, nCol, nSec);
    _data = array;
}

template <typename T>
void cArray3<T>::arrayExport(T* array)
{
    const std::string   funcName("void cArray3<T>::arrayExport(T* array)");

    require(array == NULL, funcName + ": input array must be NULL");

    array_new(array, getNelement());

    array_memcpy(array, _data, getNelement());
}

template <typename T>
cArray3<T> operator-(T value, const cArray3<T>& array)
{
    const std::string   funcName("cVector3<T> operator-(T d, const cArray3<T>& array)");

    array.requireNonEmpty(funcName);

    cArray3<T>    object(array.getSize());

    #pragma omp parallel for
    for (size_t i = 0; i < array.getNelement(); i++) {
        object[i] = value - array[i];
    }

    return object;
}

// instantiation
template class cArray3<int32_t >;
template class cArray3<uint32_t>;
template class cArray3<int64_t >;
template class cArray3<uint64_t>;
template class cArray3<float   >;
template class cArray3<double  >;

template cArray3<int32_t > operator-(int32_t  value, const cArray3<int32_t >& array);
template cArray3<uint32_t> operator-(uint32_t value, const cArray3<uint32_t>& array);
template cArray3<int64_t > operator-(int64_t  value, const cArray3<int64_t >& array);
template cArray3<uint64_t> operator-(uint64_t value, const cArray3<uint64_t>& array);
template cArray3<float   > operator-(float    value, const cArray3<float   >& array);
template cArray3<double  > operator-(double   value, const cArray3<double  >& array);

} // namespace gem
