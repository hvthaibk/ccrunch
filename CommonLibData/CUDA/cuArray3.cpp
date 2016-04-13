/***********************************************************************
 *  File:       cuArray3.cpp
 *
 *  Purpose:    Implementation of a 3D array class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuArray3.hpp"

#include "config.hpp"
#include "macro.hpp"
#include "array.cuh"

namespace gem {

template <typename T>
cuArray3<T>::cuArray3(const cSize3& size)
{
    _size = size;
    _data = NULL;

    if (getNelement() > 0) {
        cuda_arrayDev_new_zero(_data, getNelement());
    }
    else {
#ifndef NDEBUG
        //WARNING("cuArray3::cuArray3", "an empty array is requested");
#endif
    }
}

template <typename T>
cuArray3<T>::cuArray3(const cuArray3& other)
{
    _size = other._size;
    _data = NULL;

    if (getNelement() > 0) {
        cuda_arrayDev_new(_data, getNelement());
        cuda_array_memcpy_d2d(_data, other._data, getNelement());
    }
    else {
#ifndef NDEBUG
        WARNING("cuArray3::cuArray3", "an empty array is requested");
#endif
    }
}

template <typename T>
cuArray3<T>::cuArray3(const cArray3<T>& other)
{
    _size = other.getSize();
    _data = NULL;

    if (getNelement() > 0) {
        cuda_arrayDev_new(_data, getNelement());
        cuda_array_memcpy_h2d(_data, other.getAddrData(), getNelement());
    }
    else {
#ifndef NDEBUG
        WARNING("cuArray3::cuArray3", "an empty array is requested");
#endif
    }
}

template <typename T>
cuArray3<T>::~cuArray3()
{
    memFree();
}

template <typename T>
unsigned int cuArray3<T>::getDimension(void) const
{
    unsigned int    ndim = (_size[0]>1)+(_size[1]>1)+(_size[2]>1);

    if (_size[1] > 1 && ndim < 2)  ndim = 2;
    if (_size[2] > 1 && ndim < 3)  ndim = 3;

    return ndim;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator=(T value)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator=(T value)");

    cuda_array_setval(_data, value, getNelement());

    return *this;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator=(const cuArray3& other)
{
    require(this != &other, "cuArray3<T>& cuArray3<T>::operator=(const cuArray3& other): two objects are the same");

    memReAlloc(other._size);

    if (getNelement() > 0) {
        cuda_array_memcpy_d2d(_data, other._data, getNelement());
    }

    return *this;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator=(const cArray3<T>& other)
{
    memReAlloc(other.getSize());

    if (getNelement() > 0) {
        cuda_array_memcpy_h2d(_data, other.getAddrData(), getNelement());
    }

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator+(T value) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator+(T value) const");

    cuArray3<T>    object(_size);

    cuda_array_math_add(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator+=(T value)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator+=(T value)");

    cuda_array_math_add(_data, value, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator+(const cuArray3& other) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator+(const cuArray3& other) const");
    require(getNelement() == other.getNelement(),
            "const cuArray3<T> cuArray3<T>::operator+(const cuArray3& other) const: two objects do not have the same number of elements");

    cuArray3<T>    object(_size);

    cuda_array_math_add(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator+=(const cuArray3& other)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator+=(const cuArray3& other)");
    require(getNelement() == other.getNelement(),
            "cuArray3<T>& cuArray3<T>::operator+=(const cuArray3& other): two objects do not have the same number of elements");

    cuda_array_math_add(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator-(T value) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator-(T value) const");

    cuArray3<T>    object(_size);

    cuda_array_math_sub(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator-=(T value)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator-=(T value)");

    cuda_array_math_sub(_data, value, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator-(const cuArray3& other) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator-(const cuArray3& other) const");
    require(getNelement() == other.getNelement(),
            "const cuArray3<T> cuArray3<T>::operator-(const cuArray3& other) const: two objects do not have the same number of elements");

    cuArray3<T>    object(_size);

    cuda_array_math_sub(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator-=(const cuArray3& other)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator-=(const cuArray3& other)");
    require(getNelement() == other.getNelement(),
            "cuArray3<T>& cuArray3<T>::operator-=(const cuArray3& other): two objects do not have the same number of elements");

    cuda_array_math_sub(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator*(T value) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator*(T value) const");

    cuArray3<T>    object(_size);

    cuda_array_math_mul(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator*=(T value)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator*=(T value)");

    cuda_array_math_mul(_data, value, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator*(const cuArray3& other) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator*(const cuArray3& other) const");
    require(getNelement() == other.getNelement(),
            "const cuArray3<T> cuArray3<T>::operator*(const cuArray3& other) const: two objects do not have the same number of elements");

    cuArray3<T>    object(_size);

    cuda_array_math_mul(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator*=(const cuArray3& other)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator*=(const cuArray3& other)");
    require(getNelement() == other.getNelement(),
            "cuArray3<T>& cuArray3<T>::operator*=(const cuArray3& other): two objects do not have the same number of elements");

    cuda_array_math_mul(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator/(T value) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator/(T value) const");

    cuArray3<T>    object(_size);

    cuda_array_math_div(object._data, _data, value, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator/=(T value)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator/=(T value)");

    cuda_array_math_div(_data, value, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator/(const cuArray3& other) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator/(const cuArray3& other) const");
    require(getNelement() == other.getNelement(),
            "const cuArray3<T> cuArray3<T>::operator/(const cuArray3& other) const: two objects do not have the same number of elements");

    cuArray3<T>    object(_size);

    cuda_array_math_div(object._data, _data, other._data, getNelement());

    return object;
}

template <typename T>
cuArray3<T>& cuArray3<T>::operator/=(const cuArray3& other)
{
    requireNonEmpty("cuArray3<T>& cuArray3<T>::operator/=(const cuArray3& other)");
    require(getNelement() == other.getNelement(),
            "cuArray3<T>& cuArray3<T>::operator/=(const cuArray3& other): two objects do not have the same number of elements");

    cuda_array_math_div(_data, other._data, getNelement());

    return *this;
}

template <typename T>
const cuArray3<T> cuArray3<T>::operator-(void) const
{
    requireNonEmpty("const cuArray3<T> cuArray3<T>::operator-(void) const");

    cuArray3<T>    object(_size);

    cuda_array_math_inv(_data, object._data, getNelement());

    return object;
}

template <typename T>
T cuArray3<T>::getMax(void) const
{
    requireNonEmpty("T cuArray3<T>::getMax(void) const");

#ifdef __GEM_USE_THRUST__
    return cuda_array_reduce_thrust(_data, getNelement(), REDUCE_MAX);
#else
    return cuda_array_reduce_max(_data, getNelement());
#endif
}

template <typename T>
T cuArray3<T>::getMin(void) const
{
    requireNonEmpty("T cuArray3<T>::getMin(void) const");

#ifdef __GEM_USE_THRUST__
    return cuda_array_reduce_thrust(_data, getNelement(), REDUCE_MIN);
#else
    return cuda_array_reduce_min(_data, getNelement());
#endif
}

template <typename T>
T cuArray3<T>::getSum(void) const
{
    requireNonEmpty("T cuArray3<T>::getSum(void) const");

#ifdef __GEM_USE_THRUST__
    return cuda_array_reduce_thrust(_data, getNelement(), REDUCE_SUM);
#else
    return cuda_array_reduce_sum(_data, getNelement());
#endif
}

template <typename T>
T cuArray3<T>::getSum2(void) const
{
    requireNonEmpty("T cuArray3<T>::getSum2(void) const");

#ifdef __GEM_USE_THRUST__
    return cuda_array_reduce_thrust(_data, getNelement(), REDUCE_SUM2);
#else
    return cuda_array_reduce_sum2(_data, getNelement());
#endif
}

template <typename T>
T cuArray3<T>::getMean(void) const
{
    requireNonEmpty("T cuArray3<T>::getMean(void) const");

    return cuda_array_reduce_mean(_data, getNelement());
}

template <typename T>
T cuArray3<T>::getStd (void) const
{
    requireNonEmpty("T cuArray3<T>::getStd (void) const");

    return cuda_array_reduce_std(_data, getNelement());
}

template <typename T>
size_t cuArray3<T>::getMaxIndx(void) const
{
    requireNonEmpty("size_t cuArray3<T>::getMaxIndx(void) const");

#ifdef __GEM_USE_THRUST__
    return (size_t) cuda_array_reduce_thrust(_data, getNelement(), REDUCE_MAX_INDEX);
#else
    return cuda_array_reduce_maxindx(_data, getNelement());
#endif
}

template <typename T>
size_t cuArray3<T>::getMinIndx(void) const
{
    requireNonEmpty("size_t cuArray3<T>::getMinIndx(void) const");

#ifdef __GEM_USE_THRUST__
    return (size_t) cuda_array_reduce_thrust(_data, getNelement(), REDUCE_MIN_INDEX);
#else
    return cuda_array_reduce_minindx(_data, getNelement());
#endif
}

template <typename T>
void cuArray3<T>::memAlloc(const cSize3& size)
{
    requireEmpty("void cuArray3<T>::memAlloc(const cSize3& size)");

    _size = size;

    if (getNelement() > 0) {
        cuda_arrayDev_new(_data, getNelement());
    }
}

template <typename T>
void cuArray3<T>::memAllocZero(const cSize3& size)
{
    requireEmpty("void cuArray3<T>::memAllocZero(const cSize3& size)");

    _size = size;

    if (getNelement() > 0) {
        cuda_arrayDev_new_zero(_data, getNelement());
    }
}

template <typename T>
void cuArray3<T>::memReAlloc(const cSize3& size)
{
    if (getNelement() != size[0]*size[1]*size[2]) {
        memFree();
        memAlloc(size);
    }
    else {
        _size = size;
    }
}

template <typename T>
void cuArray3<T>::memReAllocZero(const cSize3& size)
{
    if (getNelement() != size[0]*size[1]*size[2]) {
        memFree();
        memAllocZero(size);
    }
    else {
        _size = size;
        memSetZero();
    }
}

template <typename T>
void cuArray3<T>::memFree(void)
{
    require((getNelement() >  0 && _data != NULL) ||
            (getNelement() == 0 && _data == NULL),
            "void cArray3<T>::memFree(void): data and size mismatch");

    if (getNelement() > 0) {
        cuda_arrayDev_delete(_data);
    }

    _size = 0;
}

template <typename T>
void cuArray3<T>::memSetZero(void)
{
    requireNonEmpty("void cuArray3<T>::memSetZero(void)");

    cuda_array_memset(_data, 0, getNelement());
}

template <typename T>
void cuArray3<T>::memSetVal(T value)
{
    requireNonEmpty("void cuArray3<T>::memSetVal(T value)");

    cuda_array_setval(_data, value, getNelement());
}

template <typename T>
void cuArray3<T>::memSwap(cuArray3& other)
{
    require(this != &other, "void cuArray3<T>::memSwap(cuArray3& other): two objects are the same");

    std::swap(_size[0], other._size[0]);
    std::swap(_size[1], other._size[1]);
    std::swap(_size[2], other._size[2]);
    std::swap(_data, other._data);
}

template <typename T>
void cuArray3<T>::printSize(std::string message, int width, int precision) const
{
    if (message.size())
        std::cout << message << ":" << std::endl;

    std::cout.precision(precision);
    std::setw(width);

    std::cout << std::setw(width) << "size = " << _size
              << std::endl << std::endl;
}

template <typename T>
void cuArray3<T>::normalize8bit(void)
{
    requireNonEmpty("void cuArray3<T>::normalize8bit(void)");

    T   dataMax = getMax();
    T   dataMin = getMin();

    *this -= dataMin;
    *this /= (dataMax-dataMin) / 255;
}

template <typename T>
void cuArray3<T>::opMathAbs(void)
{
    requireNonEmpty("void cuArray3<T>::opMathAbs(void)");

    cuda_array_math_abs(_data, getNelement());
}

template <typename T>
void cuArray3<T>::opMathInv(void)
{
    requireNonEmpty("void cuArray3<T>::opMathInv(void)");

    cuda_array_math_inv(_data, getNelement());
}

template <typename T>
void cuArray3<T>::opMathNorm(void)
{
    requireNonEmpty("void cuArray3<T>::opMathNorm(void)");

    cuda_array_math_norm(_data, getNelement());
}

template <typename T>
void cuArray3<T>::opMathSqr(void)
{
    requireNonEmpty("void cuArray3<T>::opMathSqr(void)");

    cuda_array_math_sqr(_data, getNelement());
}

template <typename T>
void cuArray3<T>::opMathSqrt(void)
{
    requireNonEmpty("void cuArray3<T>::opMathSqrt(void)");

    cuda_array_math_sqrt(_data, getNelement());
}

template <typename T>
void cuArray3<T>::opCircShift(const cOffset3& offset)
{
    const std::string   funcName("void cuArray3<T>::opCircShift(const cOffset3& offset)");

    requireNonEmpty(funcName);

    cuArray3<T>    other(_size);

    switch (getDimension()) {
        case 1:
            cuda_array_circshift(_data, other._data,
                                  _size[0],
                                 offset[0]);
            break;
        case 2:
            cuda_array_circshift(_data, other._data,
                                  _size[0],  _size[1],
                                 offset[0], offset[1]);
            break;
        case 3:
            cuda_array_circshift(_data, other._data,
                                  _size[0],  _size[1],  _size[2],
                                 offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opCircShift(const cuArray3& other, const cOffset3& offset)
{
    const std::string   funcName("void cuArray3<T>::opCircShift(const cuArray3& other, const cOffset3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_array_circshift(other._data, _data,
                                  _size[0],
                                 offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_array_circshift(other._data, _data,
                                  _size[0],  _size[1],
                                 offset[0], offset[1]);
            break;
        case 3:
            cuda_array_circshift(other._data, _data,
                                  _size[0],  _size[1],  _size[2],
                                 offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opShift(const cOffset3& offset)
{
    const std::string   funcName("void cuArray3<T>::opShift(const cOffset3& offset)");

    requireNonEmpty(funcName);

    cuArray3<T>    other(_size);

    switch (getDimension()) {
        case 1:
            cuda_array_shift(_data, other._data,
                              _size[0],
                             offset[0]);
            break;
        case 2:
            cuda_array_shift(_data, other._data,
                              _size[0],  _size[1],
                             offset[0], offset[1]);
            break;
        case 3:
            cuda_array_shift(_data, other._data,
                              _size[0],  _size[1],  _size[2],
                             offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opShift(const cuArray3& other, const cOffset3& offset)
{
    const std::string   funcName("void cuArray3<T>::opShift(const cuArray3& other, const cOffset3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other, funcName + ": two objects are the same");

    memReAlloc(other._size);

    switch (getDimension()) {
        case 1:
            require(offset[1] == 0, funcName + ": unused offset != 0");
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_array_shift(other._data, _data,
                              _size[0],
                             offset[0]);
            break;
        case 2:
            require(offset[2] == 0, funcName + ": unused offset != 0");

            cuda_array_shift(other._data, _data,
                              _size[0],  _size[1],
                             offset[0], offset[1]);
            break;
        case 3:
            cuda_array_shift(other._data, _data,
                              _size[0],  _size[1],  _size[2],
                             offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opPad(const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opPad(const cSize3& size, const cSize3& offset)");

    requireNonEmpty(funcName);
    require(size >= _size,          funcName + ": requested size is not larger than current size");
    require(size >= _size + offset, funcName + ": offset value is too large");

    cuArray3<T>    other(size);

    switch (getDimension()) {
        case 1:
            cuda_array_pad(      _data,       _size[0],
                           other._data, other._size[0],
                                             offset[0]);
            break;
        case 2:
            cuda_array_pad(      _data,       _size[0],       _size[1],
                           other._data, other._size[0], other._size[1],
                                             offset[0],      offset[1]);
            break;
        case 3:
            cuda_array_pad(      _data,       _size[0],       _size[1],       _size[2],
                           other._data, other._size[0], other._size[1], other._size[2],
                                             offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opPad(const cuArray3& other, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opPad(const cuArray3& other, const cSize3& size, const cSize3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other,               funcName + ": two objects are the same");
    require(size >= other._size,          funcName + ": other object is too large");
    require(size >= other._size + offset, funcName + ": offset value is too large");

    memReAllocZero(size);

    switch (getDimension()) {
        case 1:
            cuda_array_pad(other._data, other._size[0],
                                 _data,       _size[0],
                                             offset[0]);
            break;
        case 2:
            cuda_array_pad(other._data, other._size[0], other._size[1],
                                 _data,       _size[0],       _size[1],
                                             offset[0],      offset[1]);
            break;
        case 3:
            cuda_array_pad(other._data, other._size[0], other._size[1], other._size[2],
                                 _data,       _size[0],       _size[1],       _size[2],
                                             offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opCrop(const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opCrop(const cSize3& size, const cSize3& offset)");

    requireNonEmpty(funcName);
    require(_size >= size,          funcName + ": requested size is not smaller than current size");
    require(_size >= size + offset, funcName + ": offset value is too large");

    cuArray3<T>    other(size);

    switch (getDimension()) {
        case 1:
            cuda_array_crop(      _data,       _size[0],
                            other._data, other._size[0],
                                              offset[0]);
            break;
        case 2:
            cuda_array_crop(      _data,       _size[0],       _size[1],
                            other._data, other._size[0], other._size[1],
                                              offset[0],      offset[1]);
            break;
        case 3:
            cuda_array_crop(      _data,       _size[0],       _size[1],       _size[2],
                            other._data, other._size[0], other._size[1], other._size[2],
                                              offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opCrop(const cuArray3& other, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opCrop(const cuArray3& other, const cSize3& size, const cSize3& offset)");

    other.requireNonEmpty(funcName);
    require(this != &other,               funcName + ": two objects are the same");
    require(other._size >= size,          funcName + ": other object is too small");
    require(other._size >= size + offset, funcName + ": offset value is too large");

    memReAlloc(size);

    switch (getDimension()) {
        case 1:
            cuda_array_crop(other._data, other._size[0],
                                  _data,       _size[0],
                                              offset[0]);
            break;
        case 2:
            cuda_array_crop(other._data, other._size[0], other._size[1],
                                  _data,       _size[0],       _size[1],
                                              offset[0],      offset[1]);
            break;
        case 3:
            cuda_array_crop(other._data, other._size[0], other._size[1], other._size[2],
                                  _data,       _size[0],       _size[1],       _size[2],
                                              offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opReplace(T value, const cSize3& size, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opReplace(T value, const cSize3& size, const cSize3& offset)");

    requireNonEmpty(funcName);
    require(size <= _size,        funcName + ": replacing region is too large");
    require(size+offset <= _size, funcName + ": offset value is too large");

    switch (getDimension()) {
        case 1:
            cuda_array_replace(_data,  _size[0],
                               value,   size[0],
                                      offset[0]);
            break;
        case 2:
            cuda_array_replace(_data,  _size[0],  _size[1],
                               value,   size[0],   size[1],
                                      offset[0], offset[1]);
            break;
        case 3:
            cuda_array_replace(_data,  _size[0],  _size[1],  _size[2],
                               value,   size[0],   size[1],   size[2],
                                      offset[0], offset[1], offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opReplace(const cuArray3& other, const cSize3& offset)
{
    const std::string   funcName("void cuArray3<T>::opReplace(const cuArray3& other, const cSize3& offset)");

    requireNonEmpty(funcName);
    other.requireNonEmpty(funcName);
    require(this != &other,              funcName + ": two objects are the same");
    require(other._size <= _size,        funcName + ": other object is too large");
    require(other._size+offset <= _size, funcName + ": offset value is too large");

    switch (getDimension()) {
        case 1:
            cuda_array_replace<T,T>(      _data,       _size[0],
                                    other._data, other._size[0],
                                                      offset[0]);
            break;
        case 2:
            cuda_array_replace<T,T>(      _data,       _size[0],       _size[1],
                                    other._data, other._size[0], other._size[1],
                                                      offset[0],      offset[1]);
            break;
        case 3:
            cuda_array_replace<T,T>(      _data,       _size[0],       _size[1],       _size[2],
                                    other._data, other._size[0], other._size[1], other._size[2],
                                                      offset[0],      offset[1],      offset[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuArray3<T>::opPermute(ePermute permute)
{
    const std::string   funcName("void cuArray3<T>::opPermute(ePermute permute)");

    requireNonEmpty(funcName);

    cuArray3<T>    other;

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
            cuda_array_permute(_data, other._data,
                               _size[0], _size[1]);
            break;
        case PERMUTE3D_123:
        case PERMUTE3D_132:
        case PERMUTE3D_213:
        case PERMUTE3D_231:
        case PERMUTE3D_312:
        case PERMUTE3D_321:
            cuda_array_permute(_data, other._data,
                               _size[0], _size[1], _size[2],
                               permute);
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opPermute(const cuArray3& other, ePermute permute)
{
    const std::string   funcName("void cuArray3<T>::opPermute(const cuArray3& other, ePermute permute)");

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
            cuda_array_permute(other._data, _data,
                               other._size[0], other._size[1]);
            break;
        case PERMUTE3D_123:
        case PERMUTE3D_132:
        case PERMUTE3D_213:
        case PERMUTE3D_231:
        case PERMUTE3D_312:
        case PERMUTE3D_321:
            cuda_array_permute(other._data, _data,
                               other._size[0], other._size[1], other._size[2],
                               permute);
            break;
        default:
            ERROR(funcName, "unsupported permutation mode");
    }
}

template <typename T>
void cuArray3<T>::opIndexCol2Row(void)
{
    const std::string   funcName("void cuArray3<T>::opIndexCol2Row(void)");

    requireNonEmpty(funcName);

    cuArray3<T>    other(_size);

    switch (getDimension()) {
        case 2:
            cuda_array_index_col2row(_data, other._data,
                                     _size[0], _size[1]);
            break;
        case 3:
            cuda_array_index_col2row(_data, other._data,
                                     _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::opIndexRow2Col(void)
{
    const std::string   funcName("void cuArray3<T>::opIndexRow2Col(void)");

    requireNonEmpty(funcName);

    cuArray3<T>    other(_size);

    switch (getDimension()) {
        case 2:
            cuda_array_index_row2col(_data, other._data,
                                     _size[0], _size[1]);
            break;
        case 3:
            cuda_array_index_row2col(_data, other._data,
                                     _size[0], _size[1], _size[2]);
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }

    memSwap(other);
}

template <typename T>
void cuArray3<T>::arrayImport(const T* const array, size_t nRow, size_t nCol, size_t nSec)
{
    const std::string   funcName("void cuArray3<T>::arrayImport(const T* const array, size_t nRow, size_t nCol, size_t nSec)");

    if (nRow*nCol*nSec != 0) {
        require(array != NULL, funcName + ": input array cannot be NULL");
    }
    else {
        require(array == NULL, funcName + ": input array must be NULL");
    }

    memReAlloc(cSize3(nRow, nCol, nSec));
    cuda_array_memcpy_d2d(_data, array, getNelement());
}

template <typename T>
void cuArray3<T>::arrayHandle(T* const array, size_t nRow, size_t nCol, size_t nSec)
{
    const std::string   funcName("void cuArray3<T>::arrayHandle(T* const array, size_t nRow, size_t nCol, size_t nSec)");

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
void cuArray3<T>::arrayExport(T* array)
{
    const std::string   funcName("void cuArray3<T>::arrayExport(T* array)");

    require(array == NULL, funcName + ": input array must be NULL");

    cuda_arrayDev_new(array, getNelement());

    cuda_array_memcpy_d2d(array, _data, getNelement());
}

// instantiation
template class cuArray3<int32_t >;
template class cuArray3<uint32_t>;
template class cuArray3<int64_t >;
template class cuArray3<uint64_t>;
template class cuArray3<float   >;
template class cuArray3<double  >;

} // namespace gem
