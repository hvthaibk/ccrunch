/***********************************************************************
 *  File:       array_mem.cpp
 *
 *  Purpose:    Implementation of array-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "array.hpp"

//#ifndef NDEBUG

#ifdef __GEM_USE_STD11__
    #include <mutex>

namespace gem {

    std::mutex      mutexAlloc;

} // namespace gem

#else
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wconversion"
    #include "boost/thread.hpp"
    #pragma GCC diagnostic pop

namespace gem {

    boost::mutex    mutexAlloc;

} // namespace gem

#endif

namespace gem {

static int      countAlloc1D = 0;
static int      countAlloc2D = 0;
static int      countAlloc3D = 0;

} // namespace gem

//#endif

namespace gem {

/*****************************************
 * memcpy / memset / memcheck
 ****************************************/

// memcpy
template <typename T>
void array_memcpy(T* const arrayDst, const T* const arraySrc, size_t nRow)
{
    assert(arraySrc != NULL);
    assert(arrayDst != NULL);
    assert(nRow > 0);

    memcpy(arrayDst, arraySrc, nRow*sizeof(T));
}

// instantiation
template
void array_memcpy<char    >(char*     const arrayDst, const char*     const arraySrc, size_t nRow);
template
void array_memcpy<int32_t >(int32_t*  const arrayDst, const int32_t*  const arraySrc, size_t nRow);
template
void array_memcpy<uint32_t>(uint32_t* const arrayDst, const uint32_t* const arraySrc, size_t nRow);
template
void array_memcpy<int64_t >(int64_t*  const arrayDst, const int64_t*  const arraySrc, size_t nRow);
template
void array_memcpy<uint64_t>(uint64_t* const arrayDst, const uint64_t* const arraySrc, size_t nRow);
template
void array_memcpy<float   >(float*    const arrayDst, const float*    const arraySrc, size_t nRow);
template
void array_memcpy<double  >(double*   const arrayDst, const double*   const arraySrc, size_t nRow);
template
void array_memcpy<std::complex<float > >
        (std::complex<float >* const arrayDst, const std::complex<float >* const arraySrc, size_t nRow);
template
void array_memcpy<std::complex<double> >
        (std::complex<double>* const arrayDst, const std::complex<double>* const arraySrc, size_t nRow);

// memset
template <typename T>
void array_memset(T* const array, int value, size_t nRow)
{
    assert(array != NULL);
    assert(nRow > 0);

    memset(array, value, nRow*sizeof(T));
}

// instantiation
template
void array_memset<int32_t >(int32_t*  const array, int value, size_t nRow);
template
void array_memset<uint32_t>(uint32_t* const array, int value, size_t nRow);
template
void array_memset<int64_t >(int64_t*  const array, int value, size_t nRow);
template
void array_memset<uint64_t>(uint64_t* const array, int value, size_t nRow);
template
void array_memset<float   >(float*    const array, int value, size_t nRow);
template
void array_memset<double  >(double*   const array, int value, size_t nRow);
template
void array_memset<std::complex<float > >
        (std::complex<float >* const array, int value, size_t nRow);
template
void array_memset<std::complex<double> >
        (std::complex<double>* const array, int value, size_t nRow);

// check memory leak
void array_memcheck(void)
{
    require(countAlloc1D == 0, "array_memcheck: countAlloc1D = " + num2str(countAlloc1D));
    require(countAlloc2D == 0, "array_memcheck: countAlloc2D = " + num2str(countAlloc2D));
    require(countAlloc3D == 0, "array_memcheck: countAlloc3D = " + num2str(countAlloc3D));
}

/*****************************************
 * Array allocation/deallocation
 *       w/o zero initialization
 ****************************************/

// size(array1D) = [nRow]
template <typename T>
void array_new(T*& array1D, size_t nRow)
{
    assert(array1D == NULL);
    assert(nRow > 0);

    try {
        array1D = new T [nRow];
    }
    catch (std::bad_alloc &mem_alloc_exception) {
        std::cerr << "void array_new(T*& array1D, size_t nRow): "
                  << mem_alloc_exception.what()
                  << std::endl;

        if (array1D != NULL) {
            delete [] array1D;
        }

        exit(EXIT_FAILURE);
    }

#ifndef NDEBUG
mutexAlloc.lock();
    countAlloc1D++;
mutexAlloc.unlock();
#endif
}

// instantiation
template
void array_new<bool    >(bool*&     array1D, size_t nRow);
template
void array_new<int8_t  >(int8_t*&   array1D, size_t nRow);
template
void array_new<uint8_t >(uint8_t*&  array1D, size_t nRow);
template
void array_new<int16_t >(int16_t*&  array1D, size_t nRow);
template
void array_new<uint16_t>(uint16_t*& array1D, size_t nRow);
template
void array_new<int32_t >(int32_t*&  array1D, size_t nRow);
template
void array_new<uint32_t>(uint32_t*& array1D, size_t nRow);
template
void array_new<int64_t >(int64_t*&  array1D, size_t nRow);
template
void array_new<uint64_t>(uint64_t*& array1D, size_t nRow);
template
void array_new<float   >(float*&    array1D, size_t nRow);
template
void array_new<double  >(double*&   array1D, size_t nRow);
template
void array_new<float*  >(float**&   array1D, size_t nRow);
template
void array_new<double* >(double**&  array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void array_new_zero(T*& array1D, size_t nRow)
{
    array_new(array1D, nRow);
    array_memset(array1D, 0, nRow);
}

// instantiation
template
void array_new_zero<bool    >(bool*&     array1D, size_t nRow);
template
void array_new_zero<int8_t  >(int8_t*&   array1D, size_t nRow);
template
void array_new_zero<uint8_t >(uint8_t*&  array1D, size_t nRow);
template
void array_new_zero<int16_t >(int16_t*&  array1D, size_t nRow);
template
void array_new_zero<uint16_t>(uint16_t*& array1D, size_t nRow);
template
void array_new_zero<int32_t >(int32_t*&  array1D, size_t nRow);
template
void array_new_zero<uint32_t>(uint32_t*& array1D, size_t nRow);
template
void array_new_zero<int64_t >(int64_t*&  array1D, size_t nRow);
template
void array_new_zero<uint64_t>(uint64_t*& array1D, size_t nRow);
template
void array_new_zero<float   >(float*&    array1D, size_t nRow);
template
void array_new_zero<double  >(double*&   array1D, size_t nRow);
template
void array_new_zero<float*  >(float**&   array1D, size_t nRow);
template
void array_new_zero<double* >(double**&  array1D, size_t nRow);

// size(array1D) = [nRow]
template <typename T>
void array_delete(T*& array1D)
{
    if (array1D != NULL) {
        delete [] array1D;
        array1D = NULL;

#ifndef NDEBUG
mutexAlloc.lock();
        countAlloc1D--;
mutexAlloc.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("array_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void array_delete<bool    >(bool*&     array1D);
template
void array_delete<int8_t  >(int8_t*&   array1D);
template
void array_delete<uint8_t >(uint8_t*&  array1D);
template
void array_delete<int16_t >(int16_t*&  array1D);
template
void array_delete<uint16_t>(uint16_t*& array1D);
template
void array_delete<int32_t >(int32_t*&  array1D);
template
void array_delete<uint32_t>(uint32_t*& array1D);
template
void array_delete<int64_t >(int64_t*&  array1D);
template
void array_delete<uint64_t>(uint64_t*& array1D);
template
void array_delete<float   >(float*&    array1D);
template
void array_delete<double  >(double*&   array1D);
template
void array_delete<float*  >(float**&   array1D);
template
void array_delete<double* >(double**&  array1D);

// size(array2D) = [nRow nCol]
template <typename T>
void array_new(T**& array2D, T*& array1D, size_t nRow, size_t nCol)
{
    assert(array1D == NULL);
    assert(array2D == NULL);
    assert(nRow > 0 && nCol > 0);

    size_t    iRow;

    try {
        array2D = new T *[nRow];
        array2D[0] = new T [nRow*nCol];
    }
    catch (std::bad_alloc &mem_alloc_exception) {
        std::cerr << "void array_new(T**& array2D, T*& array1D, size_t nRow, size_t nCol): "
                  << mem_alloc_exception.what()
                  << std::endl;

        if (array2D[0] != NULL) {
            delete [] array2D[0];
        }
        if (array2D != NULL) {
            delete [] array2D;
        }

        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for
    for (iRow = 0; iRow < nRow; iRow++) {
        array2D[iRow] = array2D[0] + iRow*nCol;
    }

    array1D = array2D[0];

#ifndef NDEBUG
mutexAlloc.lock();
    countAlloc2D++;
mutexAlloc.unlock();
#endif
}

// instantiation
template
void array_new<int32_t >(int32_t**&  array2D, int32_t*&  array1D, size_t nRow, size_t nCol);
template
void array_new<uint32_t>(uint32_t**& array2D, uint32_t*& array1D, size_t nRow, size_t nCol);
template
void array_new<int64_t >(int64_t**&  array2D, int64_t*&  array1D, size_t nRow, size_t nCol);
template
void array_new<uint64_t>(uint64_t**& array2D, uint64_t*& array1D, size_t nRow, size_t nCol);
template
void array_new<float   >(float**&    array2D, float*&    array1D, size_t nRow, size_t nCol);
template
void array_new<double  >(double**&   array2D, double*&   array1D, size_t nRow, size_t nCol);

// size(array2D) = [nRow nCol]
template <typename T>
void array_new_zero(T**& array2D, T*& array1D, size_t nRow, size_t nCol)
{
    array_new(array2D, array1D, nRow, nCol);
    array_memset(array1D, 0, nRow*nCol);
}

// instantiation
template
void array_new_zero<int32_t >(int32_t**&  array2D, int32_t*&  array1D, size_t nRow, size_t nCol);
template
void array_new_zero<uint32_t>(uint32_t**& array2D, uint32_t*& array1D, size_t nRow, size_t nCol);
template
void array_new_zero<int64_t >(int64_t**&  array2D, int64_t*&  array1D, size_t nRow, size_t nCol);
template
void array_new_zero<uint64_t>(uint64_t**& array2D, uint64_t*& array1D, size_t nRow, size_t nCol);
template
void array_new_zero<float   >(float**&    array2D, float*&    array1D, size_t nRow, size_t nCol);
template
void array_new_zero<double  >(double**&   array2D, double*&   array1D, size_t nRow, size_t nCol);

template <typename T>
void array_delete(T**& array2D, T*& array1D)
{
    if (array2D != NULL) {
        delete [] array2D[0];
        delete [] array2D;
        array2D = NULL;
        array1D = NULL;

#ifndef NDEBUG
mutexAlloc.lock();
        countAlloc2D--;
mutexAlloc.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("array_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void array_delete<int32_t >(int32_t**&  array2D, int32_t*&  array1D);
template
void array_delete<uint32_t>(uint32_t**& array2D, uint32_t*& array1D);
template
void array_delete<int64_t >(int64_t**&  array2D, int64_t*&  array1D);
template
void array_delete<uint64_t>(uint64_t**& array2D, uint64_t*& array1D);
template
void array_delete<float   >(float**&    array2D, float*&    array1D);
template
void array_delete<double  >(double**&   array2D, double*&   array1D);

// size(array3D) = [nRow nCol nSec]
template <typename T>
void array_new(T***& array3D, T*& array1D, size_t nRow, size_t nCol, size_t nSec)
{
    assert(array1D == NULL);
    assert(array3D == NULL);
    assert(nRow > 0 && nCol > 0 && nSec > 0);

    size_t    iRow, iCol, indxTmp;

    try {
        array3D = new T **[nRow];

        #pragma omp parallel for
        for (iRow = 0; iRow < nRow; iRow++) {
            array3D[iRow] = new T *[nCol];
        }

        array3D[0][0] = new T [nRow*nCol*nSec];
    }
    catch (std::bad_alloc &mem_alloc_exception) {
        std::cerr << "void array_new(T***& array3D, T*& array1D, size_t nRow, size_t nCol, size_t nSec): "
                  << mem_alloc_exception.what()
                  << std::endl;

        if (array3D[0][0] != NULL) {
            delete [] array3D[0][0];
        }

        #pragma omp parallel for
        for (iRow = 0; iRow < nRow; iRow++) {
            if (array3D[iRow] != NULL) {
                delete [] array3D[iRow];
            }
        }
        if (array3D != NULL) {
            delete [] array3D;
        }

        exit(EXIT_FAILURE);
    }

#ifdef __GEM_USE_OPENMP__
    size_t    nColSec = nCol*nSec;

    #pragma omp parallel for private(iRow,iCol,indxTmp)
    for (iRow = 0; iRow < nRow; iRow++) {
        indxTmp = iRow*nColSec;

        for (iCol = 0; iCol < nCol; iCol++) {
            array3D[iRow][iCol] = array3D[0][0] + indxTmp + iCol*nSec;
        }
    }
#else
    for (iRow = 0, indxTmp = 0; iRow < nRow; iRow++) {
        for (iCol = 0; iCol < nCol; iCol++, indxTmp++) {
            array3D[iRow][iCol] = array3D[0][0] + indxTmp*nSec;
        }
    }
#endif

    array1D = array3D[0][0];

#ifndef NDEBUG
mutexAlloc.lock();
    countAlloc3D++;
mutexAlloc.unlock();
#endif
}

// instantiation
template
void array_new<int32_t >(int32_t***&  array3D, int32_t*&  array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new<uint32_t>(uint32_t***& array3D, uint32_t*& array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new<int64_t >(int64_t***&  array3D, int64_t*&  array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new<uint64_t>(uint64_t***& array3D, uint64_t*& array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new<float   >(float***&    array3D, float*&    array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new<double  >(double***&   array3D, double*&   array1D, size_t nRow, size_t nCol, size_t nSec);

// size(array3D) = [nRow nCol nSec]
template <typename T>
void array_new_zero(T***& array3D, T*& array1D, size_t nRow, size_t nCol, size_t nSec)
{
    array_new(array3D, array1D, nRow, nCol, nSec);
    array_memset(array1D, 0, nRow*nCol*nSec);
}

// instantiation
template
void array_new_zero<int32_t >(int32_t***&  array3D, int32_t*&  array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new_zero<uint32_t>(uint32_t***& array3D, uint32_t*& array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new_zero<int64_t >(int64_t***&  array3D, int64_t*&  array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new_zero<uint64_t>(uint64_t***& array3D, uint64_t*& array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new_zero<float   >(float***&    array3D, float*&    array1D, size_t nRow, size_t nCol, size_t nSec);
template
void array_new_zero<double  >(double***&   array3D, double*&   array1D, size_t nRow, size_t nCol, size_t nSec);

// size(array3D) = [nRow nCol nSec]
template <typename T>
void array_delete(T***& array3D, T*& array1D,size_t nRow)
{
    if (array3D != NULL) {
        delete [] array3D[0][0];

        #pragma omp parallel for
        for (size_t iRow = 0; iRow < nRow; iRow++) {
            delete [] array3D[iRow];
        }

        delete [] array3D;
        array3D = NULL;
        array1D = NULL;

#ifndef NDEBUG
mutexAlloc.lock();
        countAlloc3D--;
mutexAlloc.unlock();
#endif
    }
    else {
#ifndef NDEBUG
        WARNING("array_delete", "attempt to delete a NULL object");
#endif
    }
}

// instantiation
template
void array_delete<int32_t >(int32_t***&  array3D, int32_t*&  array1D, size_t nRow);
template
void array_delete<uint32_t>(uint32_t***& array3D, uint32_t*& array1D, size_t nRow);
template
void array_delete<int64_t >(int64_t***&  array3D, int64_t*&  array1D, size_t nRow);
template
void array_delete<uint64_t>(uint64_t***& array3D, uint64_t*& array1D, size_t nRow);
template
void array_delete<float   >(float***&    array3D, float*&    array1D, size_t nRow);
template
void array_delete<double  >(double***&   array3D, double*&   array1D, size_t nRow);

} // namespace gem
