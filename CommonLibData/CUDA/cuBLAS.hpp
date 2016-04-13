/***********************************************************************
 *  File:       cuBLAS.hpp
 *
 *  Purpose:    Header file for a cuBLAS-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CUBLAS_HPP__
#define __GEM_CUBLAS_HPP__

#include "cuData3.hpp"

#include <cublas_v2.h>

namespace gem {

/* Incompatibility between different CUBLAS versions. This error occurs
 * with CUDA 6.5
 *
 *      invalid conversion from float** to const float**
 *
 * and the solution is described in
 *
 *      http://stackoverflow.com/questions/16371043/why-does-cublas-use-const-pointers-for-parameters
 *
 * But the above error does not exist in CUDA 5.5 and
 * this solution does not work for CUDA 5.5
 */

template <typename T>
class cuBLAS
{
private:
    cublasHandle_t      handle;
    cublasStatus_t      status;

public:
     cuBLAS(void) ;
    ~cuBLAS()     ;

    // data conversion
    void    convertDataCol2Row(const cuData3<T>& a, cuData3<T>& b);
    void    convertDataRow2Col(const cuData3<T>& a, cuData3<T>& b);

    // linear algebra
    void    inverseMat (const cuData3<T>& A, cuData3<T>& B);
    void    multiplyMat(const cuData3<T>& A, const cuData3<T>& B, cuData3<T>& C,
                        unsigned int mode = 0);
    void    solveAxb   (const cuData3<T>& A, const cuData3<T>& B, cuData3<T>& X,
                        unsigned int mode = 0);
};

} // namespace gem

#endif
