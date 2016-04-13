/***********************************************************************
 *  File:       cuBLAS.cpp
 *
 *  Purpose:    Implementation of a cuBLAS-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cuBLAS.hpp"

namespace gem {

template <typename T>
cuBLAS<T>::cuBLAS(void)
{
    const std::string   funcName("cuBLAS<T>::cuBLAS(void)");

    if (cublasCreate(&handle) != CUBLAS_STATUS_SUCCESS) {
        cudaDeviceReset();
        ERROR(funcName, "CUBLAS initialization failed");
    }
};

template <typename T>
cuBLAS<T>::~cuBLAS()
{
    cublasDestroy(handle);
};

template <typename T>
void cuBLAS<T>::convertDataCol2Row(const cuData3<T>& a, cuData3<T>& b)
{
    const std::string   funcName("void cuBLAS<T>::convertDataCol2Row("
                                    "const cuData3<T>& a, cuData3<T>& b)");

    a.requireNonEmpty(funcName);
    b.memReAlloc(a.getSize());

    switch (a.getDimension()) {
        case 2:
            cuda_array_index_col2row(a.getAddrData(), b.getAddrData(),
                                     a.getNrow(), a.getNcol());
            break;
        case 3:
            cuda_array_index_col2row(a.getAddrData(), b.getAddrData(),
                                     a.getNrow(), a.getNcol(), a.getNsec());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <typename T>
void cuBLAS<T>::convertDataRow2Col(const cuData3<T>& a, cuData3<T>& b)
{
    const std::string   funcName("void cuBLAS<T>::convertDataRow2Col("
                                    "const cuData3<T>& a, cuData3<T>& b)");

    a.requireNonEmpty(funcName);
    b.memReAlloc(a.getSize());

    switch (a.getDimension()) {
        case 2:
            cuda_array_index_row2col(a.getAddrData(), b.getAddrData(),
                                     a.getNrow(), a.getNcol());
            break;
        case 3:
            cuda_array_index_row2col(a.getAddrData(), b.getAddrData(),
                                     a.getNrow(), a.getNcol(), a.getNsec());
            break;
        default:
            ERROR(funcName, "unsupported dimension");
    }
}

template <>
void cuBLAS<float>::inverseMat(const cuData3<float>& A,
                                     cuData3<float>& B)
{
    const std::string   funcName("void cuBLAS<float>::inverseMat("
                                    "const cuData3<float>& A, "
                                    "cuData3<float>& B)");

    A.requireNonEmpty(funcName);
    require(A.getNrow() == A.getNcol() && A.getNsec() == 1,
            funcName + ": A is not a square matrix");

    int                 n = (int) A.getNrow();
    cuArray3<int>       pivot(cSize3(n,1,1)), info(cSize3(1,1,1));

    cuData3<float>      _A;
    float**             Acpu = NULL;
    float**             Agpu = NULL;

    cuData3<float>      _B;
    float**             Bcpu = NULL;
    float**             Bgpu = NULL;

    // convert data
    convertDataRow2Col(A, _A);

    // allocation
    _B.memReAlloc(A.getSize());

    array_new(Acpu, 1);
    array_new(Bcpu, 1);
    cuda_arrayDev_new(Agpu, 1);
    cuda_arrayDev_new(Bgpu, 1);

    Acpu[0] = _A.getAddrData();
    Bcpu[0] = _B.getAddrData();
    CUDA_SAFE_CALL( cudaMemcpy(Agpu, Acpu, sizeof(float*),
                               cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(Bgpu, Bcpu, sizeof(float*),
                               cudaMemcpyHostToDevice) );

    // LU decomposition
    status = cublasSgetrfBatched(handle, n, Agpu, n,
                pivot.getAddrData(), info.getAddrData(), 1);

    require(status == CUBLAS_STATUS_SUCCESS,
            funcName + ": LU decomposition error");

/*    // matrix inverion
    status = cublasSgetriBatched(handle, n, (const float**) Agpu, n,
                pivot.getAddrData(), Bgpu, n, info.getAddrData(), 1);
*/
    require(status == CUBLAS_STATUS_SUCCESS,
            funcName + ": matrix inversion error");

    // convert data
    convertDataCol2Row(_B, B);

    // deallocation
    array_delete(Acpu);
    array_delete(Bcpu);
    cuda_arrayDev_delete(Agpu);
    cuda_arrayDev_delete(Bgpu);

       _A.memFree();
       _B.memFree();
    pivot.memFree();
     info.memFree();
}

template <>
void cuBLAS<double>::inverseMat(const cuData3<double>& A,
                                      cuData3<double>& B)
{
    const std::string   funcName("void cuBLAS<double>::inverseMat("
                                    "const cuData3<double>& A, "
                                    "cuData3<double>& B)");

    A.requireNonEmpty(funcName);
    require(A.getNrow() == A.getNcol() && A.getNsec() == 1,
            funcName + ": A is not a square matrix");

    int                 n = (int) A.getNrow();
    cuArray3<int>       pivot(cSize3(n,1,1)), info(cSize3(1,1,1));

    cuData3<double>     _A;
    double**            Acpu = NULL;
    double**            Agpu = NULL;

    cuData3<double>     _B;
    double**            Bcpu = NULL;
    double**            Bgpu = NULL;

    // convert data
    convertDataRow2Col(A, _A);

    // allocation
    _B.memReAlloc(A.getSize());

    array_new(Acpu, 1);
    array_new(Bcpu, 1);
    cuda_arrayDev_new(Agpu, 1);
    cuda_arrayDev_new(Bgpu, 1);

    Acpu[0] = _A.getAddrData();
    Bcpu[0] = _B.getAddrData();
    CUDA_SAFE_CALL( cudaMemcpy(Agpu, Acpu, sizeof(double*),
                               cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL( cudaMemcpy(Bgpu, Bcpu, sizeof(double*),
                               cudaMemcpyHostToDevice) );

    // LU decomposition
    status = cublasDgetrfBatched(handle, n, Agpu, n,
                pivot.getAddrData(), info.getAddrData(), 1);

    require(status == CUBLAS_STATUS_SUCCESS,
            funcName + ": LU decomposition error");

/*    // matrix inverion
    status = cublasDgetriBatched(handle, n, (const double**) Agpu, n,
                pivot.getAddrData(), Bgpu, n, info.getAddrData(), 1);
*/
    require(status == CUBLAS_STATUS_SUCCESS,
            funcName + ": matrix inversion error");

    // convert data
    convertDataCol2Row(_B, B);

    // deallocation
    array_delete(Acpu);
    array_delete(Bcpu);
    cuda_arrayDev_delete(Agpu);
    cuda_arrayDev_delete(Bgpu);

       _A.memFree();
       _B.memFree();
    pivot.memFree();
     info.memFree();
}

template <typename T>
void cuBLAS<T>::multiplyMat(const cuData3<T>& A,
                            const cuData3<T>& B,
                                  cuData3<T>& C,
                            unsigned int mode)
{
    const std::string   funcName("void cuBLAS<float>::multiplyMat("
                                    "const cuData3<T>& A, "
                                    "const cuData3<T>& B, "
                                    "cuData3<T>& C), "
                                    "unsigned int mode");

    require(A.getDimension() == 2 && A.getNsec() == 1,
            funcName + ": A is not a matrix");
    require(B.getDimension() == 2 && B.getNsec() == 1,
            funcName + ": B is not a matrix");
    require(A.getNcol() == B.getNrow(),
            funcName + ": A and B are not compatible");

    C.memReAlloc(cSize3(A.getNrow(),B.getNcol(),1));

    cuda_array_blas_gemm(C.getAddrData(),
                         A.getAddrData(),
                         B.getAddrData(),
                         A.getNrow(), A.getNcol(), B.getNcol(),
                         mode);
}

template <typename T>
void cuBLAS<T>::solveAxb(const cuData3<T>& A,
                         const cuData3<T>& B,
                               cuData3<T>& X,
                         unsigned int mode)
{
    const std::string   funcName("void cuBLAS<T>::solveAxb("
                                    "const cuData3<T>& A, "
                                    "const cuData3<T>& B, "
                                    "cuData3<T>& X, "
                                    "unsigned int mode)");

    cuData3<T>      _Ainv;

    inverseMat(A, _Ainv);

    multiplyMat(_Ainv, B, X, mode);
}

// instantiation
template class cuBLAS<float >;
template class cuBLAS<double>;

} // namespace gem
