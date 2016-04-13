/***********************************************************************
 *  File:       cArmadillo.cpp
 *
 *  Purpose:    Implementation of an Armadillo-interface class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cArmadillo.hpp"

#ifdef __GEM_USE_ARMADILLO__

namespace gem {

template <typename T>
void cArmadillo<T>::checkExtension(const std::string& fileName) const
{
    const std::string   funcName("void cArmadillo<T>::checkExtension("
                                    "const std::string& fileName) const");

    if (fileName == "" || (fileName.find(".mat") > fileName.length() &&
                           fileName.find(".MAT") > fileName.length())) {
        ERROR(funcName, "unsupported expansion (accept *.mat, *.MAT)");
    }
}

template <typename T>
void cArmadillo<T>::convertData(const cData3<T>& a, arma::Mat<T>& b)
{
    const std::string   funcName("void cArmadillo<T>::convertData("
                                    "const cData3<T>& a, arma::Mat<T>& b)");

    a.requireNonEmpty(funcName);
    require(a.getNsec() == 1, funcName + ": a is not a matrix");

    b.resize(a.getNrow(),a.getNcol());

    for (size_t i = 0, k = 0; i < a.getNrow(); i++) {
        for (size_t j = 0; j < a.getNcol(); j++, k++) {
            b(i,j) = a[k];
        }
    }
}

template <typename T>
void cArmadillo<T>::convertData(const arma::Mat<T>& a, cData3<T>& b)
{
    const std::string   funcName("void cArmadillo<T>::convertData("
                                    "const arma::Mat<T>& a, cData3<T>& b)");

    require(a.size() > 0, funcName + ": a is not a matrix");

    b.memReAlloc(cSize3(a.n_rows,a.n_cols,1));

    for (size_t i = 0, k = 0; i < b.getNrow(); i++) {
        for (size_t j = 0; j < b.getNcol(); j++, k++) {
            b[k] = a(i,j);
        }
    }
}

template <typename T>
void cArmadillo<T>::load(const std::string& fileName, cData3<T>& data)
{
    const std::string   funcName("void cArmadillo<T>::load("
                                    "const std::string& fileName, "
                                    "cData3<T>& data)");

    checkExtension(fileName);

    arma::Mat<T>    _data;

    require(_data.load(fileName),
            funcName + ": cannot load MATLAB file " + fileName);

    convertData(_data, data);

    require(data.getNelement() > 0, funcName + ": reading data is empty");
}

template <typename T>
void cArmadillo<T>::save(const cData3<T>& data, const std::string& fileName)
{
    const std::string   funcName("void cArmadillo<T>::save("
                                    "const cData3<T>& data, "
                                    "const std::string& fileName)");

    checkExtension(fileName);

    require(data.getNelement() > 0, funcName + ": writing data is empty");

    arma::Mat<T>    _data;

    convertData(data, _data);

    require(_data.save(fileName, arma::raw_ascii),
            funcName + ": cannot save MATLAB file " + fileName);
}

template <typename T>
void cArmadillo<T>::inverseMat(const cData3<T>& A, cData3<T>& B)
{
    const std::string   funcName("void cArmadillo<T>::inverseMat("
                                    "const cData3<T>& A, "
                                    "cData3<T>& B)");

    A.requireNonEmpty(funcName);
    require(A.getNrow() == A.getNcol() && A.getNsec() == 1,
            funcName + ": A is not a square matrix");

    arma::Mat<T>    _A, _B;

    convertData(A,_A);

    _B = arma::inv(_A);

    convertData(_B,B);
}

template <typename T>
void cArmadillo<T>::multiplyMat(const cData3<T>& A, const cData3<T>& B, cData3<T>& C)
{
    const std::string   funcName("void cArmadillo<T>::multiplyMat("
                                    "const cData3<T>& A, "
                                    "const cData3<T>& B, "
                                    "cData3<T>& C)");

    require(A.getDimension() == 2 && A.getNsec() == 1,
            funcName + ": A is not a matrix");
    require(B.getDimension() == 2 && B.getNsec() == 1,
            funcName + ": B is not a matrix");
    require(A.getNcol() == B.getNrow(),
            funcName + ": A and B are not compatible");

    arma::Mat<T>    _A, _B, _C;

    convertData(A,_A);
    convertData(B,_B);

    _C = _A * _B;

    convertData(_C,C);
}

template <typename T>
void cArmadillo<T>::solveAxb(const cData3<T>& A, const cData3<T>& B, cData3<T>& X)
{
    const std::string   funcName("void cArmadillo<T>::solveAxb("
                                    "const cData3<T>& A, "
                                    "const cData3<T>& B, "
                                    "cData3<T>& X)");

    require(A.getDimension() == 2 && A.getNsec() == 1,
            funcName + ": A is not a matrix");
    require(A.getNcol() == B.getNrow(),
            funcName + ": A and B are not compatible");

    arma::Mat<T>    _A, _B, _X;

    convertData(A,_A);
    convertData(B,_B);

    _X = arma::solve(_A,_B);

    convertData(_X,X);
}

// instantiation
template class cArmadillo<float >;
template class cArmadillo<double>;

} // namespace gem

#endif  // __GEM_USE_ARMADILLO__
