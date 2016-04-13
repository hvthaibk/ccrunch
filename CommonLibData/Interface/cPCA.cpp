/***********************************************************************
 *  File:       cPCA.cpp
 *
 *  Purpose:    Implementation of PCA classes
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cPCA.hpp"

#include <iostream>
#include <iterator>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#include <Eigen/SVD>
#pragma GCC diagnostic pop

namespace gem {

cPCA2::cPCA2(void)
{
    _nrow = 0;
    _ncol = 0;

    // variables will be scaled by default
    _is_center = true;
    _is_scale  = true;

    // SVD will be used by default
    _method  = "svd";
    _is_corr = false;

    _kaiser   = 0;
    _thresh95 = 1;
}

cPCA2::~cPCA2(void)
{
    _xXf.resize(0,0);
    _x.clear();
}


int cPCA2::computePC(std::vector<float>& x,
                     size_t nrow,
                     size_t ncol,
                     bool is_center,
                     bool is_scale,
                     bool is_corr)
{
    _ncol     = ncol;
    _nrow     = nrow;
    _is_center = is_center;
    _is_scale  = is_scale;
    _is_corr   = is_corr;

    if (x.size() != _nrow*_ncol)     { return -1; }
    if ((1 == _ncol) || (1 == nrow)) { return -1; }

    // convert vector to Eigen 2-dimensional matrix
    _xXf.resize(_nrow, _ncol);

    for (size_t i = 0; i < _nrow; ++i) {
        for (size_t j = 0; j < _ncol; ++j) {
            _xXf(i, j) = x[j + i*_ncol];
        }
    }

    // mean and standard deviation for each column
    Eigen::VectorXf     mean_vector(_ncol),
                        sd_vector(_ncol);
    size_t              zero_sd_num = 0;
    float               denom = static_cast<float>((_nrow > 1) ? _nrow - 1 : 1);

    mean_vector = _xXf.colwise().mean();

    Eigen::VectorXf     curr_col;

    for (size_t i = 0; i < _ncol; ++i) {
        curr_col = Eigen::VectorXf::Constant(_nrow, mean_vector(i));    // mean(x) for column x
        curr_col = _xXf.col(i) - curr_col;                              // x - mean(x)
        curr_col = curr_col.array().square();                           // (x-mean(x))^2

        sd_vector(i) = std::sqrt((curr_col.sum())/denom);

        if (0 == sd_vector(i)) {
            zero_sd_num++;
        }
    }

    // if colums with sd == 0 are too many, don't continue calculation
    if (1 > _ncol-zero_sd_num) {
        return -1;
    }

    // delete columns with sd == 0
    Eigen::MatrixXf     tmp(_nrow, _ncol-zero_sd_num);
    Eigen::VectorXf     tmp_mean_vector(_ncol-zero_sd_num);

    size_t              curr_col_num = 0;

    for (size_t i = 0; i < _ncol; ++i) {
        if (0 != sd_vector(i)) {
            tmp.col(curr_col_num) = _xXf.col(i);
            tmp_mean_vector(curr_col_num) = mean_vector(i);
            curr_col_num++;
        }
        else {
            _eliminated_columns.push_back(i);
        }
    }

    _ncol      -= zero_sd_num;
    _xXf        = tmp;
    mean_vector = tmp_mean_vector;
    tmp.resize(0, 0);
    tmp_mean_vector.resize(0);

    // shift to zero
    if (true == _is_center) {
        for (size_t i = 0; i < _ncol; ++i) {
            _xXf.col(i) -= Eigen::VectorXf::Constant(_nrow, mean_vector(i));
        }
    }

    // scale to unit variance
    if ( (false == _is_corr) || (true == _is_scale)) {
        for (size_t i = 0; i < _ncol; ++i) {
            _xXf.col(i) /= std::sqrt(_xXf.col(i).array().square().sum()/denom);
        }
    }

#ifndef NDEBUG
    std::cout << "\nScaled matrix:\n";
    std::cout << _xXf << std::endl;
    std::cout << "\nMean before scaling:\n" << mean_vector.transpose();
    std::cout << "\nStandard deviation before scaling:\n" << sd_vector.transpose();
#endif

    // when _nrow < _ncol then svd will be used
    // if corr is true and _nrow > _ncol then correlation matrix will be used
    // (TODO): What about covariance?
    if ((_nrow < _ncol) || (false == _is_corr)) {
        _method = "svd";

        Eigen::JacobiSVD<Eigen::MatrixXf>   svd(_xXf, Eigen::ComputeThinV);

        Eigen::VectorXf     eigen_singular_values = svd.singularValues();
        Eigen::VectorXf     tmp_vec = eigen_singular_values.array().square();
        float               tmp_sum = tmp_vec.sum();
        size_t              lim = (_nrow < _ncol)? _nrow : _ncol;

        tmp_vec /= tmp_sum;

        // PC's standard deviation and
        // PC's proportion of variance
        _kaiser = 0;
        for (size_t i = 0; i < lim; ++i) {
            _sd.push_back(eigen_singular_values(i)/std::sqrt(denom));

            if (_sd[i] >= 1) {
                _kaiser = (unsigned int) i + 1;
            }

            _prop_of_var.push_back(tmp_vec(i));
        }

        tmp_vec.resize(0);

#ifndef NDEBUG
        std::cout << "\n\nStandard deviations for PCs:\n";
        copy(_sd.begin(), _sd.end(),std::ostream_iterator<float>(std::cout," "));
        std::cout << "\n\nKaiser criterion: PC #" << _kaiser << std::endl;
#endif

        // PC's cumulative proportion
        _thresh95 = 1;
        _cum_prop.push_back(_prop_of_var[0]);

        for (size_t i = 1; i < _prop_of_var.size(); ++i) {
            _cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);

            if (_cum_prop[i] < 0.95) {
                _thresh95 = (unsigned int) i + 1;
            }
        }

#ifndef NDEBUG
        std::cout << "\nCumulative proportion:\n";
        copy(_cum_prop.begin(), _cum_prop.end(),std::ostream_iterator<float>(std::cout," "));
        std::cout << "\n\nThresh95 criterion: PC #" << _thresh95 << std::endl;
#endif

        // scores
        Eigen::MatrixXf     eigen_scores = _xXf * svd.matrixV();

#ifndef NDEBUG
        std::cout << "\n\nRotated values (scores):\n" << eigen_scores;
#endif

        _scores.reserve(lim*lim);

        for (size_t i = 0; i < lim; ++i) {
            for (size_t j = 0; j < lim; ++j) {
                _scores.push_back(eigen_scores(i, j));
            }
        }

        eigen_scores.resize(0, 0);

#ifndef NDEBUG
        std::cout << "\n\nScores in vector:\n";
        copy(_scores.begin(), _scores.end(),std::ostream_iterator<float>(std::cout," "));
        std::cout << "\n";
#endif

    }
    else {    // COR OR COV MATRICES ARE HERE
        _method = "cor";

        // calculate covariance matrix
        Eigen::MatrixXf     eigen_cov; // = MatrixXf::Zero(_ncol, _ncol);
        Eigen::VectorXf     sds;

        // (TODO) should be weighted cov matrix, even if is_center == false
        eigen_cov = (1.0f /((float) _nrow/*-1*/)) * _xXf.transpose() * _xXf;
        sds = eigen_cov.diagonal().array().sqrt();
        Eigen::MatrixXf outer_sds = sds * sds.transpose();
        eigen_cov = eigen_cov.array() / outer_sds.array();
        outer_sds.resize(0, 0);

        // ?if data matrix is scaled, covariance matrix is equal to correlation matrix
        Eigen::EigenSolver<Eigen::MatrixXf>     edc(eigen_cov);
        Eigen::VectorXf                         eigen_eigenvalues = edc.eigenvalues().real();
        Eigen::MatrixXf                         eigen_eigenvectors = edc.eigenvectors().real();

#ifndef NDEBUG
        std::cout << eigen_cov << std::endl;
        std::cout << std::endl << eigen_eigenvalues.transpose() << std::endl;
        std::cout << std::endl << eigen_eigenvectors << std::endl;
#endif

        // the eigenvalues and eigenvectors are not sorted
        // so, we should sort them
        typedef std::pair<float,int>    eigen_pair;
        std::vector<eigen_pair>         ep;

        for (size_t i = 0 ; i < _ncol; ++i) {
            ep.push_back(std::make_pair(eigen_eigenvalues(i), i));
        }

        sort(ep.begin(), ep.end());     // ascending order by default

        // sort them all in descending order
        Eigen::MatrixXf     eigen_eigenvectors_sorted = Eigen::MatrixXf::Zero(eigen_eigenvectors.rows(), eigen_eigenvectors.cols());
        Eigen::VectorXf     eigen_eigenvalues_sorted  = Eigen::VectorXf::Zero(_ncol);
        int                 colnum = 0;

        for (int i = (int) ep.size()-1; i > -1; i--) {
            eigen_eigenvalues_sorted(colnum)         = ep[i].first;
            eigen_eigenvectors_sorted.col(colnum++) += eigen_eigenvectors.col(ep[i].second);
        }

#ifndef NDEBUG
        std::cout << std::endl << eigen_eigenvalues_sorted.transpose() << std::endl;
        std::cout << std::endl << eigen_eigenvectors_sorted << std::endl;
#endif

        // we don't need not sorted arrays anymore
        eigen_eigenvalues.resize(0);
        eigen_eigenvectors.resize(0, 0);

        _sd.clear();
        _prop_of_var.clear();
        _kaiser = 0;

        float       tmp_sum = eigen_eigenvalues_sorted.sum();

        for (size_t i = 0; i < _ncol; ++i) {
            _sd.push_back(std::sqrt(eigen_eigenvalues_sorted(i)));

            if (_sd[i] >= 1) {
                _kaiser = (unsigned int) i + 1;
            }

            _prop_of_var.push_back(eigen_eigenvalues_sorted(i)/tmp_sum);
        }

#ifndef NDEBUG
        std::cout << "\nStandard deviations for PCs:\n";
        copy(_sd.begin(), _sd.end(), std::ostream_iterator<float>(std::cout," "));
        std::cout << "\nProportion of variance:\n";
        copy(_prop_of_var.begin(), _prop_of_var.end(), std::ostream_iterator<float>(std::cout," "));
        std::cout << "\nKaiser criterion: PC #" << _kaiser << std::endl;
#endif

        // PC's cumulative proportion
        _cum_prop.clear();
        _thresh95 = 1;
        _cum_prop.push_back(_prop_of_var[0]);

        for (size_t i = 1; i < _prop_of_var.size(); ++i) {
            _cum_prop.push_back(_cum_prop[i-1]+_prop_of_var[i]);

            if (_cum_prop[i] < 0.95) {
                _thresh95 = (unsigned int) i + 1;
            }
        }

#ifndef NDEBUG
        std::cout << "\n\nCumulative proportions:\n";
        copy(_cum_prop.begin(), _cum_prop.end(), std::ostream_iterator<float>(std::cout," "));
        std::cout << "\n\n95% threshold: PC #" << _thresh95 << std::endl;
#endif

        // scores for PCA with correlation matrix
        // scale before calculating new values

        for (size_t i = 0; i < _ncol; ++i) {
            _xXf.col(i) /= sds(i);
        }

        sds.resize(0);
        Eigen::MatrixXf     eigen_scores = _xXf * eigen_eigenvectors_sorted;

#ifndef NDEBUG
        std::cout << "\n\nRotated values (scores):\n" << eigen_scores;
#endif

        _scores.clear();
        _scores.reserve(_ncol*_nrow);

        for (size_t i = 0; i < _nrow; ++i) {
            for (size_t j = 0; j < _ncol; ++j) {
                _scores.push_back(eigen_scores(i, j));
            }
        }

        eigen_scores.resize(0, 0);

#ifndef NDEBUG
        std::cout << "\n\nScores in vector:\n";
        copy(_scores.begin(), _scores.end(), std::ostream_iterator<float>(std::cout," "));
        std::cout << "\n";
#endif
    }

    return 0;
}

cPCA3::cPCA3(void)
{
    _nrow = 0;
    _ncol = 0;
    _nsec = 0;

    // variables will be scaled by default
    _is_center = true;
    _is_scale  = true;

    // SVD will be used by default
    _method  = "svd";
    _is_corr = false;

    _kaiser   = 0;
    _thresh95 = 1;
}

cPCA3::~cPCA3(void)
{
    //_xXf.resize(0,0,0);
    _x.clear();
}

/*int cPCA3::computePC(std::vector<float>& x,
                     size_t nrow,
                     size_t ncol,
                     size_t nsec,
                     bool is_center,
                     bool is_scale,
                     bool is_corr);*/

} // namespace gem
