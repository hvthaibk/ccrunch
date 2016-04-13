/***********************************************************************
 *  File:       registration_TPS.cu
 *
 *  Purpose:    Implementation of registration functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "registration.cuh"

namespace gem {

/**************************
 * Thin-plate splines
 *************************/

// 2D
template <typename T> __global__
void dev_regTPS2_compute_mapping(const T* const Xo,   const T* const Yo,
                                 const T* const Xd,   const T* const Yd,
                                 const T* const X,    const T* const Y,
                                       T* const Xmap,       T* const Ymap,
                                 size_t npoint,
                                 size_t npixel)
{
    size_t      ipixel = blockDim.x * blockIdx.x + threadIdx.x;
    size_t      ipoint;
    T           valTmp, tpsK;

    if (ipixel < npixel) {

        Xmap[ipixel] = Xo[npoint] + X[ipixel]*Xo[npoint+1] + Y[ipixel]*Xo[npoint+2];
        Ymap[ipixel] = Yo[npoint] + X[ipixel]*Yo[npoint+1] + Y[ipixel]*Yo[npoint+2];

        for (ipoint = 0; ipoint < npoint; ipoint ++) {

            valTmp = (Xd[ipoint]-X[ipixel]) * (Xd[ipoint]-X[ipixel]) +
                     (Yd[ipoint]-Y[ipixel]) * (Yd[ipoint]-Y[ipixel]);

            if (valTmp != 0) { tpsK = valTmp * std::log(valTmp); }
            else             { tpsK = 0;                         }

            Xmap[ipixel] += tpsK * Xo[ipoint];
            Ymap[ipixel] += tpsK * Yo[ipoint];
        }
    }
}

template <typename T>
void cuda_regTPS2_compute_mapping(const T* const Xo,   const T* const Yo,
                                  const T* const Xd,   const T* const Yd,
                                  const T* const X,    const T* const Y,
                                        T* const Xmap,       T* const Ymap,
                                  size_t npoint,
                                  size_t npixel)
{
    assert(Xo != NULL && Xd != NULL && X != NULL && Xmap != NULL);
    assert(Yo != NULL && Yd != NULL && Y != NULL && Ymap != NULL);
    assert(npoint > 0 && npixel > 0);

    dev_regTPS2_compute_mapping<<<iDivUp(npixel, BLOCK_1D_NROW), BLOCK_1D_NROW>>>
        (Xo, Yo, Xd, Yd, X, Y, Xmap, Ymap, npoint, npixel);
}

// instantiation
template
void cuda_regTPS2_compute_mapping<float >(const float*  const Xo,   const float*  const Yo,
                                          const float*  const Xd,   const float*  const Yd,
                                          const float*  const X,    const float*  const Y,
                                                float*  const Xmap,       float*  const Ymap,
                                          size_t npoint,
                                          size_t npixel);
template
void cuda_regTPS2_compute_mapping<double>(const double* const Xo,   const double* const Yo,
                                          const double* const Xd,   const double* const Yd,
                                          const double* const X,    const double* const Y,
                                                double* const Xmap,       double* const Ymap,
                                          size_t npoint,
                                          size_t npixel);

} // namespace gem
