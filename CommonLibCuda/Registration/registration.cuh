/***********************************************************************
 *  File:       registration.cuh
 *
 *  Purpose:    Header file for registration functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_DEV_REGISTRATION_CUH__
#define __GEM_DEV_REGISTRATION_CUH__

#include "array.cuh"

namespace gem {

/**************************
 * Thin-plate splines
 *************************/

// 2D
template <typename T>
void cuda_regTPS2_compute_mapping(const T* const Xo,   const T* const Yo,
                                  const T* const Xd,   const T* const Yd,
                                  const T* const X,    const T* const Y,
                                        T* const Xmap,       T* const Ymap,
                                  size_t npoint,
                                  size_t npixel);

} // namespace gem

#endif
