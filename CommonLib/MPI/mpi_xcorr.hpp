/***********************************************************************
 *  File:       mpi_xcorr.hpp
 *
 *  Purpose:    Header file for MPI-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_LIB_MPI_HPP__
#define __GEM_LIB_MPI_HPP__

#include "config.hpp"

#ifdef __GEM_USE_MPI__

#include "array.hpp"
#include "xcorr_helper.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <mpi.h>
#pragma GCC diagnostic pop

namespace gem {

template <typename T1, typename T2>
void mpiOp_array_typecast(T1 *srcArray, T2 *dstArray, size_t nrow);

// mpiReduce_pickerV1()
void mpiOp_xcorrMergeResultGlobalV1Abs(float  *resDataAbsMax,
                                       size_t *resDataMaxInd,
                                       float  *resDataAbsMaxGlobal,
                                       size_t *resDataMaxIndGlobal,
                                       size_t resSize);
void mpiOp_xcorrMergeResultGlobalV1   (float  *resDataAbsMax,
                                       size_t *resDataMaxInd,
                                       float  *resDataAbsMaxGlobal,
                                       size_t *resDataMaxIndGlobal,
                                       size_t resSize);
// mpiReduce_pickerV2()
void mpiOp_xcorrMergeResultGlobalV2Abs(float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype *dtype);
void mpiOp_xcorrMergeResultGlobalV2   (float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype *dtype);
// mpiReduce_pickerV3()
void mpiOp_xcorrMergeResultGlobalV3Abs(float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype *dtype);
void mpiOp_xcorrMergeResultGlobalV3   (float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype *dtype);

// MPI_Send() + MPI_Recv()
void mpiReduce_pickerV1(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs);

// MPI_Reduce() with 1 element of data of size 2*resSize*sizeof(float)
void mpiReduce_pickerV2(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs);

// MPI_Reduce() with resSize elements of data of size 2*sizeof(float)
// This function produces incorrect results when resSize > 4096. It is
// observed that when the parameter count of MPI_Reduce() > 8192, the
// the data in the buffer is sent in several chunks, each has less than
// 8192 elements.
void mpiReduce_pickerV3(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs);

} // namespace gem

#endif  // __GEM_USE_MPI__

#endif  // __GEM_LIB_MPI_HPP__
