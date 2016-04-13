/***********************************************************************
 *  File:       mpi_xcorr.cpp
 *
 *  Purpose:    Implementation of MPI-related functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "mpi_xcorr.hpp"

namespace gem {

#ifdef __GEM_USE_MPI__

size_t    resSizeMPI;

template <typename T1, typename T2>
void mpiOp_array_typecast(T1 *srcArray, T2 *dstArray, size_t nrow)
{
    #pragma omp parallel for
    for (size_t i = 0; i < nrow; i++) {
        dstArray[i] = (T2) srcArray[i];
    }
}

template
void mpiOp_array_typecast<size_t,float >(size_t *srcArray, float  *dstArray, size_t nrow);
template
void mpiOp_array_typecast<float ,size_t>(float  *srcArray, size_t *dstArray, size_t nrow);

void mpiOp_xcorrMergeResultGlobalV1Abs(float  *resDataAbsMax,
                                       size_t *resDataMaxInd,
                                       float  *resDataAbsMaxGlobal,
                                       size_t *resDataMaxIndGlobal,
                                       size_t resSize)
{
    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        if (std::abs(resDataAbsMax[i]) > std::abs(resDataAbsMaxGlobal[i])) {
            resDataAbsMaxGlobal[i] = resDataAbsMax[i];
            resDataMaxIndGlobal[i] = resDataMaxInd[i];
        }
    }
}

void mpiOp_xcorrMergeResultGlobalV1(float  *resDataAbsMax,
                                    size_t *resDataMaxInd,
                                    float  *resDataAbsMaxGlobal,
                                    size_t *resDataMaxIndGlobal,
                                    size_t resSize)
{
    #pragma omp parallel for
    for (size_t i = 0; i < resSize; i++) {
        if (resDataAbsMax[i] > resDataAbsMaxGlobal[i]) {
            resDataAbsMaxGlobal[i] = resDataAbsMax[i];
            resDataMaxIndGlobal[i] = resDataMaxInd[i];
        }
    }
}

void mpiOp_xcorrMergeResultGlobalV2Abs(float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype*)
{
    if (*length == 1) {
        #pragma omp parallel for
        for (size_t i = 0; i < resSizeMPI; i++) {
            if (std::abs(invec[i]) > std::abs(inoutvec[i])) {
                inoutvec[i]            = invec[i];               // max
                inoutvec[i+resSizeMPI] = invec[i+resSizeMPI];    // indx
            }
        }
    }
    else {
        std::cout << "mpiOp_xcorrMergeResultGlobalV2Abs(): "
                  << "the vector length seems incorrect"
                  << std::endl;
    }
}

void mpiOp_xcorrMergeResultGlobalV2(float *invec,
                                    float *inoutvec,
                                    int   *length,
                                    MPI_Datatype*)
{
    if (*length == 1) {
        #pragma omp parallel for
        for (size_t i = 0; i < resSizeMPI; i++) {
            if (invec[i] > inoutvec[i]) {
                inoutvec[i]            = invec[i];               // max
                inoutvec[i+resSizeMPI] = invec[i+resSizeMPI];    // indx
            }
        }
    }
    else {
        std::cout << "mpiOp_xcorrMergeResultGlobalV2(): "
                  << "the vector length seems incorrect"
                  << std::endl;
    }
}

void mpiOp_xcorrMergeResultGlobalV3Abs(float *invec,
                                       float *inoutvec,
                                       int   *length,
                                       MPI_Datatype*)
{
    if (*length == (int) resSizeMPI) {
        #pragma omp parallel for
        for (size_t i = 0; i < resSizeMPI; i++) {
            if (std::abs(invec[i]) > std::abs(inoutvec[i])) {
                inoutvec[i]            = invec[i];               // max
                inoutvec[i+resSizeMPI] = invec[i+resSizeMPI];    // indx
            }
        }
    }
    else {
        std::cout << "mpiOp_xcorrMergeResultGlobalV3Abs(): " << resSizeMPI
                  << "the vector length seems incorrect"     << *length
                  << std::endl;
    }
}

void mpiOp_xcorrMergeResultGlobalV3(float *invec,
                                    float *inoutvec,
                                    int   *length,
                                    MPI_Datatype*)
{
    if (*length == (int) resSizeMPI) {
        #pragma omp parallel for
        for (size_t i = 0; i < resSizeMPI; i++) {
            if (invec[i] > inoutvec[i]) {
                inoutvec[i]            = invec[i];               // max
                inoutvec[i+resSizeMPI] = invec[i+resSizeMPI];    // indx
            }
        }
    }
    else {
        std::cout << "mpiOp_xcorrMergeResultGlobalV3(): " << resSizeMPI
                  << "the vector length seems incorrect"  << *length
                  << std::endl;
    }
}

void mpiReduce_pickerV1(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs)
{
    int    mpiID, mpiNumProcs;
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNumProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiID);

    float     *resDataGlobalNode = NULL;
    array_new(resDataGlobalNode, 2*resSize);
    if (mpiID != 0) {
        memcpy(resDataGlobalNode,
               resDataAbsMaxPaddedGlobal,
               resSize*sizeof(float));

        for (size_t i = 0; i < resSize; i++) {
            resDataGlobalNode[i+resSize] =
                        (float) resDataMaxIndPaddedGlobal[i];
        }

        MPI_Send(resDataGlobalNode, (int) (2*resSize),
                 MPI_FLOAT, 0, 123, MPI_COMM_WORLD);

    }
    else {
        MPI_Status    status;
        float         *resDataAbsMaxPaddedGlobalNode = NULL;
        size_t        *resDataMaxIndPaddedGlobalNode = NULL;

        array_new(resDataAbsMaxPaddedGlobalNode, resSize);
        array_new(resDataMaxIndPaddedGlobalNode, resSize);

        for (int n = 1; n < mpiNumProcs; n++){
            MPI_Recv(resDataGlobalNode, (int) (2*resSize),
                     MPI_FLOAT, n, 123, MPI_COMM_WORLD, &status);

            memcpy(resDataAbsMaxPaddedGlobalNode,
                   resDataGlobalNode,
                   resSize*sizeof(float));

            for (size_t i = 0; i < resSize; i++) {
                resDataMaxIndPaddedGlobalNode[i] =
                            (size_t) resDataGlobalNode[i+resSize];
            }

            switch (bAbs) {
                case XCORR_MERGE_NEGATIVE:
                    mpiOp_xcorrMergeResultGlobalV1Abs(resDataAbsMaxPaddedGlobalNode,
                                                      resDataMaxIndPaddedGlobalNode,
                                                      resDataAbsMaxPaddedGlobal,
                                                      resDataMaxIndPaddedGlobal,
                                                      resSize);
                    break;
                case XCORR_MERGE_POSITIVE:
                    mpiOp_xcorrMergeResultGlobalV1(resDataAbsMaxPaddedGlobalNode,
                                                   resDataMaxIndPaddedGlobalNode,
                                                   resDataAbsMaxPaddedGlobal,
                                                   resDataMaxIndPaddedGlobal,
                                                   resSize);
                    break;
                default:
                    ERROR("mpiReduce_pickerV1", "unsupported merging mode");
            }
        }

        array_delete(resDataAbsMaxPaddedGlobalNode);
        array_delete(resDataMaxIndPaddedGlobalNode);
    }
    array_delete(resDataGlobalNode);
}

void mpiReduce_pickerV2(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs)
{
    resSizeMPI = resSize;

    MPI_Datatype mpiType;
    MPI_Type_contiguous((int) (2*resSize), MPI_FLOAT, &mpiType);
    MPI_Type_commit(&mpiType);

    float     *resDataGlobalNode = NULL;
    float     *resDataGlobalNodeReduce = NULL;

    array_new(resDataGlobalNode, 2*resSize);
    array_new(resDataGlobalNodeReduce, 2*resSize);

    memcpy(resDataGlobalNode,
           resDataAbsMaxPaddedGlobal,
           resSize*sizeof(float));
    mpiOp_array_typecast(resDataMaxIndPaddedGlobal,
                         resDataGlobalNode+resSize,
                         resSize);

    MPI_Op mpiOp;

    switch (bAbs) {
        case XCORR_MERGE_NEGATIVE:
            MPI_Op_create((MPI_User_function *) mpiOp_xcorrMergeResultGlobalV2Abs,
                          1,            // commutative
                          &mpiOp);
            break;
        case XCORR_MERGE_POSITIVE:
            MPI_Op_create((MPI_User_function *) mpiOp_xcorrMergeResultGlobalV2,
                          1,            // commutative
                          &mpiOp);
            break;
        default:
            ERROR("mpiReduce_pickerV2", "unsupported merging mode");
    }

    MPI_Reduce(resDataGlobalNode,
               resDataGlobalNodeReduce,
               1,                   // 1 element of size 2*resSize*sizeof(float)
               mpiType,
               mpiOp,
               0,
               MPI_COMM_WORLD);
    MPI_Op_free(&mpiOp);

    memcpy(resDataAbsMaxPaddedGlobal,
           resDataGlobalNodeReduce,
           resSize*sizeof(float));
    mpiOp_array_typecast(resDataGlobalNodeReduce+resSize,
                         resDataMaxIndPaddedGlobal,
                         resSize);

    array_delete(resDataGlobalNode);
    array_delete(resDataGlobalNodeReduce);
    MPI_Type_free(&mpiType);
}

void mpiReduce_pickerV3(float       *resDataAbsMaxPaddedGlobal,
                        size_t      *resDataMaxIndPaddedGlobal,
                        size_t      resSize,
                        eXCorrMerge bAbs)
{
    resSizeMPI = resSize;

    MPI_Datatype mpiType;
    MPI_Type_contiguous((int) 2, MPI_FLOAT, &mpiType);
    MPI_Type_commit(&mpiType);

    float     *resDataGlobalNode = NULL;
    float     *resDataGlobalNodeReduce = NULL;

    array_new(resDataGlobalNode, 2*resSize);
    array_new(resDataGlobalNodeReduce, 2*resSize);

    memcpy(resDataGlobalNode,
           resDataAbsMaxPaddedGlobal,
           resSize*sizeof(float));
    mpiOp_array_typecast(resDataMaxIndPaddedGlobal,
                         resDataGlobalNode+resSize,
                         resSize);

    MPI_Op mpiOp;

    switch (bAbs) {
        case XCORR_MERGE_NEGATIVE:
            MPI_Op_create((MPI_User_function *) mpiOp_xcorrMergeResultGlobalV3Abs,
                          1,            // commutative
                          &mpiOp);
            break;
        case XCORR_MERGE_POSITIVE:
            MPI_Op_create((MPI_User_function *) mpiOp_xcorrMergeResultGlobalV3,
                          1,            // commutative
                          &mpiOp);
            break;
        default:
            ERROR("mpiReduce_pickerV3", "unsupported merging mode");
    }
    MPI_Reduce(resDataGlobalNode,
               resDataGlobalNodeReduce,
               (int) resSize,       // resSize elements of size 2*sizeof(float)
               mpiType,
               mpiOp,
               0,
               MPI_COMM_WORLD);
    MPI_Op_free(&mpiOp);

    memcpy(resDataAbsMaxPaddedGlobal,
           resDataGlobalNodeReduce,
           resSize*sizeof(float));
    mpiOp_array_typecast(resDataGlobalNodeReduce+resSize,
                         resDataMaxIndPaddedGlobal,
                         resSize);

    array_delete(resDataGlobalNode);
    array_delete(resDataGlobalNodeReduce);
    MPI_Type_free(&mpiType);
}

#endif  // __GEM_USE_MPI__

} // namespace gem
