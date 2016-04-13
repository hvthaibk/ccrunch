/***********************************************************************
 *  File:       cRegTrackPoint2.cpp
 *
 *  Purpose:    Implementation of a 2D point tracking class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cRegTrackPoint2.hpp"

#include "cNDgrid.hpp"
#include "cRegCorr.hpp"

namespace gem {

template <typename T>
void cRegTrackPoint2<T>::memFree(void)
{
     _ptsOriginal.memFree();
    _ptsDistorted.memFree();
       _corrValue.memFree();
}

template <typename T>
cData3<T>& cRegTrackPoint2<T>::getPtsOriginal(void)
{
    const std::string   funcName("cData3<T>& cRegTrackPoint2<T>::"
                                    "getPtsOriginal(void)");

    _ptsOriginal.requireNonEmpty(funcName);

    return _ptsOriginal;
}

template <typename T>
cData3<T>& cRegTrackPoint2<T>::getPtsDistorted(void)
{
    const std::string   funcName("cData3<T>& cRegTrackPoint2<T>::"
                                    "getPtsDistorted(void)");

    _ptsDistorted.requireNonEmpty(funcName);

    return _ptsDistorted;
}

template <typename T>
cData3<T>& cRegTrackPoint2<T>::getCorrValues(void)
{
    const std::string   funcName("cData3<T>& cRegTrackPoint2<T>::"
                                    "getCorrValues(void)");

    _corrValue.requireNonEmpty(funcName);

    return _corrValue;
}

template <typename T>
void cRegTrackPoint2<T>::setPointsOriginal(const cData3<T>& points)
{
    const std::string   funcName("void cRegTrackPoint2<T>::setPointsOriginal("
                                    "const cData3<T>& points)");

    points.requireNonEmpty(funcName);
    require(points.getNrow() > 1 && points.getNcol() == 2,
            funcName + ": input data is not 2D points");

    _ptsOriginal = points;
}

template <typename T>
void cRegTrackPoint2<T>::genPointsRectangular(const cSize3& size,
                                              const cSize3& border,
                                              const cSize3& step)
{
    const std::string   funcName("void cRegTrackPoint2<T>::genPointsRectangular("
                                    "const cSize3& size, "
                                    "const cSize3& border, "
                                    "const cSize3& step)");

    require(size[0] > 1 && size[1] > 1 && size[2] == 1, funcName + ": size is not 2D");

    cData3<T>   rVec, cVec, R, C;

    rVec.linspace((T) (border[0] - 1), (T) (size[0] - border[0] - 1),
                  (size_t) std::floor((T) (size[0] - 2*border[0]) / (T) step[0]));
    cVec.linspace((T) (border[1] - 1), (T) (size[1] - border[1] - 1),
                  (size_t) std::floor((T) (size[1] - 2*border[1]) / (T) step[1]));
    rVec.round();
    cVec.round();

    cNDgrid<T>().genRC(rVec, cVec, R, C);

    require(R.getNelement() > 0, funcName + ": point set is empty");

    _ptsOriginal.memReAlloc(cSize3(R.getNelement(),2,1));

    #pragma omp parallel for
    for (size_t i = 0; i < R.getNelement(); i++) {
        _ptsOriginal[2*i]   = R[i];
        _ptsOriginal[2*i+1] = C[i];
    }
}

template <typename T>
void cRegTrackPoint2<T>::trackPointsUsingCorr(const cData3<T>& dataOriginal,
                                              const cData3<T>& dataDistorted,
                                              const cSize3&    tplSizeHalf,
                                              const cSize3&    searchRange,
                                              eCCmode mode)
{
    const std::string   funcName("void cRegTrackPoint2<T>::trackPointsUsingCorr("
                                    "const cData3<T>& dataOriginal, "
                                    "const cData3<T>& dataDistorted, "
                                    "const cSize3& tplSizeHalf, "
                                    "eCCmode mode");

    require(dataOriginal.getSize() == dataDistorted.getSize(),
            funcName + ": input images do not have the same size");

    cRegCorr<T>         objCorr;
    cData3XCorr<T>      dataTpl, dataRef, dataRes;
    cSize3              tplSize(2*tplSizeHalf[0]+1,
                                2*tplSizeHalf[1]+1,
                                1);
    cSize3              refSize(2*tplSizeHalf[0]+2*searchRange[0]+1,
                                2*tplSizeHalf[1]+2*searchRange[1]+1,
                                1);
    cSize3              tplOffset(0,0,0), refOffset(0,0,0), resDelta(0,0,0);

    _ptsDistorted.memReAlloc(_ptsOriginal.getSize());
       _corrValue.memReAlloc(cSize3(_ptsOriginal.getNrow(),1,1));

    #pragma omp parallel for private(dataTpl,dataRef,dataRes,tplOffset,refOffset,resDelta)
    for (size_t i = 0; i < _ptsOriginal.getNrow(); i++) {

        tplOffset[0] = (size_t) _ptsOriginal[2*i]   - tplSizeHalf[0];
        tplOffset[1] = (size_t) _ptsOriginal[2*i+1] - tplSizeHalf[1];
        refOffset[0] = (size_t) _ptsOriginal[2*i]   - tplSizeHalf[0] - searchRange[0];
        refOffset[1] = (size_t) _ptsOriginal[2*i+1] - tplSizeHalf[1] - searchRange[1];

        dataTpl.opCrop(dataOriginal,  tplSize, tplOffset);
        dataRef.opCrop(dataDistorted, refSize, refOffset);

        switch (mode) {
            case CC_MODE_ECORR:
                objCorr.computeECC(dataTpl, dataRef, dataRes, XCORR_RES_VALID, false);
                break;
            case CC_MODE_NCORR:
                objCorr.computeNCC(dataTpl, dataRef, dataRes, XCORR_RES_VALID, false);

                dataRes.opLimitCorrValue((T) -1.05, (T) 1.05);
                //require(dataRes.isFinite(), funcName + ": invalid correlation results");
                break;
            default:
                ERROR(funcName, "unsupported correlation mode");
        }

        resDelta.ind2sub(dataRes.getSize(), dataRes.getMaxIndx());

        _ptsDistorted[2*i]   = _ptsOriginal[2*i]   - (T) resDelta[0] + (T) searchRange[0];
        _ptsDistorted[2*i+1] = _ptsOriginal[2*i+1] - (T) resDelta[1] + (T) searchRange[1];

        _corrValue[i] = dataRes.getMax();
    }
}

// instantiation
template class cRegTrackPoint2<float >;
template class cRegTrackPoint2<double>;

} // namespace gem
