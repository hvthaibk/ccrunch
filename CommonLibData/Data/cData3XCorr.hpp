/***********************************************************************
 *  File:       cData3XCorr.hpp
 *
 *  Purpose:    Header file for a correlation-data class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CDATA3_XCORR_HPP__
#define __GEM_CDATA3_XCORR_HPP__

#include "cData3.hpp"

namespace gem {

template <typename T>
class cData3XCorr : public cData3<T>
{
private:
    cData3<T>       _mask, _wghtTpl, _wghtRef;
    size_t          _id;
    std::string     _fileName;
    size_t          _mskArea;

public:
     cData3XCorr(const cSize3& size = cSize3(0,0,0));
     cData3XCorr(const cData3XCorr& other);
    ~cData3XCorr() {};

    void        memFree(void) { cArray3<T>::memFree(); _mask.memFree(); _wghtTpl.memFree(); _wghtRef.memFree(); }

    // ID
    void        setID(size_t id)  { _id = id;   }
    size_t      getID(void) const { return _id; }

    // file name
    void        setFileName  (const std::string& fileName) { _fileName = fileName;                                     }
    std::string getFileName  (void) const                  { return _fileName;                                         }
    void        printFileName(void) const                  { std::cout << "data source  = " << _fileName << std::endl; }

    // operation
    void        opFFTSizeChange (      cSize3& sizeOri, const cSize3& sizeMin);
    void        opFFTSizeRestore(const cSize3& sizeOri);

    void        opNormalizeTplECC (bool useMask = true);
    void        opNormalizeTplNCC (bool useMask = true);
    void        opNormalizeTplWCC1(bool useMask = true);
    void        opNormalizeTplWCC2(bool useMask = true);

    void        opLimitCorrValue   (T valMin, T valMax);
    void        opClearBorderMargin(size_t marginSize);
    void        opReplaceUsingMask (const cData3XCorr& mask, T value, const cVector3<ptrdiff_t>& offset, bool maskInside = true);
    void        opCopyUsingMask    (const cData3XCorr& mask, cData3XCorr& output, const cVector3<ptrdiff_t>& offset);

    // mask and weight
    cData3<T>           maskGet        (void)       { return  _mask;               }
    cData3<T>*          maskGetRef     (void)       { return &_mask;               }
    T*                  maskGetAddrData(void) const { return  _mask.getAddrData(); }
    cSize3              maskGetSize    (void) const { return  _mask.getSize();     }
    size_t              maskGetArea    (void) const { return  _mskArea;            }

    cData3<T>           wghtTplGet        (void)       { return  _wghtTpl;               }
    cData3<T>*          wghtTplGetRef     (void)       { return &_wghtTpl;               }
    T*                  wghtTplGetAddrData(void) const { return  _wghtTpl.getAddrData(); }
    cSize3              wghtTplGetSize    (void) const { return  _wghtTpl.getSize();     }

    cData3<T>           wghtRefGet        (void)       { return  _wghtRef;               }
    cData3<T>*          wghtRefGetRef     (void)       { return &_wghtRef;               }
    T*                  wghtRefGetAddrData(void) const { return  _wghtRef.getAddrData(); }
    cSize3              wghtRefGetSize    (void) const { return  _wghtRef.getSize();     }

    void        maskCreateByCopy    (const cData3XCorr& object);                                                            // for gEMpicker
    void        maskCreateByRotation(const cData3XCorr& object, const cVector3<T>& angle, eInter inter = INTER_LINEAR);     // for gEMpicker
    void        maskCreateByOnes    (void);
    void        maskCreateByThresholding(T value);                                                                          // for gEMfitter
    void        maskCheckValue      (void);

    void        maskGetIndexLimit(size_t &indRowMin, size_t &indRowMax,
                                  size_t &indColMin, size_t &indColMax);
    void        maskShrinkToLimit(size_t  indRowMin, size_t  indRowMax,
                                  size_t  indColMin, size_t  indColMax);

    void        wghtTplCreate(unsigned int mode = 0);
    void        wghtRefCreateUsingWghtTpl(void);
};

} // namespace gem

#endif
