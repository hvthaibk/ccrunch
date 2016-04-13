/***********************************************************************
 *  File:       cFitterData.cpp
 *
 *  Purpose:    Implementation of cFitter class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cFitter.hpp"

void cFitter::dataReadInput(void)
{
    const std::string   funcName("void cFitter::dataReadInput(void)");

    // REFERENCE
    mpiMasterPrint(mpiID, "Reading input target map...\n");
    cCCP4<float>().read(fileNameRef, mapRef, &mapRefHeader);
    require(mapRef.getDimension() == 3, "target data is not 3D");
    if (mpiID == 0) mapRef.printSize();

    if (bNyquist) {

        if (mapRefHeader.getVoxelSpacing() < 0.2*resoRef) {

            float               voxelSpacingNyquist = 0.25f * (float) resoRef;
            float               factor = (float) mapRefHeader.getVoxelSpacing() / voxelSpacingNyquist;
            cVector3<float>     mapRefOrigin;

            mpiMasterPrint(mpiID, "Nyquist-rate resampling: " + num2str(mapRefHeader.getVoxelSpacing(),4) +
                                         " -> " + num2str(voxelSpacingNyquist,4) + "\n");

            mapRefOrigin = mapRefHeader.getMapOrigin();
            cTransSimilarity<float>().scale(mapRef, mapRefOrigin, factor, INTER_LINEAR);

            mapRefHeader.updateAfterScaling(factor);

            if (mpiID == 0) mapRef.printSize();
        }
    }

    // TEMPLATE
    mpiMasterPrint(mpiID, "Reading search pdb and converting it to map...\n");
    pdbTpl.read(fileNameTpl);
    pdbTpl.opPDB2Map(resoRef, mapRefHeader.getVoxelSpacing(), PDB2VOl_USE_HETATM_NO, mapTpl, &mapTplHeader);
    if (mpiID == 0) mapTpl.printSize();

    // Map extension
    mpiMasterPrint(mpiID, "Extending search map...\n");
    cTransSimilarity<float>().extendRotate(mapTpl, EXP_ROT_VALID);
    if (mpiID == 0) mapTpl.printSize();

    mpiMasterPrint(mpiID, "Extending target map...\n");
    mapRef.opFFTSizeChange(sizeRefOriginal, mapTpl.getSize());
    //mapRef.opPad(cVector3<size_t>( 64, 64, 64), cVector3<size_t>(0,0,0));
    //mapRef.opPad(cVector3<size_t>(128,128,128), cVector3<size_t>(0,0,0));
    //mapRef.opPad(cVector3<size_t>(256,256,256), cVector3<size_t>(0,0,0));
    if (mpiID == 0) mapRef.printSize();

#ifdef __GEM_LIMIT_FUNC__
    require(mapRef.getSize().getMax() <= 256, funcName + "target data size is too-large");
    require(mCorr <= 2, funcName + "unsupported correlation mode");
#endif

    // MASK and WEIGHT
    if (fileNameMsk != "none") {
        mpiMasterPrint(mpiID, "Reading mask from file...\n");
        cCCP4<float>().read(fileNameMsk, *mapTpl.maskGetRef(), &mapMskHeader);
        require(mapTplHeader.getMapDim()    == mapMskHeader.getMapDim(),    "tplSize != mskSize  ");
        require(mapTplHeader.getMapOrigin() == mapMskHeader.getMapOrigin(), "tplOrigin != mskOrigin");
    }
    else {
        mpiMasterPrint(mpiID, "Creating mask from search data...\n");
        if (threshMsk > 0) {
            mapTpl.maskCreateByThresholding((float) threshMsk);
        }
        else if (threshMsk == 0) {
            mapTpl.maskCreateByThresholding(mapTpl.getMin()-1.0f);
        }
        else {
            mapTpl.maskCreateByThresholding(mapTpl.getMax()/std::abs((float) threshMsk));
        }
    }

    mpiMasterPrint(mpiID, "\n");

    switch (mCorr) {
        case CC_MODE_XCORR:
            mapTpl *= mapTpl.maskGet();     // because mask is not used in XCORR computation
            break;
        case CC_MODE_ECORR:
            break;
        case CC_MODE_NCORR:
            break;
        case CC_MODE_WCORR:                 // prepare the weights
            mpiMasterPrint(mpiID, "Creating weight from mWCC...\n");
            mapTpl.wghtTplCreate(mWCC);
            mapTpl.wghtRefCreateUsingWghtTpl();
            mpiMasterPrint(mpiID, "\n");
            break;
        case CC_MODE_SCORR:                 // prepare the surfaces
            mpiMasterPrint(mpiID, "Creating surface maps...\n");
            cTransMorpho<float>().dilate(mapTpl.maskGet(), mapTpl, 5);
            mapTpl -= mapTpl.maskGet();
            mapTpl /= mapTpl.getSum();

            mapRef.maskCreateByThresholding(mapRef.getMax()/std::abs((float) threshMsk));
            cTransMorpho<float>().dilate(mapRef.maskGet(), mapRef, 5);
            mapRef -= mapRef.maskGet();
            mpiMasterPrint(mpiID, "\n");
            break;
        case CC_MODE_CCORR:
            // prepare the weights
            mpiMasterPrint(mpiID, "Creating weight from mWCC...\n");
            mapTpl.wghtTplCreate(mWCC);
            mapTpl.wghtRefCreateUsingWghtTpl();
            mpiMasterPrint(mpiID, "\n");

            // prepare the surfaces
            mpiMasterPrint(mpiID, "Creating surface maps...\n");
            cTransMorpho<float>().dilate(mapTpl.maskGet(), mapTplSurface, 3);
            mapTplSurface -= mapTpl.maskGet();
            mapTplSurface /= mapTplSurface.getSum();

            //mapRef.maskCreateByThresholding(mapRef.getMax()/std::abs((float) threshMsk));
            mapRef.maskCreateByThresholding(mapRef.getMax()/std::abs((float) 4));
            cTransMorpho<float>().dilate(mapRef.maskGet(), mapRefSurface, 3);
            mapRefSurface -= mapRef.maskGet();
            mpiMasterPrint(mpiID, "\n");

            // prepare the penalty
            mpiMasterPrint(mpiID, "Creating penalty maps...\n");
            mapTplInner.memReAllocZero(mapTpl.getSize());
            mapTplInner += mapTpl.maskGet();
            mapTplInner /= mapTplInner.getSum();

            mapRefOuter.memReAllocZero(mapRef.getSize());
            mapRefOuter.memSetVal(1);
            mapRefOuter -= mapRef.maskGet();
            mpiMasterPrint(mpiID, "\n");
            break;
        default:
            ERROR(funcName, "unsupported correlation mode");
    }

#ifdef NDEBUG
    mpiMasterPrint(mpiID, "Writing debugging data to disk...\n");
    dataWriteTmp();
    mpiMasterPrint(mpiID, "\n");
#endif

    // PREPROCESSING
    if (bLaplacian) {
        mpiMasterPrint(mpiID, "Laplacian filtering target and search maps...\n");

        cFilterSpatial<float>   objFilter;

        objFilter.kernelSpatialLaplacian(3.0/16, 0, 3, 3);
        objFilter.filterSpatial(mapRef, XCORR_RES_SAME);
        objFilter.filterSpatial(mapTpl, XCORR_RES_SAME);
    }

    switch (mCorr) {
        case CC_MODE_XCORR:
            break;
        case CC_MODE_ECORR:
        case CC_MODE_NCORR:
        case CC_MODE_WCORR:
        case CC_MODE_CCORR:
            mpiMasterPrint(mpiID, "Adding noise to target map...\n");
            cRandom<float>().randZero(mapRef, 0.0f, mapRef.getMax()/1000.0f, mapRef.getMax()/100.0f);
            mpiMasterPrint(mpiID, "\n");
            break;
        case CC_MODE_SCORR:
            break;
        default:
            ERROR(funcName, "unsupported correlation mode");
    }
}

void cFitter::dataReadCorr(void)
{
    if (mpiID == 0) {
        cCCP4<float>().read(fileNameResCorr, mapCor);
        cCCP4<float>().read(fileNameResIndx, mapInd);
    }
}

void cFitter::dataWriteCorr(void)
{
    if (mpiID == 0) {
        cCCP4<float>().write(mapCor, fileNameResCorr, &mapRefHeader);
        cCCP4<float>().write(mapInd, fileNameResIndx, &mapRefHeader);
    }
}

void cFitter::dataWriteTmp(void)
{
    const std::string   funcName("void cFitter::dataWriteTmp(void)");

    bool            needDebugData = true;

    cMapHeader      mapRefHeaderExt = mapRefHeader;
    cMapHeader      mapTplHeaderExt = mapTplHeader;

    mapRefHeaderExt.updateAfterPadding  (mapRef.getSize(), cVector3<size_t>(0,0,0));
    mapTplHeaderExt.updateAfterExtension(mapTpl.getSize());

    if (mpiID == 0) {
        switch (mCorr) {
            case CC_MODE_XCORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapTpl.map").c_str());
                }
                break;
            case CC_MODE_ECORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);
                    cCCP4<float>().write(*mapTpl.maskGetRef(), dirNameFit + "mapMsk.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapTpl.map").c_str());
                    remove((dirNameFit + "mapMsk.map").c_str());
                }
                break;
            case CC_MODE_NCORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);
                    cCCP4<float>().write(*mapTpl.maskGetRef(), dirNameFit + "mapMsk.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapTpl.map").c_str());
                    remove((dirNameFit + "mapMsk.map").c_str());
                }
                break;
            case CC_MODE_WCORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);
                    cCCP4<float>().write(*mapTpl.maskGetRef(), dirNameFit + "mapMsk.map", &mapTplHeaderExt);
                    cCCP4<float>().write(*mapTpl.wghtTplGetRef(), dirNameFit + "mapWghtTpl.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapTpl.map").c_str());
                    remove((dirNameFit + "mapMsk.map").c_str());
                    remove((dirNameFit + "mapWghtTpl.map").c_str());
                }
                break;
            case CC_MODE_SCORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapRef, dirNameFit + "mapRef.map", &mapRefHeaderExt);
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapRef.map").c_str());
                    remove((dirNameFit + "mapTpl.map").c_str());
                }
                break;
            case CC_MODE_CCORR:
                if (needDebugData) {
                    cCCP4<float>().write(mapRef, dirNameFit + "mapRef.map", &mapRefHeaderExt);
                    cCCP4<float>().write(mapTpl, dirNameFit + "mapTpl.map", &mapTplHeaderExt);

                    cCCP4<float>().write(mapRefSurface, dirNameFit + "mapRefSurface.map", &mapRefHeaderExt);
                    cCCP4<float>().write(mapTplSurface, dirNameFit + "mapTplSurface.map", &mapTplHeaderExt);

                    cCCP4<float>().write(mapRefOuter, dirNameFit + "mapRefOuter.map", &mapRefHeaderExt);
                    cCCP4<float>().write(mapTplInner, dirNameFit + "mapTplInner.map", &mapTplHeaderExt);
                }
                else {
                    remove((dirNameFit + "mapRef.map").c_str());
                    remove((dirNameFit + "mapTpl.map").c_str());

                    remove((dirNameFit + "mapRefSurface.map").c_str());
                    remove((dirNameFit + "mapTplSurface.map").c_str());

                    remove((dirNameFit + "mapRefOuter.map").c_str());
                    remove((dirNameFit + "mapTplInner.map").c_str());
                }
                break;
            default:
                ERROR(funcName, "unsupported correlation mode");
        }
    }
}
