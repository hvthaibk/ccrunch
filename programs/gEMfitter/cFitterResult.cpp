/***********************************************************************
 *  File:       cFitterResult.cpp
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

#include "pft_vertices.hpp"

extern cData3<double>      powRef,     powTpl,     powMsk,     powWghtRef,     powWghtTpl;
#ifdef __GEM_USE_CUDA__
extern cuData3<double>     powRefCUDA, powTplCUDA, powMskCUDA, powWghtRefCUDA, powWghtTplCUDA;
#endif

void cFitter::computeResult(void)
{
    const std::string   funcName("void cFitter::computeResult(void)");

    float                   maxValAbs;
    size_t                  maxInd = 0, rotInd;
    cData3<float>           mapOcc, mapSub;
    cArray3<size_t>         mapWatershed;
    cFilterSpatial<float>   objFilter;

    cVector3<float>         angleFloat;
    cVector3<double>        angleDouble, transDouble;
    cVector3<double>        tplSize, tplOrigin, tplOffset, maxShift;
    cVector3<size_t>        maxPos;
    cResultFitter<double>   resultFitterVec;
    cString                 numStr;

    // original size
    mapRef.opFFTSizeRestore(sizeRefOriginal);

    tplSize   = mapTplHeader.getMapDim();
    tplOrigin = mapTplHeader.getMapOrigin();
    tplOffset = mapTplHeader.getMapOrigin() - mapRefHeader.getMapOrigin();

    // Gaussian smoothing correlation maps
    // clear correlation values that correspond to intensity < threshold
    // clear correlation values near boundary
    cVector3<size_t>    boundary((size_t) (tplSize.getMin()/3),
                                 (size_t) (tplSize.getMin()/3),
                                 (size_t) (tplSize.getMin()/3));

    mapOcc.memAllocZero(mapRef.getSize());
    mapOcc.opReplace(1, mapRef.getSize()-(size_t)2*boundary, boundary);

    objFilter.kernelSpatialGaussian(1, 3, 3);
    objFilter.filterSpatial(mapCor, XCORR_RES_SAME);
    mapCor *= mapOcc;

    // watershed transform for peak detection and suppression
    cTransDistance<float>().watershed(mapCor, mapWatershed, mapCor.getMax()/2, (float) tplSize.getMin()/5, WATER_DIR_FALL);

    std::cout << "ref size         : " << mapRefHeader.getMapDim()    << std::endl;
    std::cout << "ref origin       : " << mapRefHeader.getMapOrigin() << std::endl;
    std::cout << "tpl size         : " << mapTplHeader.getMapDim()    << std::endl;
    std::cout << "tpl origin       : " << mapTplHeader.getMapOrigin() << std::endl;
    std::cout << "tpl offset       : " << tplOffset                   << std::endl;
    std::cout << "tpl size extented: " << mapTpl.getSize()            << std::endl;
    std::cout << std::endl;

    // POWELL optimization
    if (bRefine) {
        switch (mCorr) {
            case CC_MODE_WCORR:
            case CC_MODE_CCORR:
                powWghtRef.opTypeCast(mapTpl.wghtTplGet());
                powWghtTpl.opTypeCast(mapTpl.wghtTplGet());
            case CC_MODE_XCORR:
            case CC_MODE_ECORR:
            case CC_MODE_NCORR:
            case CC_MODE_SCORR:
                powTpl.opTypeCast(mapTpl);
                powMsk.opTypeCast(mapTpl.maskGet());
                break;
            default:
                ERROR(funcName, "unsupported correlation mode");
        }

#ifdef __GEM_USE_CUDA__
        if (nGPU > 0) {
            switch (mCorr) {
                case CC_MODE_WCORR:
                case CC_MODE_CCORR:
                    powWghtRefCUDA = powWghtRef;
                    powWghtTplCUDA = powWghtTpl;
                case CC_MODE_XCORR:
                case CC_MODE_ECORR:
                case CC_MODE_NCORR:
                case CC_MODE_SCORR:
                    powTplCUDA = powTpl;
                    powMskCUDA = powMsk;
                    break;
                default:
                    ERROR(funcName, "unsupported correlation mode");
            }
        }
#endif
    }

    for (size_t iCom = 0; iCom < nComponent; iCom++) {
        // get the index of the current maxima
        maxInd    = mapCor.getMaxIndx();
        maxValAbs = mapCor[maxInd];
        rotInd    = (size_t) mapInd[maxInd] - 1;

        // get the rotation angle
        pft_angle(rotInd, nNode, nRotZ, angleFloat);
        angleDouble = angleFloat;

        // get the coordinates of the current maxima
        maxPos.ind2sub(mapCor.getSize(), maxInd);

        // refinement
        transDouble = 0;
        if (bRefine) {
            std::cout << "Off-latice refinement for iCom = " << iCom+1 << std::endl;
            std::cout << "     initial rotation   : " << angleDouble << std::endl;
            std::cout << "             translation: " << transDouble << std::endl;
            std::cout << "             maxValAbs  : " << maxValAbs   << std::endl;

            // get the sub-volume
            cVector3<ptrdiff_t>     offset((ptrdiff_t) maxPos[0] - (ptrdiff_t) (mapTpl.getNrow()-1)/2,
                                           (ptrdiff_t) maxPos[1] - (ptrdiff_t) (mapTpl.getNcol()-1)/2,
                                           (ptrdiff_t) maxPos[2] - (ptrdiff_t) (mapTpl.getNsec()-1)/2);
            //mapSub.opCrop(mapRef, mapTpl.getSize(), maxPos-offset);
            mapSub.opCopy(mapRef, mapTpl.getSize(), offset);

            powRef.opTypeCast(mapSub);
#ifdef __GEM_USE_CUDA__
            if (nGPU > 0) { powRefCUDA = powRef; }
#endif

            // powell optimization
            //maxValAbs = (float) -refineNR3(transDouble, angleDouble);
            maxValAbs = (float) -refineGSL(transDouble, angleDouble);

            std::cout << "     final   rotation   : " << angleDouble << std::endl;
            std::cout << "             translation: " << transDouble << std::endl;
            std::cout << "             maxValAbs  : " << maxValAbs   << std::endl;
            std::cout << std::endl;
        }

        // get the translation distance
        maxShift  = maxPos;
        maxShift -= cVector3<double>((double) halfUpper((size_t) tplSize[0]),
                                     (double) halfUpper((size_t) tplSize[1]),
                                     (double) halfUpper((size_t) tplSize[2])) + tplOffset - 1;
        maxShift += transDouble;
        maxShift *= mapRefHeader.getVoxelSpacing();

        // apply transformation on template
        // NOTE: remember to reverse the order of rotation
        cVector3<double>        distToOrigin = (tplSize/2+tplOrigin-0.5)*mapRefHeader.getVoxelSpacing();
        cPDB                    pdb;

        numStr = num2str(iCom+1);
        numStr.insertBegin('0',3);

        pdb.read(fileNameTpl);
        pdb.opTransRotTrans(-distToOrigin,
                            cVector3<double>(angleDouble[2],
                                             angleDouble[1],
                                             angleDouble[0]),
                            distToOrigin+maxShift);
        pdb.write(dirNameFit + "gEMfitter__best__" + numStr + ".pdb");

        /*// Use cPDB as above to get rid of ESBTL
        cESBTL().pdbTransRotTrans(fileNameTpl,
                        dirNameFit + "gEMfitter__best__" + numStr + ".pdb",
                        -distToOrigin,
                        cVector3<double>(angleDouble[2],
                                         angleDouble[1],
                                         angleDouble[0]),
                        distToOrigin+maxShift);*/

        // mask out the correlation map
        for (size_t ind = 0; ind < mapWatershed.getNelement(); ind++) {
            if (mapWatershed[ind] == (iCom+1)) {
                mapCor[ind] = 0;
            }
        }

        // store coordinates
        sFitterResult<double>   resultFitter(numStr, maxShift, angleDouble, (double) maxValAbs);
        resultFitterVec.push_back(resultFitter);
    }

    // save coordinates to file
    resultFitterVec.write(dirNameFit + "gEMfitter_fitting_result.txt");

    if (bRefine) {
        powRef.memFree();
        powTpl.memFree();
        powMsk.memFree();
        powWghtRef.memFree();
        powWghtTpl.memFree();

#ifdef __GEM_USE_CUDA__
        if (nGPU > 0) {
            powRefCUDA.memFree();
            powTplCUDA.memFree();
            powMskCUDA.memFree();
            powWghtRefCUDA.memFree();
            powWghtTplCUDA.memFree();
        }
#endif
    }
}

void cFitter::computeRMSD(void)
{
    const std::string   funcName("void cFitter::computeRMSD(void)");

    double                  rmsd;

    if (chainName == "all") {
        require(nComponent == 1,
                funcName + ": #components != 1");

        cPDB    pdb1, pdb2;

        pdb1.read(fileNameMdl);
        pdb2.read(dirNameFit + "gEMfitter__best__001.pdb");

        rmsd = pdb1.opComputeRMSD(pdb2, true);
    }
    else {
        require(chainName.length() == nComponent,
                funcName + ": #chains != #components");

        std::vector<std::string>    fileNameGT (nComponent);
        std::vector<std::string>    fileNameEst(nComponent);
        std::string                 chainID = "";
        cString                     numStr;

        for (size_t iCom = 0; iCom < nComponent; iCom++) {
            chainID = chainName[iCom];

            cESBTL().pdbExtractChain(fileNameMdl,
                                     dirNameFit + chainName[iCom] + ".pdb",
                                     chainID);

            numStr = num2str(iCom+1);
            numStr.insertBegin('0',3);

            fileNameGT[iCom]  = dirNameFit + chainName[iCom] + ".pdb";
            fileNameEst[iCom] = dirNameFit + "gEMfitter__best__" + numStr + ".pdb";

            //PRINT(fileNameGT[iCom]);
            //PRINT(fileNameEst[iCom]);
        }

        std::vector< cSetPoint3<double> >   atomListGT(nComponent);
        std::vector< cSetPoint3<double> >   atomListEst(nComponent);

        for (size_t iCom = 0; iCom < nComponent; iCom++) {
            cESBTL().pdbReadAtom(fileNameGT[iCom],  atomListGT[iCom],  "CA");
            cESBTL().pdbReadAtom(fileNameEst[iCom], atomListEst[iCom], "CA");

            //PRINT(atomListGT[iCom].size());
            //PRINT(atomListEst[iCom].size());
        }

        cArray3<double>         sqrDist(cVector3<size_t>(nComponent,nComponent,1));
        double                  sqrDistMax;

        for (size_t iCom = 0; iCom < nComponent; iCom++) {
            for (size_t jCom = 0; jCom < nComponent; jCom++) {
                sqrDist[iCom*nComponent+jCom] = atomListGT[iCom].computeSD(atomListEst[jCom]);
            }
        }
        sqrDistMax = sqrDist.getMax();
        //sqrDist.printData("sqrDist matrix", 15, 4);

        size_t                  minIndx;
        std::vector<size_t>     setMinRow(nComponent), setMinCol(nComponent);
        double                  sqrDistAll = 0;
        size_t                  numAtomAll = nComponent * atomListGT[0].size();

        for (size_t iCom = 0; iCom < nComponent; iCom++) {
            minIndx = sqrDist.getMinIndx();

            sqrDistAll += sqrDist[minIndx];

            setMinRow[iCom] = minIndx / sqrDist.getNcol();
            setMinCol[iCom] = minIndx - setMinRow[iCom]*sqrDist.getNcol();

            for (size_t jCom = 0; jCom < nComponent; jCom++) {
                sqrDist[setMinRow[iCom]*nComponent + jCom] = sqrDistMax;
                sqrDist[jCom*nComponent + setMinCol[iCom]] = sqrDistMax;
            }
        }

        rmsd = std::sqrt(sqrDistAll / (double) numAtomAll);
    }

    std::cout << "     rmsd = " << rmsd << "\n\n";

    std::ofstream       outfile;

    outfile.open("rmsd.txt", std::ofstream::app |
                             std::ofstream::out);
    outfile << fileNameResCorr << std::fixed << std::setw(20) << rmsd << std::endl;

    outfile.close();
}
