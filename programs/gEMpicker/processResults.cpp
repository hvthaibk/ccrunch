/***********************************************************************
 *  File:       processResults.cpp
 *
 *  Purpose:    Implementation of post-processing functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "gEMpicker.hpp"

void pickParticle(sPickerParams& params)
{
    assert(params.threshold >= 0 && params.threshold <= 1);

    std::string                 dirName;
    std::string                 fileNameTmp, fileNameRef, fileNameTpl;
    float                       peakVal;
    size_t                      peakInd = 0, peakIndRow = 0, peakIndCol = 0;
    cData3XCorr<float>          imRef, imTpl, imMskRef, imMskTpl, imOut, imOutWithTpl, imCor, imInd;
    size_t                      tplInd;
    size_t                      indMskRowMin, indMskRowMax, indMskColMin, indMskColMax;
    int                         indRefRowMin, indRefRowMax, indRefColMin, indRefColMax;
    bool                        bChangeMsk, bOutOfBorder;
    unsigned int                nPicked;
    cString                     numStr;
    size_t                      iRef, iTpl, iRot;
    cVector3<float>             angle = 0.0f;
    cResultPicker<float>        resultPickerVec;
    size_t                      nRemovedParticle;

    //for (iRef = 0; iRef <= 0; iRef++) {
    for (iRef = 0; iRef < params.fileListRef.size(); iRef++) {
        std::cout << "Processing micrograph " << iRef+1 << "/"
                  << params.fileListRef.size() << ": "
                  << params.fileListRef[iRef]
                  << std::endl;

        fileNameRef = cFileSystem(params.fileListRef[iRef]).getFileRoot();

        nPicked = 0;
        resultPickerVec.clear();

        // read reference image (original micrograph)
        cMRC<float>().read(params.fileListRef[iRef], imRef);

        // read xcorr and index images (results from previous xcorr computation)
        fileNameTmp = params.dirResStr + "xcorr/" + fileNameRef + "__corr.mrc";
        cMRC<float>().read(fileNameTmp, imCor);

        fileNameTmp = params.dirResStr + "xcorr/" + fileNameRef + "__indx.mrc";
        cMRC<float>().read(fileNameTmp, imInd);

        imCor.opClearBorderMargin(5);   // this is critical to avoid error
                                        // with the peaks at the boundary

        // reset the tpl counting array to zero for each reference image
        std::vector<size_t>     numTplSelected(params.fileListTpl.size(), 0);

        // mask to correlation maps
        fileNameTmp = params.dirMskRefStr + fileNameRef + ".tif";
        if (cFileSystem(fileNameTmp).isExist()) {
            cTIFF<float>().read(fileNameTmp,imMskRef);
            imMskRef /= 255;
            imCor *= imMskRef;
        }

        while (1) {
            // get the index of the current maxima
            if (params.contrast) {      // NEGATIVELY STAINED
                peakInd = imCor.getMaxIndx();
            }
            else {                      // CRYOGENIC
                peakInd = imCor.getMinIndx();
            }
            require(peakInd > 1, "something wrong with peakInd in gEMpicker");
            // I cannot explain why we need to use (peakInd-1)
            // Maybe due to the dimension definition in Fourier transform (ncol/2+1) ?
            // This part needs to be further investigated so that we can avoid calling
            //     imCor.opClearBorderMargin(5);

            peakIndRow = (peakInd-1) / imCor.getNcol();
            peakIndCol = (peakInd-1) % imCor.getNcol();

            peakVal    = imCor[peakInd];
            tplInd     = (size_t) imInd[peakInd] - 1;

            iTpl = tplInd / params.numRot2D;
            iRot = tplInd % params.numRot2D;
            angle[0] = deg2rad((float) iRot * params.angle2D);

            if (std::abs(peakVal) < std::abs(params.threshold)) break;

            // read the corresponding template
            cMRC<float>().read(params.fileListTpl[iTpl], imTpl);

            // read the corresponding mask and compute its extends
            // need to re-read even for a single mask since it may be
            // cropped in the previous iteration
            if (params.bMskSingle) { cTIFF<float>().read(params.fileListMsk[0],    imMskTpl); }
            else                   { cTIFF<float>().read(params.fileListMsk[iTpl], imMskTpl); }

            cTransSimilarity<float>().rotate(imTpl,    angle);
            cTransSimilarity<float>().rotate(imMskTpl, angle);

            imMskTpl.maskGetIndexLimit(indMskRowMin, indMskRowMax,
                                       indMskColMin, indMskColMax);

            // compute the coordinates of mask's extends in the reference
            indRefRowMin = (int) peakIndRow - (int) (imMskTpl.getNrow()/2 - indMskRowMin);
            indRefRowMax = (int) peakIndRow + (int) (indMskRowMax - imMskTpl.getNrow()/2);
            indRefColMin = (int) peakIndCol - (int) (imMskTpl.getNcol()/2 - indMskColMin);
            indRefColMax = (int) peakIndCol + (int) (indMskColMax - imMskTpl.getNcol()/2);

            // adjust the mask' extends to avoid out-of-border
            bChangeMsk = false;
            if (indRefRowMin < 0) {
                //indMskRowMin = imMskTpl.getNrow()/2 - peakIndRow - 1;
                indMskRowMin = imMskTpl.getNrow()/2 - peakIndRow;
                indRefRowMin = 0;
                bChangeMsk = true;
            }
            if (indRefRowMax >= (int) imRef.getNrow()) {
                //indMskRowMax = imRef.getNrow() - peakIndRow + imMskTpl.getNrow()/2 - 2;
                indMskRowMax = imRef.getNrow() - peakIndRow + imMskTpl.getNrow()/2 - 1;
                indRefRowMax = (int) imRef.getNrow() - 1;
                bChangeMsk = true;
            }
            if (indRefColMin < 0) {
                //indMskColMin = imMskTpl.getNcol()/2 - peakIndCol - 1;
                indMskColMin = imMskTpl.getNcol()/2 - peakIndCol;
                indRefColMin = 0;
                bChangeMsk = true;
            }
            if (indRefColMax >= (int) imRef.getNcol()) {
                //indMskColMax = imRef.getNcol() - peakIndCol + imMskTpl.getNcol()/2 - 2;
                indMskColMax = imRef.getNcol() - peakIndCol + imMskTpl.getNcol()/2 - 1;
                indRefColMax = (int) imRef.getNcol() - 1;
                bChangeMsk = true;
            }

            if (bChangeMsk) {
                imMskTpl.maskShrinkToLimit(indMskRowMin, indMskRowMax,
                                           indMskColMin, indMskColMax);
            }

            // calculate boxSize if not provided
            if (params.boxSize == 0) {
                params.boxSize = imTpl.getNrow();
            }

            // mask out the region around the current maxima in imCor
            /*imCor.opReplaceUsingMask(imMskTpl, 0, cVector3<ptrdiff_t>(peakIndRow + 1 - imMskTpl.getNrow()/2,
                                                                      peakIndCol + 1 - imMskTpl.getNcol()/2,
                                                                      0));*/
            imCor.opReplaceUsingMask(imMskTpl, 0, cVector3<ptrdiff_t>(peakIndRow - imMskTpl.getNrow()/2,
                                                                      peakIndCol - imMskTpl.getNcol()/2,
                                                                      0));

            // check out-of-border
            bOutOfBorder = false;
            if ((peakIndRow + 1 < params.boxSize/2 + params.boxBorder) ||
                (peakIndCol + 1 < params.boxSize/2 + params.boxBorder) ||
                (peakIndRow + params.boxSize/2 + params.boxBorder >= imRef.getNrow()) ||
                (peakIndCol + params.boxSize/2 + params.boxBorder >= imRef.getNcol())) {
                bOutOfBorder = true;
            }

            if (!bChangeMsk && !bOutOfBorder) {
                nPicked++;

                fileNameTpl = cFileSystem(params.fileListTpl[iTpl]).getFileRoot();

                // save the selected particle
                require(nPicked <= 9999, "gEMpicker have already picked 9999 "
                                         "particles from a single microgaph, "
                                         "you should increase your threshold value");

                numStr = num2str(nPicked);
                numStr.insertBegin('0',4);

                // increase the corresponding element in the tpl counting array
                numTplSelected[iTpl]++;

                // store coordinates
                sPickerResult<float>    resultPicker(numStr, peakIndCol, imCor.getNrow() - 1 - peakIndRow ,
                                                     params.boxSize, iTpl+1,
                                                     cFileSystem(params.fileListTpl[iTpl]).getFileRoot(),
                                                     cFileSystem(params.bMskSingle ? params.fileListMsk[0] : params.fileListMsk[iTpl]).getFileRoot(),
                                                     (float) iRot*params.angle2D, peakVal);
                resultPickerVec.push_back(resultPicker);

                if (nPicked >= params.nPickMax) break;
            }
        }

        // refine picking results
        if (params.boxDist == 0) {
            params.boxDist = imTpl.getNrow();
        }

        nRemovedParticle = 0;
        if (resultPickerVec.size() > 0) {
            nRemovedParticle += resultPickerVec.refineUsingDistance(params.boxDist);
        }
        if (resultPickerVec.size() > 0) {
            nRemovedParticle += resultPickerVec.refineUsingHighThreshold(params.thresholdHigh);
        }

        // save particles
        if (params.mode == 2 && resultPickerVec.size() > 0) {

            // create the corresponding subdirectory to save image if necessary
            dirName = params.dirResStr + "pik_ext/" + fileNameRef;
            if (cFileSystem(dirName).isExist()) {
                cFileSystem(dirName).dirDelete();
            }
            cFileSystem(dirName).dirCreate();

            imOut.memReAlloc(imMskTpl.getSize());
            imOutWithTpl.memReAllocZero(cVector3<size_t>(std::max(imMskTpl.getNrow(),params.boxSize),
                                                         imMskTpl.getNcol()+params.boxSize,
                                                         1));

            sPickerResult<float>    result;

            for (size_t iPik = 0; iPik < resultPickerVec.size(); iPik++) {

                result = resultPickerVec.at(iPik);

                peakIndRow = imCor.getNrow() - 1 - result._y;
                peakIndCol = result._x;

                // extract the particle from imRef and store it in imOut
                imOut.opCrop(imRef, cVector3<size_t>(params.boxSize,
                                                     params.boxSize,
                                                     1),
                                    cVector3<size_t>(peakIndRow + 1 - params.boxSize/2,
                                                     peakIndCol + 1 - params.boxSize/2,
                                                     0));

                fileNameTmp = dirName + "/" + result._numStr
                                    + "__" + result._tplName
                                    + "__" + num2str(result._tplIdx)
                                    + "__" + num2str(result._angle2D);

                cMRC<float>().write(imOut, fileNameTmp + ".mrc");

                // put correlation value onto the combined image = picked + tpl
                // the two images imOut and imTpl need to be normalized before combining
                // -> they have to be reloaded in the next iteration

                cMRC<float>().read(params.fileListTpl[result._tplIdx-1], imTpl);
                cTransSimilarity<float>().rotate(imTpl, deg2rad((float) result._angle2D));

                imOut.normalize8bit();
                imTpl.normalize8bit();
                imOutWithTpl.opReplace(imOut, cVector3<size_t>(0,0,0));
                imOutWithTpl.opReplace(imTpl, cVector3<size_t>(0,imOut.getNcol(),0));

#ifdef __GEM_USE_OPENCV__
                cOpenCV<float>().putText(imOutWithTpl, fileNameTmp + ".tif",
                                         num2str(result._corrVal, 4),
                                         COLOR_BLUE, 1.0f, NORM8BIT_FALSE);
#else
                cTIFF<float>().write(imOutWithTpl, fileNameTmp + ".tif", NORM8BIT_FALSE);
#endif
            }
        }

        // save coordinates to file
        fileNameTmp = params.dirResStr + "pik_coord/" + fileNameRef + ".txt";
        resultPickerVec.write(fileNameTmp);
        fileNameTmp = params.dirResStr + "pik_box/"   + fileNameRef + ".box";
        resultPickerVec.writeBoxEMAN(fileNameTmp, imRef.getNrow());

        // message
        if (nPicked >= params.nPickMax) {
            std::cout << "     nPickMax reached,  nPicked  = " << resultPickerVec.size() << ", peakVal = "  << peakVal<< "\n";
            std::cout << "                        nRemoved = " << nRemovedParticle << "\n\n";
        }
        else {
            std::cout << "     threshold reached, nPicked  = " << resultPickerVec.size() << "\n";
            std::cout << "                        nRemoved = " << nRemovedParticle << "\n\n";
        }
    }

    mpiMasterPrint(params.mpiID, "\n");
}
