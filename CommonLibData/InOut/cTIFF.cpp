/***********************************************************************
 *  File:       cTIFF.cpp
 *
 *  Purpose:    Implementation of a TIFF class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cTIFF.hpp"

#include "macro.hpp"

#include "tiffio.h"

namespace gem {

template <typename T>
void cTIFF<T>::checkExtension(const std::string& fileName) const
{
    const std::string   funcName("void cTIFF::checkExtension(const std::string& fileName) const");

    if (fileName == "" || (fileName.find(".tif" ) > fileName.length() &&
                           fileName.find(".tiff") > fileName.length() &&
                           fileName.find(".TIF" ) > fileName.length() &&
                           fileName.find(".TIFF") > fileName.length())) {
        ERROR(funcName, "unsupported expansion (accept *.tif, *.tiff, *.TIF and *.TIFF)");
    }
}

template <typename T>
void cTIFF<T>::getInfo(const std::string& fileName) const
{
    const std::string   funcName("void cTIFF::getInfo(const std::string& fileName) const");

    checkExtension(fileName);

    TIFF*   fileHandle = TIFFOpen(fileName.c_str(), "r");

    require(fileHandle, funcName + ": could not open TIFF file " + fileName);

    uint32  imageWidth, imageLength;
    uint16  samplePerPixel, bitPerSample, sampleFormat;
    uint32  tileWidth, tileLength;

    require(TIFFGetField         (fileHandle, TIFFTAG_IMAGEWIDTH,      &imageWidth)     == TIFF_SUCCESS,
            "cannot get field TIFFTAG_IMAGEWIDTH in " + funcName);
    require(TIFFGetField         (fileHandle, TIFFTAG_IMAGELENGTH,     &imageLength)    == TIFF_SUCCESS,
            "cannot get field TIFFTAG_IMAGELENGTH in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_SAMPLESPERPIXEL, &samplePerPixel) == TIFF_SUCCESS,
            "cannot get field TIFFTAG_SAMPLESPERPIXEL in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_BITSPERSAMPLE,   &bitPerSample)   == TIFF_SUCCESS,
            "cannot get field TIFFTAG_BITSPERSAMPLE in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_SAMPLEFORMAT,    &sampleFormat)   == TIFF_SUCCESS,
            "cannot get field TIFFTAG_SAMPLEFORMAT in " + funcName);

    std::cout << "\nTIFF HEADER INFO: "   << fileName       << std::endl;
    std::cout << "     imageWidth     = " << imageWidth     << std::endl;
    std::cout << "     imageLength    = " << imageLength    << std::endl;
    std::cout << "     samplePerPixel = " << samplePerPixel << std::endl;
    std::cout << "     bitPerSample   = " << bitPerSample   << std::endl;
    switch (sampleFormat) {
        case SAMPLEFORMAT_UINT:
            std::cout << "     sampleFormat   = UINT"  << std::endl;
            break;
        case SAMPLEFORMAT_INT:
            std::cout << "     sampleFormat   = INT"   << std::endl;
            break;
        case SAMPLEFORMAT_IEEEFP:
            std::cout << "     sampleFormat   = FLOAT" << std::endl;
            break;
        default:
            std::string     errorMsg("unsupported sample format = ");
            errorMsg.append(num2str(sampleFormat));
            ERROR(funcName, errorMsg.c_str());
    }

    if (TIFFGetField(fileHandle, TIFFTAG_TILEWIDTH,  &tileWidth)  == TIFF_SUCCESS)
        std::cout << "     tileWidth      = " << tileWidth      << std::endl;
    if (TIFFGetField(fileHandle, TIFFTAG_TILELENGTH, &tileLength) == TIFF_SUCCESS)
        std::cout << "     tileLength     = " << tileLength     << std::endl;

    std::cout << std::endl;

    TIFFClose(fileHandle);
}

template <typename T>
void cTIFF<T>::read(const std::string& fileName, cArray3<T>& object) const
{
    const std::string   funcName("void cTIFF::read(const std::string& fileName, "
                                    "cArray3<T>& object) const");

    checkExtension(fileName);

    TIFF*   fileHandle = TIFFOpen(fileName.c_str(), "r");

    require(fileHandle, funcName + ": could not open TIFF file " + fileName);

    uint32  nRow, nCol;
    uint16  samplePerPixel, bitPerSample, sampleFormat;

    require(TIFFGetField(fileHandle, TIFFTAG_IMAGELENGTH, &nRow) == TIFF_SUCCESS,
            "cannot get field TIFFTAG_IMAGELENGTH in " + funcName);
    require(TIFFGetField(fileHandle, TIFFTAG_IMAGEWIDTH,  &nCol) == TIFF_SUCCESS,
            "cannot get field TIFFTAG_IMAGEWIDTH in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_SAMPLESPERPIXEL, &samplePerPixel) == TIFF_SUCCESS,
            "cannot get field TIFFTAG_SAMPLESPERPIXEL in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_BITSPERSAMPLE,   &bitPerSample)   == TIFF_SUCCESS,
            "cannot get field TIFFTAG_BITSPERSAMPLE in " + funcName);
    require(TIFFGetFieldDefaulted(fileHandle, TIFFTAG_SAMPLEFORMAT,    &sampleFormat)   == TIFF_SUCCESS,
            "cannot get field TIFFTAG_SAMPLEFORMAT in " + funcName);

    require(samplePerPixel == 1, funcName + " supports single channel data only");

    // allocate new memory
    object.memReAlloc(cVector3<size_t>(nRow, nCol, 1));

    unsigned char*      bufferData   = NULL;
    size_t              bufferLength = nRow*nCol*bitPerSample/8;

    array_new(bufferData, bufferLength);

    for (uint32_t iRow = 0; iRow < nRow; iRow++) {
        require(TIFFReadScanline(fileHandle, bufferData+iRow*nCol*bitPerSample/8, iRow) == TIFF_SUCCESS,
                "TIFFReadScanline() failed in " + funcName);
    }

    switch (sampleFormat) {
        case SAMPLEFORMAT_UINT:
            switch (bitPerSample) {
                case 8:     array_typecast((uint8_t*)  bufferData, object.getNelement(), object.getAddrData());     break;
                case 16:    array_typecast((uint16_t*) bufferData, object.getNelement(), object.getAddrData());     break;
                case 32:    array_typecast((uint32_t*) bufferData, object.getNelement(), object.getAddrData());     break;
                case 64:    array_typecast((uint64_t*) bufferData, object.getNelement(), object.getAddrData());     break;
                default:    ERROR(funcName, "unsupported bitPerSample value");
            }
            break;
        case SAMPLEFORMAT_INT:
            switch (bitPerSample) {
                case 8:     array_typecast((int8_t*)   bufferData, object.getNelement(), object.getAddrData());     break;
                case 16:    array_typecast((int16_t*)  bufferData, object.getNelement(), object.getAddrData());     break;
                case 32:    array_typecast((int32_t*)  bufferData, object.getNelement(), object.getAddrData());     break;
                case 64:    array_typecast((int64_t*)  bufferData, object.getNelement(), object.getAddrData());     break;
                default:    ERROR(funcName, "unsupported bitPerSample value");
            }
            break;
        case SAMPLEFORMAT_IEEEFP:
            switch (bitPerSample) {
                case 32:    array_typecast((float*)    bufferData, object.getNelement(), object.getAddrData());     break;
                case 64:    array_typecast((double*)   bufferData, object.getNelement(), object.getAddrData());     break;
                default:    ERROR(funcName, "unsupported bitPerSample value");
            }
            break;
        default:
            ERROR(funcName, "unsupported sampleFormat value");
    }

    array_delete(bufferData);

    TIFFClose(fileHandle);
}

template <typename T>
void cTIFF<T>::write(const cArray3<T>& object, const std::string& fileName, eNorm8bit norm8bit) const
{
    const std::string   funcName("void cTIFF::write(""const cArray3<T>& object, "
                                    "const std::string& fileName, eNorm8bit norm8bit) const");

    require(object.getAddrData()  != NULL, funcName + ": writing an empty object to a TIFF file");
    require(object.getDimension() == 2,    funcName + ": writing a non-image object to a TIFF file");

    checkExtension(fileName);

    TIFF*   fileHandle = TIFFOpen(fileName.c_str(), "w");

    require(fileHandle, funcName + ": could not open TIFF file " + fileName);

    TIFFSetField(fileHandle, TIFFTAG_IMAGELENGTH,     object.getNrow());
    TIFFSetField(fileHandle, TIFFTAG_IMAGEWIDTH,      object.getNcol());
    TIFFSetField(fileHandle, TIFFTAG_SAMPLEFORMAT,    SAMPLEFORMAT_UINT);
    TIFFSetField(fileHandle, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(fileHandle, TIFFTAG_BITSPERSAMPLE,   8);
    TIFFSetField(fileHandle, TIFFTAG_ORIENTATION,     ORIENTATION_TOPLEFT);
    TIFFSetField(fileHandle, TIFFTAG_PLANARCONFIG,    PLANARCONFIG_CONTIG);
    TIFFSetField(fileHandle, TIFFTAG_PHOTOMETRIC,     PHOTOMETRIC_MINISBLACK);
    TIFFSetField(fileHandle, TIFFTAG_PLANARCONFIG,    PLANARCONFIG_CONTIG);
    TIFFSetField(fileHandle, TIFFTAG_COMPRESSION,     COMPRESSION_DEFLATE);

    // convert data to 8bit
    std::vector<uint8_t>    dataTIFF;
    cArray3<T>              dataSave;

    switch (norm8bit) {
        case NORM8BIT_TRUE:
            dataSave = object;
            dataSave.normalize8bit();

            for(size_t i = 0; i < object.getNelement(); i++) {
                dataTIFF.push_back((uint8_t) dataSave[i]);
            }
            break;
        case NORM8BIT_FALSE:
            for(size_t i = 0; i < object.getNelement(); i++) {
                dataTIFF.push_back((uint8_t) object[i]);
            }
            break;
        default:
            ERROR(funcName, "unsupported norm8bit mode");
    }

    // write the information to the file
    TIFFWriteEncodedStrip(fileHandle, 0, &dataTIFF[0], (tsize_t) object.getNelement());

    TIFFClose(fileHandle);
}

// instantiation
template class cTIFF<int32_t >;
template class cTIFF<uint32_t>;
template class cTIFF<float   >;
template class cTIFF<double  >;

} // namespace gem
