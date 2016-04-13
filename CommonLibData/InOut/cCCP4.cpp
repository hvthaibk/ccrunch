/***********************************************************************
 *  File:       cCCP4.cpp
 *
 *  Purpose:    Implementation of a CCP4 class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cCCP4.hpp"

#include "edian.hpp"

#include "gpp4/ccp4/cmaplib.h"

namespace gem {

template <typename T>
void cCCP4<T>::checkExtension(const std::string& fileName) const
{
    const std::string   funcName("void cCCP4::checkExtension("
                                    "const std::string& fileName) const");

    if (fileName == "" || (fileName.find(".map") > fileName.length() &&
                           fileName.find(".MAP") > fileName.length())) {
        ERROR(funcName, "unsupported expansion (accept *.map, *.MAP)");
    }
}

template <typename T>
void cCCP4<T>::readHeader(const std::string& fileName, cMapHeader& header)
{
    const std::string   funcName("void cCCP4::readHeader("
                                    "const std::string& fileName, "
                                    "cMapHeader& header)");

    checkExtension(fileName);

    CMap_io::CMMFile *fileHandle = (CMap_io::CMMFile *)
                CMap_io::ccp4_cmap_open(fileName.c_str(), O_RDONLY);

    require(fileHandle, funcName + ": could not open CCP4 file " + fileName);

    // save header information
    header.file_name = fileName;

    header.data_mode  = fileHandle->data_mode;
    header.spacegroup = fileHandle->spacegroup;

    header.cell_dim[0] = fileHandle->cell[0];
    header.cell_dim[1] = fileHandle->cell[1];
    header.cell_dim[2] = fileHandle->cell[2];

    header.cell_angle[0] = fileHandle->cell[3];
    header.cell_angle[1] = fileHandle->cell[4];
    header.cell_angle[2] = fileHandle->cell[5];

    header.map_dim[0] = fileHandle->map_dim[0];
    header.map_dim[1] = fileHandle->map_dim[1];
    header.map_dim[2] = fileHandle->map_dim[2];

    header.origin[0] = fileHandle->origin[0];
    header.origin[1] = fileHandle->origin[1];
    header.origin[2] = fileHandle->origin[2];

    header.cell_grid[0] = fileHandle->cell_grid[0];
    header.cell_grid[1] = fileHandle->cell_grid[1];
    header.cell_grid[2] = fileHandle->cell_grid[2];

    header.axes_order[0] = fileHandle->axes_order[0];
    header.axes_order[1] = fileHandle->axes_order[1];
    header.axes_order[2] = fileHandle->axes_order[2];

    ccp4_cmap_get_mapstats(fileHandle, &header.min, &header.max, &header.mean, &header.rms);
    header.sym_byte = 80 * ccp4_cmap_num_symop(fileHandle);

    ccp4_cmap_close(fileHandle);
}

template <typename T>
void cCCP4<T>::readData(const std::string& fileName, cData3<T>& object, cMapHeader& header)
{
    const std::string   funcName("void cCCP4::readData(const std::string& fileName, "
                                    "cData3<T>& object, cMapHeader& header)");

    checkExtension(fileName);

    CMap_io::CMMFile *fileHandle = (CMap_io::CMMFile *)
                CMap_io::ccp4_cmap_open(fileName.c_str(), O_RDONLY);

    require(fileHandle, funcName + ": could not open CCP4 file " + fileName);

    // allocate new memory
    object.memReAlloc(cVector3<size_t>((size_t) fileHandle->map_dim[0],
                                       (size_t) fileHandle->map_dim[1],
                                       (size_t) fileHandle->map_dim[2]));

    // seek to data starting point
    ccp4_cmap_seek_data(fileHandle,
                        (CCP4_HEADER_SIZE + header.sym_byte) / 4,
                        SEEK_SET);
    /*ccp4_cmap_seek_data(fileHandle,
                        256+20*ccp4_cmap_num_symop(fileHandle),
                        SEEK_SET);*/

    int8_t      *dataInt8  = NULL;
    int16_t     *dataInt16 = NULL;
    float       *dataFloat = NULL;
    switch (fileHandle->data_mode) {
        case 0:            // int8
            array_new(dataInt8, object.getNelement());

            ccp4_cmap_read_data(fileHandle, (void *) dataInt8,
                                            (int)    object.getNelement());

            array_index_col2row(dataInt8, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            array_delete(dataInt8);

            break;
        case 1:            // int16
            array_new(dataInt16, object.getNelement());

            ccp4_cmap_read_data(fileHandle, (void *) dataInt16,
                                            (int)    object.getNelement());

            array_index_col2row(dataInt16, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            array_delete(dataInt16);

            break;
        case 2:            // float32
            array_new(dataFloat, object.getNelement());

            ccp4_cmap_read_data(fileHandle, (void *) dataFloat,
                                            (int)    object.getNelement());

            array_index_col2row(dataFloat, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            array_delete(dataFloat);

            break;
        default:
            ERROR(funcName, "unsupported data mode in reading " + fileName);
    }

    ccp4_cmap_close(fileHandle);
}

template <typename T>
void cCCP4<T>::read(const std::string& fileName, cData3<T>& object, cMapHeader* header)
{
    if (header) {
        readHeader(fileName, *header);
        readData  (fileName, object, *header);
    }
    else {
        cMapHeader  headerTmp;

        readHeader(fileName, headerTmp);
        readData  (fileName, object, headerTmp);
    }
}

template <typename T>
void cCCP4<T>::write(const cData3<T>& object, const std::string& fileName, cMapHeader* header)
{
    const std::string   funcName("void cCCP4::write(const cData3<T>& object, "
                                    "const std::string& fileName, cMapHeader* header)");

    require(object.getAddrData() != NULL, funcName + ":writing an empty object");

    checkExtension(fileName);

    CMap_io::CMMFile *fileHandle = (CMap_io::CMMFile *)
                CMap_io::ccp4_cmap_open(fileName.c_str(), O_WRONLY);

    require(fileHandle, funcName + ": could not open CCP4 file " + fileName);

    // write map header
    if (header) {
        if (header->map_dim[0] == (int) object.getNrow() &&
            header->map_dim[1] == (int) object.getNcol() &&
            header->map_dim[2] == (int) object.getNsec()) {

            fileHandle->data_mode  = header->data_mode;
            fileHandle->spacegroup = header->spacegroup;

            fileHandle->cell[0] = header->cell_dim[0];
            fileHandle->cell[1] = header->cell_dim[1];
            fileHandle->cell[2] = header->cell_dim[2];

            fileHandle->cell[3] = header->cell_angle[0];
            fileHandle->cell[4] = header->cell_angle[1];
            fileHandle->cell[5] = header->cell_angle[2];

            fileHandle->map_dim[0] = header->map_dim[0];
            fileHandle->map_dim[1] = header->map_dim[1];
            fileHandle->map_dim[2] = header->map_dim[2];

            fileHandle->origin[0] = header->origin[0];
            fileHandle->origin[1] = header->origin[1];
            fileHandle->origin[2] = header->origin[2];

            fileHandle->cell_grid[0] = header->cell_grid[0];
            fileHandle->cell_grid[1] = header->cell_grid[1];
            fileHandle->cell_grid[2] = header->cell_grid[2];

            fileHandle->axes_order[0] = header->axes_order[0];
            fileHandle->axes_order[1] = header->axes_order[1];
            fileHandle->axes_order[2] = header->axes_order[2];
        }
        else {
            ERROR(funcName, "header information does not match in writing " + fileName);
        }
    }
    else {
        fileHandle->map_dim[0] = (int) object.getNrow();
        fileHandle->map_dim[1] = (int) object.getNcol();
        fileHandle->map_dim[2] = (int) object.getNsec();

        fileHandle->data_mode  = 2;
        fileHandle->spacegroup = 0;

        fileHandle->cell[0] = (float) fileHandle->map_dim[0];
        fileHandle->cell[1] = (float) fileHandle->map_dim[1];
        fileHandle->cell[2] = (float) fileHandle->map_dim[2];

        fileHandle->cell[3] = 90;
        fileHandle->cell[4] = 90;
        fileHandle->cell[5] = 90;

        fileHandle->origin[0] = 0;
        fileHandle->origin[1] = 0;
        fileHandle->origin[2] = 0;

        fileHandle->cell_grid[0] = fileHandle->map_dim[0];
        fileHandle->cell_grid[1] = fileHandle->map_dim[1];
        fileHandle->cell_grid[2] = fileHandle->map_dim[2];

        fileHandle->axes_order[0] = 1;
        fileHandle->axes_order[1] = 2;
        fileHandle->axes_order[2] = 3;
    }

    // write map data, assume no symmetry information
    ccp4_cmap_seek_data(fileHandle, CCP4_HEADER_SIZE / 4, SEEK_SET);

    int8_t      *dataInt8  = NULL;
    int16_t     *dataInt16 = NULL;
    float       *dataFloat = NULL;
    switch (fileHandle->data_mode) {
        case 0:            // int8
            array_new(dataInt8, object.getNelement());
            array_index_row2col(object.getAddrData(), dataInt8,
                                object.getNrow(), object.getNcol(), object.getNsec());
            ccp4_cmap_write_data(fileHandle, (void *) dataInt8,
                                             (int)    object.getNelement());
            array_delete(dataInt8);
        case 1:            // int16
            array_new(dataInt16, object.getNelement());
            array_index_row2col(object.getAddrData(), dataInt16,
                                object.getNrow(), object.getNcol(), object.getNsec());
            ccp4_cmap_write_data(fileHandle, (void *) dataInt16,
                                             (int)    object.getNelement());
            array_delete(dataInt16);
        case 2:            // float32
            array_new(dataFloat, object.getNelement());
            array_index_row2col(object.getAddrData(), dataFloat,
                                object.getNrow(), object.getNcol(), object.getNsec());
            ccp4_cmap_write_data(fileHandle, (void *) dataFloat,
                                             (int)    object.getNelement());
            array_delete(dataFloat);
            break;
        default:
            ERROR(funcName, "unsupported data mode in writing " + fileName);
    }

    ccp4_cmap_close(fileHandle);
}

// instantiation
template class cCCP4<float >;
template class cCCP4<double>;

} // namespace gem
