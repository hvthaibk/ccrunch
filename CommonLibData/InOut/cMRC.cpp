/***********************************************************************
 *  File:       cMRC.cpp
 *
 *  Purpose:    Implementation of an MRC class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cMRC.hpp"

#include "edian.hpp"

namespace gem {

template <typename T>
void cMRC<T>::checkExtension(const std::string& fileName) const
{
    const std::string   funcName("void cMRC::checkExtension("
                                    "const std::string& fileName) const");

    if (fileName == "" || (fileName.find(".mrc") > fileName.length() &&
                           fileName.find(".MRC") > fileName.length())) {
        ERROR(funcName, "unsupported expansion (accept *.mrc, *.MRC)");
    }
}

template <typename T>
void cMRC<T>::readHeader(const std::string& fileName, cMapHeader& header)
{
    const std::string   funcName("void cMRC::readHeader("
                                    "const std::string& fileName, "
                                    "cMapHeader& header)");

    checkExtension(fileName);

    std::ifstream       fileHandle;

    fileHandle.open(fileName.c_str(), std::ifstream::binary);
    assure(fileHandle, fileName.c_str());

    int             map_dim[3], data_mode, origin[3], cell_grid[3], axes_order[3];
    float           cell_dim[3], cell_angle[3];
    float           density_min, density_max, density_mean;
    int             spacegroup, sym_byte;
    char            mapStamp[] = "    ";

    fileHandle.read((char*) map_dim,       3*sizeof(int));         // NX, NY, NZ
    fileHandle.read((char*) &data_mode,      sizeof(int));         // MODE
    fileHandle.read((char*) origin,        3*sizeof(int));         // NXSTART, NYSTART, NZSTART
    fileHandle.read((char*) cell_grid,     3*sizeof(int));         // MX, MY, MZ
    fileHandle.read((char*) cell_dim,      3*sizeof(float));       // CELLA
    fileHandle.read((char*) cell_angle,    3*sizeof(float));       // CELLB
    fileHandle.read((char*) axes_order,    3*sizeof(int));         // MAPC, MAPR, MAPS
    fileHandle.read((char*) &density_min,    sizeof(float));       // DMIN
    fileHandle.read((char*) &density_max,    sizeof(float));       // DMAX
    fileHandle.read((char*) &density_mean,   sizeof(float));       // DMEAN
    fileHandle.read((char*) &spacegroup,     sizeof(int));         // ISPG
    fileHandle.read((char*) &sym_byte,       sizeof(int));         // NSYMBT

    fileHandle.seekg(208, std::ios::beg);
    fileHandle.read((char*) mapStamp,      4*sizeof(char));        // MAP

    // check for the need of byte swapping
    header.edian_swap = false;
    if (map_dim[0] < 0 || map_dim[0] > 100000 ||
        map_dim[1] < 0 || map_dim[1] > 100000 ||
        map_dim[2] < 0 || map_dim[2] > 100000) {
        header.edian_swap = true;
    }

    if (header.edian_swap) {
        swap4_aligned(map_dim, 3);
        swap4_aligned(&data_mode, 1);
        swap4_aligned(origin, 3);
        swap4_aligned(cell_grid, 3);
        swap4_aligned(cell_dim, 3);
        swap4_aligned(cell_angle, 3);
        swap4_aligned(axes_order, 3);
        swap4_aligned(&density_min, 1);
        swap4_aligned(&density_max, 1);
        swap4_aligned(&density_mean, 1);
        swap4_aligned(&spacegroup, 1);
        swap4_aligned(&sym_byte, 1);
        swap4_aligned(mapStamp, 1);
    }

    // check data mode and 'MAP' stamp
    require(data_mode >=0 || data_mode <= 2,
            funcName + ": non-supported data type (data_mode != 0,1,2) in" + fileName);
    if (map_dim[2] != 1) {
        require(strcmp(mapStamp, "MAP ") == 0,
                funcName + ": cannot detect the 'MAP' string in " + fileName);
    }

    // check symmetry validity
    size_t      fileSize, dataOffset,
                nelement = (size_t) map_dim[0] * (size_t) map_dim[1] * (size_t) map_dim[2];

    fileHandle.seekg(0, std::ios::end);
    fileSize   = fileHandle.tellg();

    switch (data_mode) {
        case 0:            // int8
            dataOffset = fileSize - 1 * nelement;
            break;
        case 1:            // int16
            dataOffset = fileSize - 2 * nelement;
            break;
        case 2:            // float32
            dataOffset = fileSize - 4 * nelement;
            break;
        case 6:            // uint16
            dataOffset = fileSize - 2 * nelement;
            break;
        default:
            ERROR(funcName, "unsupported data mode in reading " + fileName);
    }

    if ((int) dataOffset != (CCP4_HEADER_SIZE + sym_byte)) {
        if (dataOffset < CCP4_HEADER_SIZE) {
            ERROR(funcName, "file " + fileName + " has size smaller than expected");
        }
        else if (dataOffset == CCP4_HEADER_SIZE) {
            ERROR(funcName, "the symmetry record in " + fileName + " is bogus");
        }
        else {
            ERROR(funcName, "file " + fileName + " has size larger than expected");
        }
    }

    // save header information
    header.file_name = fileName;

    header.data_mode  = data_mode;
    header.spacegroup = spacegroup;
    header.sym_byte   = sym_byte;

    header.cell_dim[0] = cell_dim[0];
    header.cell_dim[1] = cell_dim[1];
    header.cell_dim[2] = cell_dim[2];

    header.cell_angle[0] = cell_angle[0];
    header.cell_angle[1] = cell_angle[1];
    header.cell_angle[2] = cell_angle[2];

    header.map_dim[0] = map_dim[0];
    header.map_dim[1] = map_dim[1];
    header.map_dim[2] = map_dim[2];

    header.origin[0] = origin[0];
    header.origin[1] = origin[1];
    header.origin[2] = origin[2];

    header.cell_grid[0] = cell_grid[0];
    header.cell_grid[1] = cell_grid[1];
    header.cell_grid[2] = cell_grid[2];

    header.axes_order[0] = axes_order[0];
    header.axes_order[1] = axes_order[1];
    header.axes_order[2] = axes_order[2];

    header.min  = density_min;
    header.max  = density_max;
    header.mean = density_mean;

    fileHandle.close();
}

template <typename T>
void cMRC<T>::readData(const std::string& fileName, cData3<T>& object, cMapHeader& header)
{
    const std::string   funcName("void cMRC::readData(const std::string& fileName, "
                                    "cData3<T>& object, cMapHeader& header)");

    checkExtension(fileName);

    std::ifstream       fileHandle;

    fileHandle.open(fileName.c_str(), std::ifstream::binary);
    assure(fileHandle, fileName.c_str());

    // allocate new memory
    object.memReAlloc(cVector3<size_t>((size_t) header.map_dim[0],
                                       (size_t) header.map_dim[1],
                                       (size_t) header.map_dim[2]));

    // seek to data starting point
    fileHandle.seekg(CCP4_HEADER_SIZE + header.sym_byte, std::ios::beg);

    unsigned char*      bufferData   = NULL;
    size_t              bufferLength = (size_t) std::pow(2,header.data_mode) * header.getNelement();

    array_new(bufferData, bufferLength);
    fileHandle.read((char*) bufferData, bufferLength);

    switch (header.data_mode) {
        case 0:            // int8
            array_index_col2row((int8_t*) bufferData, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 1:            // int16
            if (header.edian_swap) { swap2_aligned(bufferData, header.getNelement()); }

            array_index_col2row((int16_t*) bufferData, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 2:            // float32
            if (header.edian_swap) { swap4_aligned(bufferData, header.getNelement()); }

            array_index_col2row((float*) bufferData, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 6:            // uint16
            if (header.edian_swap) { swap2_aligned(bufferData, header.getNelement()); }

            array_index_col2row((uint16_t*) bufferData, object.getAddrData(),
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        default:
            ERROR(funcName, "unsupported data mode in reading " + fileName);
    }

    array_delete(bufferData);

    fileHandle.close();
}

template <typename T>
void cMRC<T>::read(const std::string& fileName, cData3<T>& object, cMapHeader* header)
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
void cMRC<T>::write(const cData3<T>& object, const std::string& fileName, cMapHeader* header)
{
    const std::string funcName("void cMRC::write(const cData3<T>& object, "
                                    "const std::string& fileName, cMapHeader* header)");

    require(object.getAddrData() != NULL, funcName + "writing an empty object");

    checkExtension(fileName);

    std::ofstream       fileHandle;

    fileHandle.open(fileName.c_str(), std::ofstream::binary);
    assure(fileHandle, fileName.c_str());

    int             map_dim[3], data_mode, origin[3], cell_grid[3], axes_order[3];
    float           cell_dim[3], cell_angle[3];
    float           density_min, density_max, density_mean;
    int             spacegroup, sym_byte = 0;
    char            mapStamp[] = "MAP ";

    // write map header
    if (header) {
        if (header->map_dim[0] == (int) object.getNrow() &&
            header->map_dim[1] == (int) object.getNcol() &&
            header->map_dim[2] == (int) object.getNsec()) {

            data_mode  = header->data_mode;
            spacegroup = header->spacegroup;

            cell_dim[0] = header->cell_dim[0];
            cell_dim[1] = header->cell_dim[1];
            cell_dim[2] = header->cell_dim[2];

            cell_angle[0] = header->cell_angle[0];
            cell_angle[1] = header->cell_angle[1];
            cell_angle[2] = header->cell_angle[2];

            map_dim[0] = header->map_dim[0];
            map_dim[1] = header->map_dim[1];
            map_dim[2] = header->map_dim[2];

            origin[0] = header->origin[0];
            origin[1] = header->origin[1];
            origin[2] = header->origin[2];

            cell_grid[0] = header->cell_grid[0];
            cell_grid[1] = header->cell_grid[1];
            cell_grid[2] = header->cell_grid[2];

            axes_order[0] = header->axes_order[0];
            axes_order[1] = header->axes_order[1];
            axes_order[2] = header->axes_order[2];
        }
        else {
            ERROR(funcName, "header information does not match in writing " + fileName);
        }
    }
    else {
        map_dim[0] = (int) object.getNrow();
        map_dim[1] = (int) object.getNcol();
        map_dim[2] = (int) object.getNsec();

        data_mode  = 2;
        spacegroup = 0;

        cell_dim[0] = (float) map_dim[0];
        cell_dim[1] = (float) map_dim[1];
        cell_dim[2] = (float) map_dim[2];

        cell_angle[0] = 90;
        cell_angle[1] = 90;
        cell_angle[2] = 90;

        origin[0] = 0;
        origin[1] = 0;
        origin[2] = 0;

        cell_grid[0] = map_dim[0];
        cell_grid[1] = map_dim[1];
        cell_grid[2] = map_dim[2];

        axes_order[0] = 1;
        axes_order[1] = 2;
        axes_order[2] = 3;
    }

    fileHandle.write((char*) map_dim,       3*sizeof(int));         // NX, NY, NZ
    fileHandle.write((char*) &data_mode,      sizeof(int));         // MODE
    fileHandle.write((char*) origin,        3*sizeof(int));         // NXSTART, NYSTART, NZSTART
    fileHandle.write((char*) cell_grid,     3*sizeof(int));         // MX, MY, MZ
    fileHandle.write((char*) cell_dim,      3*sizeof(float));       // CELLA
    fileHandle.write((char*) cell_angle,    3*sizeof(float));       // CELLB
    fileHandle.write((char*) axes_order,    3*sizeof(int));         // MAPC, MAPR, MAPS
    fileHandle.write((char*) &density_min,    sizeof(float));       // DMIN
    fileHandle.write((char*) &density_max,    sizeof(float));       // DMAX
    fileHandle.write((char*) &density_mean,   sizeof(float));       // DMEAN
    fileHandle.write((char*) &spacegroup,     sizeof(int));         // ISPG
    fileHandle.write((char*) &sym_byte,       sizeof(int));         // NSYMBT

    fileHandle.seekp(208, std::ios::beg);
    fileHandle.write((char*) mapStamp,      4*sizeof(char));        // MAP

    // write map data, assume no symmetry information
    fileHandle.seekp(CCP4_HEADER_SIZE, std::ios::beg);

    unsigned char*      bufferData   = NULL;
    size_t              bufferLength,
                        nelement     = (size_t) (map_dim[0]*map_dim[1]*map_dim[2]);

    switch (data_mode) {
        case 0:            // int8
            bufferLength = 1 * nelement;
            break;
        case 1:            // int16
            bufferLength = 2 * nelement;
            break;
        case 2:            // float32
            bufferLength = 4 * nelement;
            break;
        case 6:            // uint16
            bufferLength = 2 * nelement;
            break;
        default:
            ERROR(funcName, "unsupported data mode in reading " + fileName);
    }

    array_new(bufferData, bufferLength);

    switch (data_mode) {
        case 0:            // int8
            array_index_row2col(object.getAddrData(), (int8_t*) bufferData,
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 1:            // int16
            array_index_row2col(object.getAddrData(), (int16_t*) bufferData,
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 2:            // float32
            array_index_row2col(object.getAddrData(), (float*) bufferData,
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        case 6:            // uint16
            array_index_row2col(object.getAddrData(), (uint16_t*) bufferData,
                                object.getNrow(), object.getNcol(), object.getNsec());

            break;
        default:
            ERROR(funcName, "unsupported data mode in writing " + fileName);
    }

    fileHandle.write((char*) bufferData, bufferLength);

    array_delete(bufferData);

    fileHandle.close();
}

// instantiation
template class cMRC<float >;
template class cMRC<double>;

} // namespace gem
