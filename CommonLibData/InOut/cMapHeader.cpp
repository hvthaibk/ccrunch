/***********************************************************************
 *  File:       cMapHeader.cpp
 *
 *  Purpose:    Implementation of a map header structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cMapHeader.hpp"
#include "transform.hpp"

#include <iostream>

namespace gem {

size_t cMapHeader::getNelement(void) const
{
    return (size_t) map_dim[0] * (size_t) map_dim[1] * (size_t) map_dim[2];
}

cVector3<size_t> cMapHeader::getMapDim(void) const
{
    return cVector3<size_t>((size_t) map_dim[0],
                            (size_t) map_dim[1],
                            (size_t) map_dim[2]);
}

cVector3<ptrdiff_t> cMapHeader::getMapOrigin(void) const
{
    return cVector3<ptrdiff_t>((ptrdiff_t) origin[0],
                               (ptrdiff_t) origin[1],
                               (ptrdiff_t) origin[2]);
}

double cMapHeader::getVoxelSpacing(void) const
{
    const std::string   funcName("double cMapHeader::getVoxelSpacing(void) const");

    cVector3<double>    voxelSpacing((double) cell_dim[0] / (double) map_dim[0],
                                     (double) cell_dim[1] / (double) map_dim[1],
                                     (double) cell_dim[2] / (double) map_dim[2]);

    require(std::abs(voxelSpacing[0]-voxelSpacing[1]) < GEM_FLOATING_TOLERANCE &&
            std::abs(voxelSpacing[0]-voxelSpacing[2]) < GEM_FLOATING_TOLERANCE,
            funcName + ": voxel spacing in three dimensions are not equal");

    return voxelSpacing[0];
}

double cMapHeader::getVoxelVolumn(void) const
{
    return (double) cell_dim.getProduct() / (double) map_dim.getProduct();
}

void cMapHeader::simulate(const cVector3<size_t>&    size,
                          const cVector3<ptrdiff_t>& pos,
                          float spacing)
{
    data_mode  = 2;
    map_dim    = size;
    origin     = pos;
    cell_grid  = map_dim;
    axes_order = cVector3<int>(1,2,3);
    cell_dim   = map_dim;
    cell_dim  *= spacing;
    cell_angle = cVector3<float>(90,90,90);
    spacegroup = 0;
}

void cMapHeader::updateAfterExtension(const cVector3<size_t>& newsize)
{
    float   spacing = (float) getVoxelSpacing();

    for (size_t i = 0; i < 3; i++) {
        origin[i]    = origin[i] - ((int) newsize[i] - map_dim[i]) / 2;
        map_dim[i]   = (int) newsize[i];
        cell_grid[i] = map_dim[i];
        cell_dim[i]  = (float) map_dim[i] * spacing;
    }
}

void cMapHeader::updateAfterPadding(const cVector3<size_t>& newsize, const cVector3<size_t>& offset)
{
    float   spacing = (float) getVoxelSpacing();

    for (size_t i = 0; i < 3; i++) {
        origin[i]    = origin[i] - (int) offset[i];
        map_dim[i]   = (int) newsize[i];
        cell_grid[i] = map_dim[i];
        cell_dim[i]  = (float) map_dim[i] * spacing;
    }
}

void cMapHeader::updateAfterCropping(const cVector3<size_t>& newsize, const cVector3<size_t>& offset)
{
    float   spacing = (float) getVoxelSpacing();

    for (size_t i = 0; i < 3; i++) {
        origin[i]    = origin[i] + (int) offset[i];
        map_dim[i]   = (int) newsize[i];
        cell_grid[i] = map_dim[i];
        cell_dim[i]  = (float) map_dim[i] * spacing;
    }
}

void cMapHeader::updateAfterScaling(float factor)
{
    const std::string   funcName("void cMapHeader::updateAfterScaling(float factor)");

    float   spacing = (float) getVoxelSpacing() / factor;

    for (size_t i = 0; i < 3; i++) {
        origin[i]    = (int) transform_centerCoord((float) origin[i], factor);
        map_dim[i]   = (int) transform_scaleSize((size_t) map_dim[i], factor);
        cell_grid[i] = map_dim[i];
        cell_dim[i]  = (float) map_dim[i] * spacing;
    }

    WARNING(funcName, "this functions is not perfect due to the problem of placing map's origin after scaling, use with precaution");
}

void cMapHeader::print(void) const
{
    std::cout << "\n";
    std::cout << "Map HEADER INFO: "     << file_name  << "\n";
    std::cout << "     Data mode     : " << data_mode  << "\n";
    std::cout << "     Map dimension : " << map_dim    << "\n";
    std::cout << "     Map origin    : " << origin     << "\n";
    std::cout << "     Cell dimension: " << cell_dim   << "\n";
    std::cout << "     Cell angle    : " << cell_angle << "\n";
    std::cout << "     Cell grid     : " << cell_grid  << "\n";
    std::cout << "     Axis order    : " << axes_order << "\n";
    std::cout << "     Space group   : " << spacegroup << "\n";
    std::cout << "     Min           : " << min        << "\n";
    std::cout << "     Max           : " << max        << "\n";
    std::cout << "     Mean          : " << mean       << "\n";
    std::cout << "     RMS           : " << rms        << "\n";
    std::cout << "     Symmetry space: " << sym_byte   << "\n";
    std::cout << "     Byte swapping : " << edian_swap << "\n";
    std::cout << "\n";
}

} // namespace gem
