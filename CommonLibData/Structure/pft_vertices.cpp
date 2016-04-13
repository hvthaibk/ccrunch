/***********************************************************************
 *  File:       pft_vertices.cpp
 *
 *  Purpose:    Implementation of vertices functions
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "macro.hpp"
#include "pft_vertices.hpp"

#include <cstdlib>
#include <iostream>
#include <iomanip>

namespace gem {

template <class T>
void pft_angle(size_t index, size_t nNode, size_t nRotZ, cVector3<T>& angle)
{
    const float     *vertices = pft_vertices_getref(nNode);

    size_t  iNode = index / nRotZ;
    size_t  iRotZ = index - iNode*nRotZ;

    angle[0] = (T) vertices[2*iNode];
    angle[1] = (T) vertices[2*iNode+1];
    angle[2] = (T) iRotZ * (T) 2 * (T) M_PI / (T) nRotZ;
}

// instantiation
template
void pft_angle<float >(size_t index, size_t nNode, size_t nRotZ, cVector3<float >& angle);
template
void pft_angle<double>(size_t index, size_t nNode, size_t nRotZ, cVector3<double>& angle);

float pft_vertices_getangle(size_t nNode)
{

    switch (nNode) {
        case 1:     return pft_angle_1;       break;
        case 12:    return pft_angle_12;      break;
        case 20:    return pft_angle_20;      break;
        case 32:    return pft_angle_32;      break;
        case 42:    return pft_angle_42;      break;
        case 60:    return pft_angle_60;      break;
        case 80:    return pft_angle_80;      break;
        case 92:    return pft_angle_92;      break;
        case 162:   return pft_angle_162;     break;
        case 180:   return pft_angle_180;     break;
        case 252:   return pft_angle_252;     break;
        case 320:   return pft_angle_320;     break;
        case 362:   return pft_angle_362;     break;
        case 492:   return pft_angle_492;     break;
        case 500:   return pft_angle_500;     break;
        case 642:   return pft_angle_642;     break;
        case 720:   return pft_angle_720;     break;
        case 980:   return pft_angle_980;     break;
        case 1280:  return pft_angle_1280;    break;
        case 2562:  return pft_angle_2562;    break;
        case 10242: return pft_angle_10242;   break;
        case 40962: return pft_angle_40962;   break;
        default:    ERROR("pft_vertices_getangle", "unsupported nNode");
    }
}

void pft_vertices_check(size_t nNode)
{
    switch (nNode) {
        case 1:     break;
        case 12:    break;
        case 20:    break;
        case 32:    break;
        case 42:    break;
        case 60:    break;
        case 80:    break;
        case 92:    break;
        case 162:   break;
        case 180:   break;
        case 252:   break;
        case 320:   break;
        case 362:   break;
        case 492:   break;
        case 500:   break;
        case 642:   break;
        case 720:   break;
        case 980:   break;
        case 1280:  break;
        case 2562:  break;
        case 10242: break;
        case 40962: break;
        default:
            std::cerr << "\nnNode = " << nNode << " is unsupported! \n"
                      << std::setprecision(2) << std::fixed
                      << "The supporteds are:     1 (" << std::setw(6) << pft_angle_1     << " degree) \n"
                      << "                       12 (" << std::setw(6) << pft_angle_12    << " degree) \n"
                      << "                       32 (" << std::setw(6) << pft_angle_32    << " degree) \n"
                      << "                       42 (" << std::setw(6) << pft_angle_42    << " degree) \n"
                      << "                       60 (" << std::setw(6) << pft_angle_60    << " degree) \n"
                      << "                       80 (" << std::setw(6) << pft_angle_80    << " degree) \n"
                      << "                       92 (" << std::setw(6) << pft_angle_92    << " degree) \n"
                      << "                      162 (" << std::setw(6) << pft_angle_162   << " degree) \n"
                      << "                      180 (" << std::setw(6) << pft_angle_180   << " degree) \n"
                      << "                      252 (" << std::setw(6) << pft_angle_252   << " degree) \n"
                      << "                      320 (" << std::setw(6) << pft_angle_320   << " degree) \n"
                      << "                      362 (" << std::setw(6) << pft_angle_362   << " degree) \n"
                      << "                      492 (" << std::setw(6) << pft_angle_492   << " degree) \n"
                      << "                      500 (" << std::setw(6) << pft_angle_500   << " degree) \n"
                      << "                      642 (" << std::setw(6) << pft_angle_642   << " degree) \n"
                      << "                      720 (" << std::setw(6) << pft_angle_720   << " degree) \n"
                      << "                      980 (" << std::setw(6) << pft_angle_980   << " degree) \n"
                      << "                     1280 (" << std::setw(6) << pft_angle_1280  << " degree) \n"
                      << "                     2562 (" << std::setw(6) << pft_angle_2562  << " degree) \n"
                      << "                    10242 (" << std::setw(6) << pft_angle_10242 << " degree) \n"
                      << "                    40962 (" << std::setw(6) << pft_angle_40962 << " degree) \n\n";
            exit(EXIT_FAILURE);
    }
}

const float* pft_vertices_getref(size_t nNode)
{
    const float     *vertices = NULL;

    switch (nNode) {
        case 1:     vertices = pft_vertices_1;       break;
        case 12:    vertices = pft_vertices_12;      break;
        case 20:    vertices = pft_vertices_20;      break;
        case 32:    vertices = pft_vertices_32;      break;
        case 42:    vertices = pft_vertices_42;      break;
        case 60:    vertices = pft_vertices_60;      break;
        case 80:    vertices = pft_vertices_80;      break;
        case 92:    vertices = pft_vertices_92;      break;
        case 162:   vertices = pft_vertices_162;     break;
        case 180:   vertices = pft_vertices_180;     break;
        case 252:   vertices = pft_vertices_252;     break;
        case 320:   vertices = pft_vertices_320;     break;
        case 362:   vertices = pft_vertices_362;     break;
        case 492:   vertices = pft_vertices_492;     break;
        case 500:   vertices = pft_vertices_500;     break;
        case 642:   vertices = pft_vertices_642;     break;
        case 720:   vertices = pft_vertices_720;     break;
        case 980:   vertices = pft_vertices_980;     break;
        case 1280:  vertices = pft_vertices_1280;    break;
        case 2562:  vertices = pft_vertices_2562;    break;
        case 10242: vertices = pft_vertices_10242;   break;
        case 40962: vertices = pft_vertices_40962;   break;
        default:    ERROR("pft_vertices", "unsupported nNode");
    }

    return vertices;
}

} // namespace gem
