/***********************************************************************
 *  File:       cMapHeader.hpp
 *
 *  Purpose:    Header file for a map header structure
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CMAPHEADER_HPP__
#define __GEM_CMAPHEADER_HPP__

#include "cVector3.hpp"

#include <cstddef>

namespace gem {

/*

 1      NC              # of Columns    (fastest changing in map)
 2      NR              # of Rows
 3      NS              # of Sections   (slowest changing in map)
 4      MODE            Data type
                          0 = envelope stored as signed bytes (from
                              -128 lowest to 127 highest)
                          1 = Image     stored as Integer*2
                          2 = Image     stored as Reals
                          3 = Transform stored as Complex Integer*2
                          4 = Transform stored as Complex Reals
                          5 == 0

                          Note: Mode 2 is the normal mode used in
                                the CCP4 programs. Other modes than 2 and 0
                                may NOT WORK

 5      NCSTART         Number of first COLUMN  in map
 6      NRSTART         Number of first ROW     in map
 7      NSSTART         Number of first SECTION in map
 8      NX              Number of intervals along X
 9      NY              Number of intervals along Y
10      NZ              Number of intervals along Z
11      X length        Cell Dimensions (Angstroms)
12      Y length                     "
13      Z length                     "
14      Alpha           Cell Angles     (Degrees)
15      Beta                         "
16      Gamma                        "
17      MAPC            Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
18      MAPR            Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
19      MAPS            Which axis corresponds to Sects. (1,2,3 for X,Y,Z)
20      AMIN            Minimum density value
21      AMAX            Maximum density value
22      AMEAN           Mean    density value    (Average)
23      ISPG            Space group number
24      NSYMBT          Number of bytes used for storing symmetry operators
25      LSKFLG          Flag for skew transformation, =0 none, =1 if foll
26-34   SKWMAT          Skew matrix S (in order S11, S12, S13, S21 etc) if
                        LSKFLG .ne. 0.
35-37   SKWTRN          Skew translation t if LSKFLG .ne. 0.
                        Skew transformation is from standard orthogonal
                        coordinate frame (as used for atoms) to orthogonal
                        map frame, as

                                Xo(map) = S * (Xo(atoms) - t)

38      future use       (some of these are used by the MSUBSX routines
 .          "              in MAPBRICK, MAPCONT and FRODO)
 .          "   (all set to zero by default)
 .          "
52          "

53      MAP             Character string 'MAP ' to identify file type
54      MACHST          Machine stamp indicating the machine type
                        which wrote file
55      ARMS            Rms deviation of map from mean density
56      NLABL           Number of labels being used
57-256  LABEL(20,10)    10  80 character text labels (ie. A4 format)

*/

#define CCP4_HEADER_SIZE    1024

class cMapHeader        // mimic _CMMFile from gpp4
{
public:
    std::string         file_name;
    int                 data_mode;
    cVector3<int>       map_dim;
    cVector3<int>       origin;
    cVector3<int>       cell_grid;
    cVector3<int>       axes_order;
    cVector3<float>     cell_dim;
    cVector3<float>     cell_angle;
    int                 spacegroup;
    float               min;
    float               max;
    double              mean;       // gpp4 use double instead of float
    double              rms;        // gpp4 use double instead of float
    int                 sym_byte;
    bool                edian_swap;

    size_t              getNelement    (void) const;
    cVector3<size_t>    getMapDim      (void) const;
    cVector3<ptrdiff_t> getMapOrigin   (void) const;
    double              getVoxelSpacing(void) const;
    double              getVoxelVolumn (void) const;

    void                updateAfterExtension(const cVector3<size_t>& newsize);
    void                updateAfterPadding  (const cVector3<size_t>& newsize, const cVector3<size_t>& offset);
    void                updateAfterCropping (const cVector3<size_t>& newsize, const cVector3<size_t>& offset);
    void                updateAfterScaling  (float factor);

    void                simulate(const cVector3<size_t>&    size,
                                 const cVector3<ptrdiff_t>& pos,
                                 float spacing);
    void                print   (void) const;
};

} // namespace gem

#endif
