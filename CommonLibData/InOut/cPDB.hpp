/***********************************************************************
 *  File:       cPDB.hpp
 *
 *  Purpose:    Header file for a PDB class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#ifndef __GEM_CPDB_HPP__
#define __GEM_CPDB_HPP__

#include "cVector3.hpp"
#include "cData3.hpp"
#include "cMapHeader.hpp"

#include <vector>

namespace gem {

enum ePDB2VolUseHetatm
{
    PDB2VOl_USE_HETATM_NO  = 0,
    PDB2VOl_USE_HETATM_YES = 1
};

/* Source: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

COLUMNS        DATA TYPE       CONTENTS
--------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  " / "HETATM"
 7 - 11        Integer         Atom serial number
13 - 16        Atom            Atom name
17             Character       Alternate location indicator
18 - 20        Residue name    Residue name
22             Character       Chain identifier
23 - 26        Integer         Residue sequence number
27             AChar           Code for insertion of residues
31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms
39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms
47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms
55 - 60        Real(6.2)       Occupancy
61 - 66        Real(6.2)       Temperature factor (Default = 0.0)
73 - 76        LString(4)      Segment identifier, left-justified
77 - 78        LString(2)      Element symbol, right-justified
79 - 80        LString(2)      Charge on the atom

*/

struct sAtom
{
    std::string         _recordName;            //  1 -  6
    int                 _serialNum;             //  7 - 11
    std::string         _atomName;              // 13 - 16
    std::string         _location;              //      17
    std::string         _residueName;           // 18 - 20
    std::string         _chainName;             //      22
    int                 _sequenceNum;           // 23 - 26
    std::string         _insertionCode;         //      27
    cVector3<double>    _coord;                 // 31 - 38
                                                // 39 - 46
                                                // 47 - 54
    double              _occupancy;             // 55 - 60
    double              _temperature;           // 61 - 66
    std::string         _segmentID;             // 73 - 76
    std::string         _elementSym;            // 77 - 78
    std::string         _charge;                // 79 - 80

    // internal use
    double              _mass;

    friend std::ostream& operator<<(std::ostream& os, const sAtom& atom);
};

class cPDB : public std::vector<sAtom>
{
protected:

public:
     cPDB(void) {};
    ~cPDB()     {};

    // empty checking
    void    requireEmpty   (const std::string& funcName) const { require( size() == 0, funcName + ": PDB is not empty"); }
    void    requireNonEmpty(const std::string& funcName) const { require( size()  > 0, funcName + ": PDB is empty");     }

    // operator overloading
    cPDB&                operator= (const cPDB& other);

    friend std::ostream& operator<<(std::ostream& os, const cPDB& pdb);

    // properties
    double  getMaxX(void) const;
    double  getMinX(void) const;
    double  getMaxY(void) const;
    double  getMinY(void) const;
    double  getMaxZ(void) const;
    double  getMinZ(void) const;

    // transformation
    void    opRotate       (const cVector3<double>& angle  );
    void    opScale        (const cVector3<double>& factor );
    void    opTranslate    (const cVector3<double>& offset );
    void    opTransRotTrans(const cVector3<double>& offset1,
                            const cVector3<double>& angle,
                            const cVector3<double>& offset2);

    void    opExtractChain(             const std::string& chains);
    void    opExtractChain(cPDB& other, const std::string& chains) const;
    void    opMergeAtom   (const std::vector<cPDB>& others);
    void    opMergeChain  (const std::vector<cPDB>& others);

    void    opPDB2Map    (double reso, double spacing, ePDB2VolUseHetatm hetatm,
                          cData3<float>& map, cMapHeader* header = NULL) const;
    double  opComputeRMSD(const cPDB& other, bool CAonly = true) const;

    double  opEstimateVolume(void) const;

    // input & output
    void    read       (const std::string& fileName);
    void    write      (const std::string& fileName) const;
    void    writeModels(const std::string& fileName,
                        const std::vector<cPDB>& others) const;
};

} // namespace gem

#endif
