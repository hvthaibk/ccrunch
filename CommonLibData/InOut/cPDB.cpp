/***********************************************************************
 *  File:       cPDB.cpp
 *
 *  Purpose:    Implementation of a PDB class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cPDB.hpp"

#include "cFilterSpatial.hpp"
#include "cMatrix33.hpp"
#include "cFileSystem.hpp"

namespace gem {

const char* atom_name[] = {
    " H  ",  " HE ",  " LI ",  " BE ",  " B  ",  " C  ",  " N  ",  " O  ",  " F  ",  " NE ",
    " NA ",  " MG ",  " AL ",  " SI ",  " P  ",  " S  ",  " CL ",  " AR ",  " K  ",  " CA ",
    " SC ",  " TI ",  " V  ",  " CR ",  " MN ",  " FE ",  " CO ",  " NI ",  " CU ",  " ZN "
};

const double atom_mass[30] = {
     1.008,  4.003, 6.940,   9.012, 10.810, 12.011, 14.007, 16.000, 18.998, 20.179,
    22.990, 24.305, 26.982, 28.086, 30.974, 32.060, 35.453, 39.948, 39.098, 40.080,
    44.956, 47.880, 50.942, 51.996, 54.938, 55.847, 58.933, 58.710, 63.546, 65.380
};

const size_t atom_frequent[6] = {0, 5, 6, 7, 14, 15};

std::ostream& operator<<(std::ostream& os, const sAtom& atom)
{
    std::cout << std::fixed
              << atom._recordName
              << std::setw(5) << atom._serialNum
              << " "
              << atom._atomName
              << atom._location
              << atom._residueName
              << " "
              << atom._chainName
              << std::setw(4) << atom._sequenceNum
              << atom._insertionCode
              << "   "
              << std::setw(8) << std::setprecision(3) << atom._coord[0]
              << std::setw(8) << std::setprecision(3) << atom._coord[1]
              << std::setw(8) << std::setprecision(3) << atom._coord[2]
              << std::setw(6) << std::setprecision(2) << atom._occupancy
              << std::setw(6) << std::setprecision(2) << atom._temperature
              << "      "
              << atom._segmentID
              << atom._elementSym
              << atom._charge
              << std::setw(8) << std::setprecision(3) << atom._mass
              << std::endl;

    return os;
}

cPDB& cPDB::operator=(const cPDB& other)
{
    const std::string   funcName("cPDB& cPDB::operator=(const cPDB& other)");

    other.requireNonEmpty(funcName);

    resize(other.size());

    #pragma omp parallel for
    for (size_t i = 0; i < other.size(); i++) {
        at(i) = other.at(i);
    }

    return *this;
}

std::ostream& operator<<(std::ostream& os, const cPDB& pdb)
{
    const std::string   funcName("std::ostream& operator<<(std::ostream& os, const cPDB& pdb)");

    pdb.requireNonEmpty(funcName);

    for (size_t i = 0; i < pdb.size(); i++) {
        std::cout << pdb.at(i);
    }

    return os;
}

double cPDB::getMaxX(void) const
{
    const std::string   funcName("double cPDB::getMaxX(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[0];

    for (size_t i = 1; i < size(); i++) {
        value = (value < at(i)._coord[0]) ? at(i)._coord[0] : value;
    }

    return value;
}

double cPDB::getMinX(void) const
{
    const std::string   funcName("double cPDB::getMinX(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[0];

    for (size_t i = 1; i < size(); i++) {
        value = (value > at(i)._coord[0]) ? at(i)._coord[0] : value;
    }

    return value;
}

double cPDB::getMaxY(void) const
{
    const std::string   funcName("double cPDB::getMaxY(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[1];

    for (size_t i = 1; i < size(); i++) {
        value = (value < at(i)._coord[1]) ? at(i)._coord[1] : value;
    }

    return value;
}

double cPDB::getMinY(void) const
{
    const std::string   funcName("double cPDB::getMinY(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[1];

    for (size_t i = 1; i < size(); i++) {
        value = (value > at(i)._coord[1]) ? at(i)._coord[1] : value;
    }

    return value;
}

double cPDB::getMaxZ(void) const
{
    const std::string   funcName("double cPDB::getMaxZ(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[2];

    for (size_t i = 1; i < size(); i++) {
        value = (value < at(i)._coord[2]) ? at(i)._coord[2] : value;
    }

    return value;
}

double cPDB::getMinZ(void) const
{
    const std::string   funcName("double cPDB::getMinZ(void) const");

    requireNonEmpty(funcName);

    double value = at(0)._coord[2];

    for (size_t i = 1; i < size(); i++) {
        value = (value > at(i)._coord[2]) ? at(i)._coord[2] : value;
    }

    return value;
}

void cPDB::opRotate(const cVector3<double>& angle)
{
    const std::string   funcName("void cPDB::opRotate(const cVector3<double>& angle)");

    requireNonEmpty(funcName);

    cVector3<double>    point3;
    cMatrix33<double>   matRot;

    transform_rotmat(angle[0], angle[1], angle[2],
                     matRot[0][0], matRot[0][1], matRot[0][2],
                     matRot[1][0], matRot[1][1], matRot[1][2],
                     matRot[2][0], matRot[2][1], matRot[2][2],
                     ROT3D_RIGHT_ZYZ);

    #pragma omp parallel for private(point3)
    for (size_t i = 0; i < size(); i++) {
        point3       = matRot * at(i)._coord;
        at(i)._coord = point3;
    }
}

void cPDB::opScale(const cVector3<double>& factor)
{
    const std::string   funcName("void cPDB::opScale(const cVector3<double>& factor)");

    requireNonEmpty(funcName);
    require((factor[0] != 0) && (factor[1] != 0) && (factor[2] != 0),
            funcName + ": scale with zero factor");

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
        at(i)._coord /= factor;
    }
}

void cPDB::opTranslate(const cVector3<double>& offset)
{
    const std::string   funcName("void cPDB::opTranslate(const cVector3<double>& offset)");

    requireNonEmpty(funcName);

    #pragma omp parallel for
    for (size_t i = 0; i < size(); i++) {
        at(i)._coord += offset;
    }
}

void cPDB::opTransRotTrans(const cVector3<double>& offset1,
                           const cVector3<double>& angle,
                           const cVector3<double>& offset2)
{
    const std::string   funcName("void cPDB::opTransRotTrans("
                                    "const cVector3<double>& offset1, "
                                    "const cVector3<double>& angle, "
                                    "const cVector3<double>& offset2)");

    requireNonEmpty(funcName);

    cVector3<double>    point3;
    cMatrix33<double>   matRot;

    transform_rotmat(angle[0], angle[1], angle[2],
                     matRot[0][0], matRot[0][1], matRot[0][2],
                     matRot[1][0], matRot[1][1], matRot[1][2],
                     matRot[2][0], matRot[2][1], matRot[2][2],
                     ROT3D_RIGHT_ZYZ);

    #pragma omp parallel for private(point3)
    for (size_t i = 0; i < size(); i++) {
        point3       = matRot * (at(i)._coord + offset1);
        at(i)._coord = point3 + offset2;
    }
}


void cPDB::opExtractChain(const std::string& chains)
{
    const std::string   funcName("void cPDB::opExtractChain(const std::string& chains)");

    requireNonEmpty(funcName);

    cPDB        other;

    for (size_t i = 0; i < size(); i++) {
        for (size_t c = 0; c < chains.size(); c++) {
            if (at(i)._chainName[0] == chains[c]) {
                other.push_back(at(i));
            }
        }
    }

    *this = other;
}

void cPDB::opExtractChain(cPDB& other, const std::string& chains) const
{
    const std::string   funcName("void cPDB::opExtractChain(cPDB& other, const std::string& chains) const");

    requireNonEmpty(funcName);

    other.clear();

    for (size_t i = 0; i < size(); i++) {
        for (size_t c = 0; c < chains.size(); c++) {
            if (at(i)._chainName[0] == chains[c]) {
                other.push_back(at(i));
            }
        }
    }
}

void cPDB::opMergeAtom(const std::vector<cPDB>& others)
{
    const std::string   funcName("void cPDB::opMergeAtom(const std::vector<cPDB>& others)");

    requireEmpty(funcName);
    require(others.size() > 1, funcName + ": no need to merge");

    for (size_t i = 0; i < others.size(); i++) {

        require(others[i].size() > 0, funcName + ": PDB is empty");

        for (size_t j = 0; j < others[i].size(); j++) {
            push_back(others[i].at(j));
        }
    }
}

void cPDB::opMergeChain(const std::vector<cPDB>& others)
{
    const std::string   funcName("void cPDB::opMergeChain(const std::vector<cPDB>& others)");

    requireEmpty(funcName);
    require(others.size() > 1, funcName + ": no need to merge");

    sAtom       atom;

    for (size_t i = 0; i < others.size(); i++) {

        require(others[i].size() > 0, funcName + ": PDB is empty");

        for (size_t j = 0; j < others[i].size(); j++) {

            atom = others[i].at(j);
            atom._chainName = char(65+i);

            push_back(atom);
        }
    }
}

void cPDB::opPDB2Map(double reso, double spacing, ePDB2VolUseHetatm hetatm,
                     cData3<float>& map, cMapHeader* header) const
{
    const std::string   funcName("void cPDB::opPDB2Map(double reso, double spacing, ePDB2VolUseHetatm hetatm, cData3<float>& map, cMapHeader* header) const");

    requireNonEmpty(funcName);
    require(reso > 0.1 && spacing > 0.1,
            funcName + ": incorrect input params");

    double      mapMaxX = spacing * std::ceil (getMaxX() / spacing),
                mapMinX = spacing * std::floor(getMinX() / spacing),
                mapMaxY = spacing * std::ceil (getMaxY() / spacing),
                mapMinY = spacing * std::floor(getMinY() / spacing),
                mapMaxZ = spacing * std::ceil (getMaxZ() / spacing),
                mapMinZ = spacing * std::floor(getMinZ() / spacing);

    // project atoms into lattice
    size_t      mapMargin = 2;
    double      mapX, mapY, mapZ, fracX, fracY, fracZ;
    size_t      intX, intY, intZ, intX1, intY1, intZ1;

    map.memReAllocZero(cVector3<size_t>((size_t) std::ceil((mapMaxX-mapMinX)/spacing) + 2*mapMargin + 1,
                                        (size_t) std::ceil((mapMaxY-mapMinY)/spacing) + 2*mapMargin + 1,
                                        (size_t) std::ceil((mapMaxZ-mapMinZ)/spacing) + 2*mapMargin + 1));

    size_t      mapNcol = map.getNcol();
    size_t      mapNsec = map.getNsec();

    for (size_t i = 0; i < size(); i++) {

        if (at(i)._recordName == "ATOM  " ||
           (at(i)._recordName == "HETATM" && hetatm == PDB2VOl_USE_HETATM_YES)) {

            mapX = (double) mapMargin + (at(i)._coord[0] - mapMinX) / spacing;
            mapY = (double) mapMargin + (at(i)._coord[1] - mapMinY) / spacing;
            mapZ = (double) mapMargin + (at(i)._coord[2] - mapMinZ) / spacing;

            intX = (size_t) std::floor(mapX);
            intY = (size_t) std::floor(mapY);
            intZ = (size_t) std::floor(mapZ);

            intX1 = intX + 1;
            intY1 = intY + 1;
            intZ1 = intZ + 1;

            fracX = mapX - (double) intX;
            fracY = mapY - (double) intY;
            fracZ = mapZ - (double) intZ;

            map[sub2ind(intX ,intY ,intZ ,mapNcol,mapNsec)] += (float) (at(i)._mass * (1-fracX) * (1-fracY) * (1-fracZ));
            map[sub2ind(intX ,intY ,intZ1,mapNcol,mapNsec)] += (float) (at(i)._mass * (1-fracX) * (1-fracY) * fracZ    );
            map[sub2ind(intX ,intY1,intZ ,mapNcol,mapNsec)] += (float) (at(i)._mass * (1-fracX) * fracY     * (1-fracZ));
            map[sub2ind(intX ,intY1,intZ1,mapNcol,mapNsec)] += (float) (at(i)._mass * (1-fracX) * fracY     * fracZ    );
            map[sub2ind(intX1,intY ,intZ ,mapNcol,mapNsec)] += (float) (at(i)._mass * fracX     * (1-fracY) * (1-fracZ));
            map[sub2ind(intX1,intY ,intZ1,mapNcol,mapNsec)] += (float) (at(i)._mass * fracX     * (1-fracY) * fracZ    );
            map[sub2ind(intX1,intY1,intZ ,mapNcol,mapNsec)] += (float) (at(i)._mass * fracX     * fracY     * (1-fracZ));
            map[sub2ind(intX1,intY1,intZ1,mapNcol,mapNsec)] += (float) (at(i)._mass * fracX     * fracY     * fracZ    );
        }
    }

    // filter map with a Gaussian kernel
    float                   sigma     = (float) (reso / (2.0 * spacing * std::sqrt(3.0)));
    size_t                  kernelExt = (size_t) std::ceil(3.0f * sigma);
    cFilterSpatial<float>   objFilter;

    objFilter.kernelSpatialGaussian(sigma, 2*kernelExt+1, 3);
    objFilter.filterSpatial(map, XCORR_RES_FULL);

    // update header information
    if (header != NULL) {
        cVector3<ptrdiff_t>     mapOrigin((ptrdiff_t) (mapMinX/spacing) - kernelExt - mapMargin,
                                          (ptrdiff_t) (mapMinY/spacing) - kernelExt - mapMargin,
                                          (ptrdiff_t) (mapMinZ/spacing) - kernelExt - mapMargin);

        header->simulate(map.getSize(), mapOrigin, (float) spacing);
    }
}

double cPDB::opComputeRMSD(const cPDB& other, bool CAonly) const
{
    const std::string   funcName("double cPDB::opComputeRMSD(const cPDB& other) const");

    requireNonEmpty(funcName);
    require(size() == other.size(),
            funcName + ": PDBs have different sizes");

    double      sqrSum = 0;

    #pragma omp parallel for reduction(+:sqrSum)
    for (size_t i = 0; i < size(); i++) {
        if (CAonly) {
            if (at(i)._atomName == " CA ") {
                sqrSum += (at(i)._coord - other.at(i)._coord).getLength();
            }
            else {
                sqrSum += 0;
            }
        }
        else {
            sqrSum += (at(i)._coord - other.at(i)._coord).getLength();
        }
    }

    return std::sqrt(sqrSum / (double) size());
}

double cPDB::opEstimateVolume(void) const
{
    const std::string   funcName("double cPDB::opEstimateVolume(void) const");

    requireNonEmpty(funcName);
}

void cPDB::read(const std::string& fileName)
{
    const std::string   funcName("void cPDB::read(const std::string& fileName)");

    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".pdb",
            funcName + ": invalid file name = " + fileName);

    std::ifstream       fileHandle;
    std::string         lineStr;
    sAtom               atom;
    std::string         recordName;

    clear();

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    while (!fileHandle.eof()) {

        std::getline(fileHandle, lineStr);

        if (lineStr != "") { recordName = lineStr.substr(0,6); }
        else               { recordName = "";                  }

        if (recordName == "MODEL ") {
            MESSAGE(funcName + ": multiple models is not supported");
        }
        //require(recordName != "MODEL ", funcName + ": multiple models is not supported");

        if (recordName == "ATOM  " || recordName == "HETATM") {

            // extract atom record
            atom._recordName    = recordName;
            atom._serialNum     = std::atoi(lineStr.substr( 6,5).c_str());
            atom._atomName      =           lineStr.substr(12,4);
            atom._location      =           lineStr.substr(16,1);
            atom._residueName   =           lineStr.substr(17,3);
            atom._chainName     =           lineStr.substr(21,1);
            atom._sequenceNum   = std::atoi(lineStr.substr(22,4).c_str());
            atom._insertionCode =           lineStr.substr(26,1);
            atom._coord[0]      = std::atof(lineStr.substr(30,8).c_str());
            atom._coord[1]      = std::atof(lineStr.substr(38,8).c_str());
            atom._coord[2]      = std::atof(lineStr.substr(46,8).c_str());
            atom._occupancy     = std::atof(lineStr.substr(54,6).c_str());
            atom._temperature   = std::atof(lineStr.substr(60,6).c_str());
            atom._segmentID     =           lineStr.substr(72,4);
            atom._elementSym    =           lineStr.substr(76,2);
            atom._charge        =           lineStr.substr(78,2);

            // assign mass to each atom
            atom._mass = 0;
            for (size_t i = 0; i < 30; i++) {
                if (atom._atomName == atom_name[i]) {
                    atom._mass = atom_mass[i];
                }
            }

            if (atom._mass == 0) {
                for (size_t j = 0; j < 6; j++) {
                    if (atom._atomName[1] == atom_name[atom_frequent[j]][1]) {
                        atom._mass = atom_mass[atom_frequent[j]];
                    }
                }
            }

            if (atom._mass == 0) {
                WARNING(funcName, "cannot assign mass to atom");
                std::cout << atom << std::endl;
            }

            //require(atom._mass > 0, funcName + ": cannot assign mass to atom");

            //std::cout << lineStr << std::endl;
            //std::cout << atom    << std::endl;

            push_back(atom);
        }
    }

    fileHandle.close();
}

void cPDB::write(const std::string& fileName) const
{
    const std::string   funcName("void cPDB::write(const std::string& fileName) const");

    requireNonEmpty(funcName);
    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".pdb",
            funcName + ": invalid file name = " + fileName);

    std::ofstream       fileHandle;
    sAtom               atom;

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    for (size_t i = 0; i < size(); i++) {

        atom = at(i);

        fileHandle << std::fixed
              << atom._recordName
              << std::setw(5) << atom._serialNum
              << " "
              << atom._atomName
              << atom._location
              << atom._residueName
              << " "
              << atom._chainName
              << std::setw(4) << atom._sequenceNum
              << atom._insertionCode
              << "   "
              << std::setw(8) << std::setprecision(3) << atom._coord[0]
              << std::setw(8) << std::setprecision(3) << atom._coord[1]
              << std::setw(8) << std::setprecision(3) << atom._coord[2]
              << std::setw(6) << std::setprecision(2) << atom._occupancy
              << std::setw(6) << std::setprecision(2) << atom._temperature
              << "      "
              << atom._segmentID
              << atom._elementSym
              << atom._charge
              << std::endl;
    }

    fileHandle.close();
}

void cPDB::writeModels(const std::string& fileName,
                         const std::vector<cPDB>& others) const
{
    const std::string   funcName("void cPDB::writeModels(const std::string& fileName, const std::vector<cPDB>& others) const");

    requireEmpty(funcName);
    require(others.size() > 1, funcName + ": no need to merge");
    require(fileName != "" && cFileSystem(fileName).getFileExt() == ".pdb",
            funcName + ": invalid file name = " + fileName);

    std::ofstream       fileHandle;
    sAtom               atom;

    fileHandle.open(fileName.c_str());
    assure(fileHandle, fileName + " in " + funcName);

    for (size_t i = 0; i < others.size(); i++) {

        require(others[i].size() > 0, funcName + ": PDB is empty");

        fileHandle << "MODEL " << "    " << std::setw(4) << num2str(i+1) << std::endl;

        for (size_t j = 0; j < others[i].size(); j++) {

            atom = others[i].at(j);

            fileHandle << std::fixed
                  << atom._recordName
                  << std::setw(5) << atom._serialNum
                  << " "
                  << atom._atomName
                  << atom._location
                  << atom._residueName
                  << " "
                  << atom._chainName
                  << std::setw(4) << atom._sequenceNum
                  << atom._insertionCode
                  << "   "
                  << std::setw(8) << std::setprecision(3) << atom._coord[0]
                  << std::setw(8) << std::setprecision(3) << atom._coord[1]
                  << std::setw(8) << std::setprecision(3) << atom._coord[2]
                  << std::setw(6) << std::setprecision(2) << atom._occupancy
                  << std::setw(6) << std::setprecision(2) << atom._temperature
                  << "      "
                  << atom._segmentID
                  << atom._elementSym
                  << atom._charge
                  << std::endl;
        }

        fileHandle << "ENDMDL" << std::endl;
    }

    fileHandle.close();
}

} // namespace gem
