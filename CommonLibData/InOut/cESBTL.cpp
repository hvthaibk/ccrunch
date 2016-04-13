/***********************************************************************
 *  File:       cESBTL.cpp
 *
 *  Purpose:    Implementation of an ESBTL class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cESBTL.hpp"

#include "macro.hpp"
#include "transform.hpp"
#include "cMatrix33.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "ESBTL/default.h"
#include <ESBTL/atom_classifier.h>
#include <ESBTL/atom_selectors.h>
#include <ESBTL/weighted_atom_iterator.h>
#include <ESBTL/selected_atom_iterator.h>
#pragma GCC diagnostic pop

#include <iostream>
#include <fstream>

namespace ESBTL {
    std::ostream& operator<<(std::ostream& os, const ESBTL::Default_system::Atom& atom) {
        os << ESBTL::PDB::get_atom_pdb_format(atom);
        return os;
    }
}

namespace gem {

typedef ESBTL::Accept_all_occupancy_policy<ESBTL::PDB::Line_format<> >                              Occupancy_policy_all;
typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,ESBTL::Default_system::Atom> >       Atom_classifier;
typedef ESBTL::Selected_atom_iterator<ESBTL::Default_system::Model,ESBTL::Select_by_chainids,true>  Restrict_iterator;

void cESBTL::checkExtension(const std::string& fileName) const
{
    const std::string   funcName("void cESBTL::checkExtension(const std::string& fileName) const");

    if (fileName == "" || (fileName.find(".pdb") > fileName.length() &&
                           fileName.find(".PDB") > fileName.length())) {
        ERROR(funcName, "unsupported expansion (accept *.pdb, *.PDB)");
    }
}

void cESBTL::pdbPrintAtom(const std::string& fileIn, const std::string& atomName)
{
    checkExtension(fileIn);

    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    // iterate over all models
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;
        //std::cout << "Model " << model.model_number() << std::endl;

        for (ESBTL::Default_system::Model::Chains_const_iterator iterChain  = model.chains_begin();
                                                                 iterChain != model.chains_end();
                                                                 iterChain++) {

            const ESBTL::Default_system::Chain&     chain = *iterChain;
            //std::cout << "Model "   << model.model_number()
            //          << ", chain " << chain.chain_identifier() << std::endl;

            for (ESBTL::Default_system::Chain::Residues_const_iterator iterResidue  = chain.residues_begin();
                                                                       iterResidue != chain.residues_end();
                                                                       iterResidue++) {

                const ESBTL::Default_system::Residue&     residue = *iterResidue;
                //std::cout << "Model "     << model.model_number()
                //          << ", chain "   << chain.chain_identifier()
                //          << ", residue " << residue.residue_name() << std::endl;

                for (ESBTL::Default_system::Residue::Atoms_const_iterator iterAtom  = residue.atoms_begin();
                                                                          iterAtom != residue.atoms_end();
                                                                          iterAtom++) {
                    if (atomName == "" || iterAtom->atom_name() == atomName) {
                        std::cout << "Model "     << model.model_number()
                                  << ", chain "   << chain.chain_identifier()
                                  << ", residue " << residue.residue_name()
                                  << ", atom "    << iterAtom->atom_name();
                        std::cout << std::setprecision(4) << std::fixed
                                  << std::setw(7)  << iterAtom->atom_serial_number()
                                  << std::setw(5)  << iterAtom->charge()
                                  << std::setw(4)  << iterAtom->element()
                                  << std::setw(11) << iterAtom->x()
                                  << std::setw(11) << iterAtom->y()
                                  << std::setw(11) << iterAtom->z();
                        std::cout << std::endl;
                    }
                }
            }
        }
    }
}

void cESBTL::pdbReadAtom (const std::string& fileIn, cSetPoint3<double>& atomList, const std::string& atomName)
{
    checkExtension(fileIn);

    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            if ((atomName != "" && iterAtom->atom_name() == atomName) || atomName == "") {
                atomList.push_back(cVector3<double>(iterAtom->x(),iterAtom->y(),iterAtom->z()));
            }
        }
    }
}

void cESBTL::pdbExtractChain(const std::string& fileIn, const std::string& fileOut, const std::string& chainIDs)
{
    checkExtension(fileIn);
    checkExtension(fileOut);

    // input
    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    // output
    std::ofstream               output(fileOut.c_str());
    assure(output);

    // selection
    ESBTL::Select_by_chainids   pdbSelectorChainIDs(chainIDs);

    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (Restrict_iterator iterAtom  = ESBTL::make_selected_atom_iterator(model.atoms_begin(),pdbSelectorChainIDs);
                               iterAtom != ESBTL::make_selected_atom_iterator<ESBTL::Select_by_chainids>(model.atoms_end());
                               iterAtom++) {
            output << ESBTL::PDB::get_atom_pdb_format(*iterAtom) << "\n";
        }
        output.close();
    }
}

void cESBTL::pdbRotate(const std::string& fileIn, const std::string& fileOut, const cVector3<double>& angle)
{
    checkExtension(fileIn);
    checkExtension(fileOut);

    // input
    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    // output
    std::ofstream               output(fileOut.c_str());
    assure(output);

    // rotate
    cVector3<double>            pointRot;
    cMatrix33<double>           matRot;

    transform_rotmat(angle[0], angle[1], angle[2],
                     matRot[0][0], matRot[0][1], matRot[0][2],
                     matRot[1][0], matRot[1][1], matRot[1][2],
                     matRot[2][0], matRot[2][1], matRot[2][2],
                     ROT3D_RIGHT_ZYZ);

    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            ESBTL::Default_system::Atom     atom(*iterAtom);

            pointRot = matRot * cVector3<double>(atom.x(), atom.y(), atom.z());

            atom.x() = pointRot[0];
            atom.y() = pointRot[1];
            atom.z() = pointRot[2];

            output << ESBTL::PDB::get_atom_pdb_format(atom) << "\n";
        }
    }
}

void cESBTL::pdbTranslate(const std::string& fileIn, const std::string& fileOut, const cVector3<double>& offset)
{
    checkExtension(fileIn);
    checkExtension(fileOut);

    // input
    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    // output
    std::ofstream               output(fileOut.c_str());
    assure(output);

    // translate
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            ESBTL::Default_system::Atom     atom(*iterAtom);

            atom.x() += offset[0];
            atom.y() += offset[1];
            atom.z() += offset[2];

            output << ESBTL::PDB::get_atom_pdb_format(atom) << "\n";
        }
    }
}

void cESBTL::pdbTransRotTrans(const std::string& fileIn, const std::string& fileOut, const cVector3<double>& offset1, const cVector3<double>& angle, const cVector3<double>& offset2)
{
    checkExtension(fileIn);
    checkExtension(fileOut);

    // input
    ESBTL::PDB_line_selector                                pdbSelector;
    std::vector<ESBTL::Default_system>                      pdbSystems;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder(pdbSystems, pdbSelector.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn, pdbSelector, pdbBuilder, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems.empty(),     "No PDB system");
    require( pdbSystems.size() == 1, "More than one PDB system: pdbSystems.size() = " + num2str(pdbSystems.size()));
    require(!pdbSystems[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems[0].number_of_models() = "
                                                 + num2str(pdbSystems[0].number_of_models()));

    // output
    std::ofstream               output(fileOut.c_str());
    assure(output);

    // translate - rotate - translate
    cVector3<double>            pointRot;
    cMatrix33<double>           matRot;

    transform_rotmat(angle[0], angle[1], angle[2],
                     matRot[0][0], matRot[0][1], matRot[0][2],
                     matRot[1][0], matRot[1][1], matRot[1][2],
                     matRot[2][0], matRot[2][1], matRot[2][2],
                     ROT3D_RIGHT_ZYZ);

    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems[0].models_begin();
                                                iterModel != pdbSystems[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            ESBTL::Default_system::Atom     atom(*iterAtom);

            pointRot = matRot * cVector3<double>(atom.x() + offset1[0],
                                                 atom.y() + offset1[1],
                                                 atom.z() + offset1[2]);

            atom.x() = pointRot[0] + offset2[0];
            atom.y() = pointRot[1] + offset2[1];
            atom.z() = pointRot[2] + offset2[2];

            output << ESBTL::PDB::get_atom_pdb_format(atom) << "\n";
        }
    }
}

double cESBTL::pdbComputeRMSD(const std::string& fileIn1, const std::string& fileIn2, const std::string& atomName)
{
    checkExtension(fileIn1);
    checkExtension(fileIn2);

    // input 1
    ESBTL::PDB_line_selector                                pdbSelector1;
    std::vector<ESBTL::Default_system>                      pdbSystems1;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder1(pdbSystems1, pdbSelector1.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn1, pdbSelector1, pdbBuilder1, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems1.empty(),     "No PDB system");
    require( pdbSystems1.size() == 1, "More than one PDB system: pdbSystems1.size() = " + num2str(pdbSystems1.size()));
    require(!pdbSystems1[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems1[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems1[0].number_of_models() = "
                                                 + num2str(pdbSystems1[0].number_of_models()));

    cSetPoint3<double>      atomList1;
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems1[0].models_begin();
                                                iterModel != pdbSystems1[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            if ((atomName != "" && iterAtom->atom_name() == atomName) || atomName == "") {
                atomList1.push_back(cVector3<double>(iterAtom->x(),iterAtom->y(),iterAtom->z()));
            }
        }
    }

    // input 2
    ESBTL::PDB_line_selector                                pdbSelector2;
    std::vector<ESBTL::Default_system>                      pdbSystems2;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder2(pdbSystems2, pdbSelector2.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn2, pdbSelector2, pdbBuilder2, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems2.empty(),     "No PDB system");
    require( pdbSystems2.size() == 1, "More than one PDB system: pdbSystems2.size() = " + num2str(pdbSystems2.size()));
    require(!pdbSystems2[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems2[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems2[0].number_of_models() = "
                                                 + num2str(pdbSystems2[0].number_of_models()));

    cSetPoint3<double>      atomList2;
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems2[0].models_begin();
                                                iterModel != pdbSystems2[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            if ((atomName != "" && iterAtom->atom_name() == atomName) || atomName == "") {
                atomList2.push_back(cVector3<double>(iterAtom->x(),iterAtom->y(),iterAtom->z()));
            }
        }
    }

    return atomList1.computeRMSD(atomList2);
}

void cESBTL::pdbComputeSD(const std::string& fileIn1, const std::string& fileIn2, double& sqrDist, size_t& numAtom, const std::string& atomName)
{
    checkExtension(fileIn1);
    checkExtension(fileIn2);

    // input 1
    ESBTL::PDB_line_selector                                pdbSelector1;
    std::vector<ESBTL::Default_system>                      pdbSystems1;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder1(pdbSystems1, pdbSelector1.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn1, pdbSelector1, pdbBuilder1, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems1.empty(),     "No PDB system");
    require( pdbSystems1.size() == 1, "More than one PDB system: pdbSystems1.size() = " + num2str(pdbSystems1.size()));
    require(!pdbSystems1[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems1[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems1[0].number_of_models() = "
                                                 + num2str(pdbSystems1[0].number_of_models()));

    cSetPoint3<double>      atomList1;
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems1[0].models_begin();
                                                iterModel != pdbSystems1[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            if ((atomName != "" && iterAtom->atom_name() == atomName) || atomName == "") {
                atomList1.push_back(cVector3<double>(iterAtom->x(),iterAtom->y(),iterAtom->z()));
            }
        }
    }

    // input 2
    ESBTL::PDB_line_selector                                pdbSelector2;
    std::vector<ESBTL::Default_system>                      pdbSystems2;
    ESBTL::All_atom_system_builder<ESBTL::Default_system>   pdbBuilder2(pdbSystems2, pdbSelector2.max_nb_systems());

    require(ESBTL::read_a_pdb_file(fileIn2, pdbSelector2, pdbBuilder2, Occupancy_policy_all()), "cannot read PDB file");
    require(!pdbSystems2.empty(),     "No PDB system");
    require( pdbSystems2.size() == 1, "More than one PDB system: pdbSystems2.size() = " + num2str(pdbSystems2.size()));
    require(!pdbSystems2[0].has_no_model(),          "No model in 1st PDB system");
    require( pdbSystems2[0].number_of_models() == 1, "More than one model in 1st PDB system: pdbSystems2[0].number_of_models() = "
                                                 + num2str(pdbSystems2[0].number_of_models()));

    cSetPoint3<double>      atomList2;
    for (ESBTL::Default_system::Models_iterator iterModel  = pdbSystems2[0].models_begin();
                                                iterModel != pdbSystems2[0].models_end();
                                                iterModel++) {

        const ESBTL::Default_system::Model&     model = *iterModel;

        for (ESBTL::Default_system::Model::Atoms_const_iterator iterAtom  = model.atoms_begin();
                                                                iterAtom != model.atoms_end();
                                                                iterAtom++) {
            if ((atomName != "" && iterAtom->atom_name() == atomName) || atomName == "") {
                atomList2.push_back(cVector3<double>(iterAtom->x(),iterAtom->y(),iterAtom->z()));
            }
        }
    }

    numAtom = atomList1.size();
    sqrDist = atomList1.computeSD(atomList2);
}

} // namespace gem
