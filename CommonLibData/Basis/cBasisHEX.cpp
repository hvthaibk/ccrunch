/***********************************************************************
 *  File:       cBasisHEX.cpp
 *
 *  Purpose:    Implementation of a basis-related class
 *
 *  Author:     Thai V. Hoang
 *
 *  Contact:    hvthaibk@gmail.com
 *
 *  Copyright (C) 2012 Thai V. Hoang, INRIA
 **********************************************************************/

#include "cBasisHEX.hpp"

#ifdef __GEM_USE_HEX__
#include "HEX/hex_math.h"
#include "HEX/hex_dock.h"
#include "HEX/kids.h"
#endif

namespace gem {

#ifdef __GEM_USE_HEX__

template <typename T>
void cBasisHEX<T>::hexComputeCoeff(const cVector3<double>& origin,
                                     const cVector3<double>& spacing,
                                     unsigned int nmax, double scale,
                                     cData3<double>& coff) const
{
    const std::string   funcName("void cBasisHEX<T>::hexComputeCoeff("
                                    "const cVector3<double>& origin, "
                                    "const cVector3<double>& spacing, "
                                    "unsigned int nmax, double scale, "
                                    "cData3<double>& coff) const)");

    require(nmax > 0, funcName + ": nmax is not positive");
    require(cArray3<T>::getDimension() == 3, funcName + ": input data is not 3D");

    int             n_max, l_max, nlm_dim, nl_dim, lm_dim;
    Point3D         pt;
    Polar           pp;
    cData3<double>  rfn, yfn, rnl_norm, ylm_norm;

    n_max = nmax;
    l_max = n_max - 1;

    nlm_dim = NLM_DIM(n_max);
    nl_dim  = NL_DIM(n_max);
    lm_dim  = YLM_DIM(l_max);

    // allocate memory
    coff.memReAllocZero(cSize3(nlm_dim,1,1));

    rfn.memAlloc(cSize3(nl_dim,1,1));
    yfn.memAlloc(cSize3(lm_dim,1,1));

    rnl_norm.memAlloc(cSize3(nl_dim,1,1));
    ylm_norm.memAlloc(cSize3(lm_dim,1,1));

    // prepare HEX
    hnl_factor(20*std::pow(scale,2));
    hnl_all_norm(n_max, rnl_norm.getAddrData());
    ylm_all_norm(l_max, ylm_norm.getAddrData());

    // avoid using rnl_norm *= spacing.getProduct() for efficiency
    for (int nl = 0; nl < nl_dim; nl++) {
        rnl_norm[nl] *= spacing.getProduct();
    }

    // main loop: browsing all the voxels of the map
    for (size_t iRow = 0, mapIdx = 0; iRow < cArray3<T>::getNrow(); iRow++) {
        for (size_t iCol = 0; iCol < cArray3<T>::getNcol(); iCol++) {
            for (size_t iSec = 0; iSec < cArray3<T>::getNsec(); iSec++, mapIdx++) {
                // Cartesian to polar
                pt.x = ((double) iRow + 0.5 - origin[0]) * spacing[0];
                pt.y = ((double) iCol + 0.5 - origin[1]) * spacing[1];
                pt.z = ((double) iSec + 0.5 - origin[2]) * spacing[2];

                pp = pp_pt(pt);

                // radial functions
                hnl_all_bar(n_max, pp.r, rfn.getAddrData());

                for (int nl = 0; nl < nl_dim; nl++) {
                    rfn[nl] *= rnl_norm[nl];    // avoid using rfn *= rnl_norm for efficiency
                }

                // spherical harmonics
                ylm_all_bar(l_max, pp.theta, pp.phi, yfn.getAddrData());

                for (int lm = 0; lm < lm_dim; lm++) {
                    yfn[lm] *= ylm_norm[lm];    // avoid using yfn *= ylm_norm for efficiency
                }

                // coefficients
                for (int n = 1, nl = 0, nc = 0; n <= n_max; n++) {
                    for (int l = 0, lm = 0; l < n; l++, nl++) {
                        for (int m = -l; m <= l; m++, nc++, lm++) {
                            coff[nc] += cArray3<T>::_data[mapIdx] * rfn[nl] * yfn[lm];
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
void cBasisHEX<T>::hexReconData(const cVector3<double>& origin,
                                  const cVector3<double>& spacing,
                                  unsigned int nmax, double scale,
                                  const cData3<double>& coff)
{
    const std::string   funcName("void cBasisHEX<T>::hexReconData("
                                    "const cVector3<double>& origin, "
                                    "const cVector3<double>& spacing, "
                                    "unsigned int nmax, double scale, "
                                    "const cData3<double>& coff)");

    require(nmax > 0, funcName + ": nmax is not positive");
    require(cArray3<T>::getDimension() == 3, funcName + ": output data is not 3D");

    int             n_max, l_max, nlm_dim, nl_dim, lm_dim;
    Point3D         pt;
    Polar           pp;
    cData3<double>  rfn, yfn, rnl_norm, ylm_norm;

    n_max = nmax;
    l_max = n_max - 1;

    nlm_dim = NLM_DIM(n_max);
    nl_dim  = NL_DIM(n_max);
    lm_dim  = YLM_DIM(l_max);

    require(coff.getSize() >= cSize3(nlm_dim,1,1), funcName + ": coefficient data seems incorrect");

    // allocate memory
    cArray3<T>::memSetZero();

    rfn.memAlloc(cSize3(nl_dim,1,1));
    yfn.memAlloc(cSize3(lm_dim,1,1));

    rnl_norm.memAlloc(cSize3(nl_dim,1,1));
    ylm_norm.memAlloc(cSize3(lm_dim,1,1));

    // prepare HEX
    hnl_factor(20*std::pow(scale,2));
    hnl_all_norm(n_max, rnl_norm.getAddrData());
    ylm_all_norm(l_max, ylm_norm.getAddrData());

    // avoid using rnl_norm *= spacing.getProduct() for efficiency
    for (int nl = 0; nl < nl_dim; nl++) {
        rnl_norm[nl] *= spacing.getProduct();
    }

    // main loop: browsing all the voxels of the map
    for (size_t iRow = 0, mapIdx = 0; iRow < cArray3<T>::getNrow(); iRow++) {
        for (size_t iCol = 0; iCol < cArray3<T>::getNcol(); iCol++) {
            for (size_t iSec = 0; iSec < cArray3<T>::getNsec(); iSec++, mapIdx++) {
                // Cartesian to polar
                pt.x = ((double) iRow + 0.5 - origin[0]) * spacing[0];
                pt.y = ((double) iCol + 0.5 - origin[1]) * spacing[1];
                pt.z = ((double) iSec + 0.5 - origin[2]) * spacing[2];

                pp = pp_pt(pt);

                // radial functions
                hnl_all_bar(n_max, pp.r, rfn.getAddrData());

                for (int nl = 0; nl < nl_dim; nl++) {
                    rfn[nl] *= rnl_norm[nl];    // avoid using rfn *= rnl_norm for efficiency
                }

                // spherical harmonics
                ylm_all_bar(l_max, pp.theta, pp.phi, yfn.getAddrData());

                for (int lm = 0; lm < lm_dim; lm++) {
                    yfn[lm] *= ylm_norm[lm];    // avoid using yfn *= ylm_norm for efficiency
                }

                // coefficients
                for (int n = 1, nl = 0, nc = 0; n <= n_max; n++) {
                    for (int l = 0, lm = 0; l < n; l++, nl++) {
                        for (int m = -l; m <= l; m++, nc++, lm++) {
                            cArray3<T>::_data[mapIdx] +=  coff[nc] * rfn[nl] * yfn[lm];
                        }
                    }
                }
            }
        }
    }
}

template <typename T>
void cBasisHEX<T>::hexSetupOption(DockSpec& ds)
{
    hex_setDockSpec(&ds);               /* pick up the standard usual defaults */

    ds.fft_type             = 3;        /* 3=CPU 3D FFT on CPU  1=GPU 1D FFT on GPU*/
    ds.fft_device           = 0;        /* 0=CPU, 1=GPU, 2=BOTH (if possible) */

    ds.n_min                = 0;        /* no initial steric scan */
    ds.n_max                = 10;       /* main correlation to order k_max */

    ds.r12_range            = 6;        /* search distance range in Angstrom */
    ds.r12_step             = 0.8;      /* in steps of 0.75 Angstrom ... */
    ds.r12_substeps         = 0;        /* not used ... */
    ds.r12_guess            = 0.0;      /* both molecules start at the origin */
    ds.grid_size            = 0.6;      /* size of side of sampling grid */
    ds.integration_radius   = 40.0;     /* grid integration radius in Angstrom */

    ds.receptor_range_angle = PI;       /* receptor range (beta,gamma) in radians */
    ds.ligand_range_angle   = PI;       /* ligand range (beta,gamma) in radians */
    ds.alpha_range_angle    = (2*PI);   /* twist angle (alpha) in radians */

    ds.problem_type         = 1;        /* 1=similarity, 2=docking,3=multidock */
    ds.potential_type       = 0;        /* 0=tau, 1=sigma, 2=tau+sigma */

// FIX ME - currently, SIG is no loaded, and is never used?

    ds.refine               = 0;        /* no docking refinement here */

    ds.max_solutions        = 1000;     /* should not normally be changed */

    ds.receptor_stepsize    = 20;       /* in degrees */
    ds.ligand_stepsize      = 10;       /* in degrees */
    ds.alpha_stepsize       = 5.6;      /* in degrees */

    ds.receptor_bandlimit   = 12;       /* calculated from the stepsizes */
    ds.ligand_bandlimit     = 18;
    ds.alpha_bandlimit      = 32;

    ds.receptor_samples     = 162;
    ds.ligand_samples       = 492;
    ds.alpha_samples        = 64;
}

template <typename T>
void cBasisHEX<T>::hexMatcher(DockSpec& ds, const cBasisHEX<T>& other, Docking& docking)
{
   int          n, ns;
   double       r12, s12, score;
   Euler        eq, et;
   Matrix3D     rt;                 // target rotation
   Vector3D     vt;                 // target translation
   Molecule     *mol1, *mol2;

   mol1 = (Molecule *) kgetB_0(sizeof(Molecule));
   mol2 = (Molecule *) kgetB_0(sizeof(Molecule));

   mol1->open                = 1;
   mol2->open                = 1;

   mol1->model_type          = 0;
   mol2->model_type          = 0;

   mol1->model               = (Model *) kgetB_0(sizeof(Model));
   mol1->model->id           = 1;
   //mol1->model->tau          = tau1;
   //mol1->model->sig          = sig1;

   strcpy(mol1->model->label, "A");

   mol2->model               = (Model *) kgetB_0(sizeof(Model));
   mol2->model->id           = 1;
   //mol2->model->tau          = tau2;
   //mol2->model->sig          = sig2;

   strcpy(mol2->model->label, "B");

   hex_matcher(&ds, &docking, mol1, mol2);

// for consistence with docking, Hex outputs matching scores as negative numbers
// in the range (0 to -1000), so it is simplest to convert them back to the
// NB. These are raw unclustered scores. There will be many near-duplicates.

   for (n=0; n<docking.num_solutions; n++) {

      docking.solution[n].e_total /= - 1000;
      docking.solution[n].e_shape /= - 1000;
   }

   //hex_msg("Scores for fixed query = %s + moving target = %s\n\n", name1, name2);

   hex_msg(" Soln  Score    Target Euler Rotation  PDB Translation  PDB (SPF) Shift\n");
   hex_msg("-----  ------   ---------------------  ---------------  ---------------\n");

   //ns = min_(s_max, docking.num_solutions);

   for (n=0; n<ns; n++) {

// z-translations were calculated in SPF units, so just make a note for output

      s12 = docking.solution[n].r12;

// convert z-translations back to PDB units for consumption by other functions

      //docking.solution[n].r12 /= kpax_potScale(pot);

      score = docking.solution[n].e_total;

      eq    = docking.solution[n].reuler;   // mol1 = "receptor" = stationary QUERY
      et    = docking.solution[n].leuler;   // mol2 = "ligand"   = moving TARGET
      r12   = docking.solution[n].r12;

      Matrix3D tt = tf_tftf(tf_tftf(       // put all motion on the target
                               tf_euler(et),
                               tf_v(v_make(0.0, 0.0, r12))),
                               tf_transpose(tf_euler(eq)));

      tr_tf(tt, &vt, &rt);         // matrix to rotation plus translation

      Euler er = euler_tf(rt);     // rotation to Euler

      hex_msg("%5d %7.4f   %6.1f %6.1f %6.1f   %5.1f %5.1f %5.1f  %5.1f (%.1f)\n",
              n+1, score, rad2deg(er.alpha), rad2deg(er.beta), rad2deg(er.gamma),
              vt.x, vt.y, vt.z, r12, s12);
   }

   hex_msg("----------------------------------------------------------------------\n");

   kfree(mol1->model);
   kfree(mol2->model);

   kfree(mol1);
   kfree(mol2);
}

#endif

// instantiation
template class cBasisHEX<float >;
template class cBasisHEX<double>;

} // namespace gem
