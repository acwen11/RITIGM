//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <algorithm>
#include <cassert>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_printer.hh"
#include "thc_M1_macro.h"

#include "utils.hh"

using namespace thc;
using namespace std;

extern "C" void THC_M1_CalcOpacity(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    // Opacities are constant throught the timestep
    if (*TimeIntegratorStage != 2) {
        return;
    }

    if (verbose) {
        CCTK_INFO("THC_M1_CalcOpacity");
    }

    CCTK_REAL const dt = CCTK_DELTA_TIME;

    // Setup Printer
    thc::Printer::start(
            "[INFO|THC|THC_M1_CalcOpacity]: ",
            "[WARN|THC|THC_M1_CalcOpacity]: ",
            "[ERR|THC|THC_M1_CalcOpacity]: ",
            m1_max_num_msg, m1_max_num_msg);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_calc_opacity,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if (thc_m1_mask[ijk]) {
                for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                    abs_0[i4D] = 0.0;
                    abs_1[i4D] = 0.0;
                    eta_0[i4D] = 0.0;
                    eta_1[i4D] = 0.0;
                    scat_1[i4D] = 0.0;
                }
                continue;
            }

            assert(nspecies == 3);
            assert(ngroups == 1);

            // Get the transport opacity (absorption + scattering)
            CCTK_REAL kappa_0_loc[3], kappa_1_loc[3];
            int ierr = NeutrinoOpacity(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &kappa_0_loc[0], &kappa_0_loc[1], &kappa_0_loc[2],
                    &kappa_1_loc[0], &kappa_1_loc[1], &kappa_1_loc[2]);
            assert(!ierr);
            assert(isfinite(kappa_0_loc[0]));
            assert(isfinite(kappa_0_loc[1]));
            assert(isfinite(kappa_0_loc[2]));
            assert(isfinite(kappa_1_loc[0]));
            assert(isfinite(kappa_1_loc[1]));
            assert(isfinite(kappa_1_loc[2]));

            // Get the absorption opacity
            CCTK_REAL abs_0_loc[3], abs_1_loc[3];
            ierr = NeutrinoAbsorptionRate(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &abs_0_loc[0], &abs_0_loc[1], &abs_0_loc[2],
                    &abs_1_loc[0], &abs_1_loc[1], &abs_1_loc[2]);
            assert(!ierr);
            assert(isfinite(abs_0_loc[0]));
            assert(isfinite(abs_0_loc[1]));
            assert(isfinite(abs_0_loc[2]));
            assert(isfinite(abs_1_loc[0]));
            assert(isfinite(abs_1_loc[1]));
            assert(isfinite(abs_1_loc[2]));

            // An effective optical depth used to decide whether to compute
            // the black body function for neutrinos assuming neutrino trapping
            // or at a fixed temperature and Ye
            CCTK_REAL const tau = min(
                    sqrt(abs_1_loc[0]*kappa_1_loc[0]),
                    sqrt(abs_1_loc[1]*kappa_1_loc[1]))*dt;

            // Compute the neutrino emission rates
            CCTK_REAL eta_0_loc[3], eta_1_loc[3];
            ierr = NeutrinoEmission(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &eta_0_loc[0], &eta_0_loc[1], &eta_0_loc[2],
                    &eta_1_loc[0], &eta_1_loc[1], &eta_1_loc[2]);
            assert(!ierr);
            assert(isfinite(eta_0_loc[0]));
            assert(isfinite(eta_0_loc[1]));
            assert(isfinite(eta_0_loc[2]));
            assert(isfinite(eta_1_loc[0]));
            assert(isfinite(eta_1_loc[1]));
            assert(isfinite(eta_1_loc[2]));

            // Compute the neutrino black body functions assuming trapped neutrinos
            CCTK_REAL nudens_0_trap[3], nudens_1_trap[3];
            if (opacity_tau_trap >= 0 && tau > opacity_tau_trap) {
                CCTK_REAL temperature_trap, Y_e_trap;
                // Compute local neutrino densities (undensitized)
                CCTK_REAL const nudens_0[3] = {
                    rnnu[CCTK_VectGFIndex3D(cctkGH, i, j, k, 0)]/volform[ijk],
                    rnnu[CCTK_VectGFIndex3D(cctkGH, i, j, k, 1)]/volform[ijk],
                    rnnu[CCTK_VectGFIndex3D(cctkGH, i, j, k, 2)]/volform[ijk],
                };
                CCTK_REAL const nudens_1[3] = {
                    rJ[CCTK_VectGFIndex3D(cctkGH, i, j, k, 0)]/volform[ijk],
                    rJ[CCTK_VectGFIndex3D(cctkGH, i, j, k, 1)]/volform[ijk],
                    rJ[CCTK_VectGFIndex3D(cctkGH, i, j, k, 2)]/volform[ijk],
                };
                ierr = WeakEquilibrium(
                        rho[ijk], temperature[ijk], Y_e[ijk],
                        nudens_0[0], nudens_0[1], nudens_0[2],
                        nudens_1[0], nudens_1[1], nudens_1[2],
                        &temperature_trap, &Y_e_trap,
                        &nudens_0_trap[0], &nudens_0_trap[1], &nudens_0_trap[2],
                        &nudens_1_trap[0], &nudens_1_trap[1], &nudens_1_trap[2]);
                if (ierr) {
                    // Try to recompute the weak equilibrium using neglecting
                    // current neutrino data
                    ierr = WeakEquilibrium(
                            rho[ijk], temperature[ijk], Y_e[ijk],
                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            &temperature_trap, &Y_e_trap,
                            &nudens_0_trap[0], &nudens_0_trap[1], &nudens_0_trap[2],
                            &nudens_1_trap[0], &nudens_1_trap[1], &nudens_1_trap[2]);
                    if (ierr) {
                        ostringstream ss;
                        ss << "Could not find the weak equilibrium!" << endl;
                        ss << "Reflevel = " << ilogb(cctkGH->cctk_levfac[0]) << endl;
                        ss << "Iteration = " << cctk_iteration << endl;
                        ss << "(i, j, k) = (" << i << ", " << j << ", " << k << ")\n";
                        ss << "(x, y, z) = (" << x[ijk] << ", " << y[ijk] << ", "
                                              << z[ijk] << ")\n";
                        ss << "rho = " << rho[ijk] << endl;
                        ss << "temperature = " << temperature[ijk] << endl;
                        ss << "Y_e = " << Y_e[ijk] << endl;
                        ss << "alp = " << alp[ijk] << endl;
                        ss << "nudens_0 = " << nudens_0[0] << " " << nudens_0[1]
                                            << " " << nudens_0[2] << endl;
                        ss << "nudens_1 = " << nudens_1[0] << " " << nudens_1[1]
                                            << " " << nudens_1[2] << endl;
                        Printer::print_warn(ss.str());
                    }
                }
                assert(isfinite(nudens_0_trap[0]));
                assert(isfinite(nudens_0_trap[1]));
                assert(isfinite(nudens_0_trap[2]));
                assert(isfinite(nudens_1_trap[0]));
                assert(isfinite(nudens_1_trap[1]));
                assert(isfinite(nudens_1_trap[2]));
            }

            // Compute the neutrino black body function assuming fixed temperature and Y_e
            CCTK_REAL nudens_0_thin[3], nudens_1_thin[3];
            ierr = NeutrinoDensity(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &nudens_0_thin[0], &nudens_0_thin[1], &nudens_0_thin[2],
                    &nudens_1_thin[0], &nudens_1_thin[1], &nudens_1_thin[2]);
            assert(!ierr);

            // Correct cross-sections for incoming neutrino energy
            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                // Set the neutrino black body function
                CCTK_REAL nudens_0, nudens_1;
                if (opacity_tau_trap < 0 || tau <= opacity_tau_trap) {
                    nudens_0 = nudens_0_thin[ig];
                    nudens_1 = nudens_1_thin[ig];
                }
                else if (tau > opacity_tau_trap + opacity_tau_delta) {
                    nudens_0 = nudens_0_trap[ig];
                    nudens_1 = nudens_1_trap[ig];
                }
                else {
                    CCTK_REAL const lam = (tau - opacity_tau_trap)/opacity_tau_delta;
                    nudens_0 = lam*nudens_0_trap[ig] + (1 - lam)*nudens_0_thin[ig];
                    nudens_1 = lam*nudens_1_trap[ig] + (1 - lam)*nudens_1_thin[ig];
                }

                // Set the neutrino energies
                nueave[i4D] = nudens_1/nudens_0;

                // Correct absorption opacities for non-LTE effects
                // (kappa ~ E_nu^2)
                CCTK_REAL corr_fac = 1.0;
                corr_fac = (rJ[i4D]/rnnu[i4D])*(nudens_0/nudens_1);
                if (!isfinite(corr_fac)) {
                    corr_fac = 1.0;
                }
                corr_fac *= corr_fac;
                corr_fac = max(1.0/opacity_corr_fac_max, min(corr_fac, opacity_corr_fac_max));

                // Extract scattering opacity
                scat_1[i4D] = corr_fac*(kappa_1_loc[ig] - abs_1_loc[ig]);

                // Enforce Kirchhoff's laws.
                // . For the heavy lepton neutrinos this is implemented by
                //   changing the opacities.
                // . For the electron type neutrinos this is implemented by
                //   changing the emissivities.
                // It would be better to have emissivities and absorptivities
                // that satisfy Kirchhoff's law.
                if (ig == 2) {
                    eta_0[i4D] = corr_fac*eta_0_loc[ig];
                    eta_1[i4D] = corr_fac*eta_1_loc[ig];
                    abs_0[i4D] = (nudens_0 > rad_N_floor ? eta_0[i4D]/nudens_0 : 0);
                    abs_1[i4D] = (nudens_1 > rad_E_floor ? eta_1[i4D]/nudens_1 : 0);
                }
                else {
                    abs_0[i4D] = corr_fac*abs_0_loc[ig];
                    abs_1[i4D] = corr_fac*abs_1_loc[ig];
                    eta_0[i4D] = abs_0[i4D]*nudens_0;
                    eta_1[i4D] = abs_1[i4D]*nudens_1;
                }
            }
        } UTILS_ENDLOOP3(thc_m1_calc_opacity);
    }
    // Done with printing
    thc::Printer::stop();
}
