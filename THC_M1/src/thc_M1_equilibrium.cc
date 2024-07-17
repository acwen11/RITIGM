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


#include <cstring>
#include <sstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"
#include "thc_printer.hh"
#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"

using namespace std;
using namespace thc;
using namespace thc::m1;
using namespace utils;

#define CGS_GCC (1.619100425158886e-18)

extern "C" void THC_M1_SetToEquilibrium(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_SetToEquilibrium");
    }

    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, psi_bssn);
    tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz, fidu_w_lorentz,
            fidu_velx, fidu_vely, fidu_velz);

    // Setup Printer
    thc::Printer::start(
            "[INFO|THC|THC_M1_SetToEquilibrium]: ",
            "[WARN|THC|THC_M1_SetToEquilibrium]: ",
            "[ERR|THC|THC_M1_SetToEquilibrium]: ",
            m1_max_num_msg, m1_max_num_msg);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_equilibrium,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Skip this point if it is at low density
            if (rho[ijk] < equilibrium_rho_min*CGS_GCC) {
                continue;
            }

            tensor::metric<4> g_dd;
            tensor::inv_metric<4> g_uu;
            tensor::generic<CCTK_REAL, 4, 1> n_u;
            tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
            geom.get_metric(ijk, &g_dd);
            geom.get_inv_metric(ijk, &g_uu);
            geom.get_normal(ijk, &n_u);
            geom.get_space_proj(ijk, &gamma_ud);

            tensor::generic<CCTK_REAL, 4, 1> u_u;
            tensor::generic<CCTK_REAL, 4, 1> u_d;
            fidu.get(ijk, &u_u);
            tensor::contract(g_dd, u_u, &u_d);

            tensor::generic<CCTK_REAL, 4, 1> F_d;
            tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
            tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;

            assert(nspecies == 3);
            assert(ngroups == 1);

            //
            // Compute the optically thick weak equilibrium
            CCTK_REAL nudens_0[3], nudens_1[3];
            int ierr = NeutrinoDensity(
                    rho[ijk], temperature[ijk], Y_e[ijk],
                    &nudens_0[0], &nudens_0[1], &nudens_0[2],
                    &nudens_1[0], &nudens_1[1], &nudens_1[2]);
            assert(!ierr);

						//
						//Get det(g)
						double volform = std::pow(psi_bssn[ijk], 6);

            for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                //
                // Set to equilibrium in the fluid frame
                rHt[i4D] = 0.0;
                rHx[i4D] = 0.0;
                rHy[i4D] = 0.0;
                rHz[i4D] = 0.0;
                rJ[i4D]  = nudens_1[ig]*volform;
                chi[i4D] = 1.0/3.0;

                //
                // Compute quantities in the lab frame
                for (int a = 0; a < 4; ++a)
                for (int b = 0; b < 4; ++b) {
                    T_dd(a,b) = (4./3.) * rJ[i4D] * u_d(a) * u_d(b) +
                                (1./3.) * rJ[i4D] * g_dd(a,b);
                }
                rE[i4D] = calc_J_from_rT(T_dd, n_u);
                calc_H_from_rT(T_dd, n_u, gamma_ud, &F_d);
                apply_floor(g_uu, &rE[i4D], &F_d);
                unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
                calc_K_from_rT(T_dd, gamma_ud, &P_dd);
                unpack_P_dd(P_dd, &rPxx[i4D], &rPxy[i4D], &rPxz[i4D],
                        &rPyy[i4D], &rPyz[i4D], &rPzz[i4D]);

                //
                // Now compute neutrino number density
                rnnu[i4D] = nudens_0[ig]*volform;
                rN[i4D]   = max(fidu_w_lorentz[ijk]*rnnu[i4D], rad_N_floor);
            }
        } UTILS_ENDLOOP3(thc_m1_equilibrium);
    }
    // Done with printing
    thc::Printer::stop();
}
