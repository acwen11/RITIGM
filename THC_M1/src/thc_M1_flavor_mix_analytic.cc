//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2023, David Radice <david.radice@psu.edu>
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
#include "utils_macro.h"
#include "utils_tensor.hh"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"

using namespace utils;
using namespace thc;
using namespace std;
using namespace thc::m1;

extern "C" void THC_M1_FlavorMixAnalytic(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_FlavorMixAnalytic");
    }

    // Setup Printer
    thc::Printer::start(
            "[INFO|THC|THC_M1_MaxFlavorMix]: ",
            "[WARN|THC|THC_M1_MaxFlavorMix]: ",
            "[ERR|THC|THC_M1_MaxFlavorMix]: ",
            m1_max_num_msg, m1_max_num_msg);

    assert(ngroups == 1);
    assert(nspecies == 4);

    int flavor_mix_type = -1;
    if (CCTK_Equals(flavor_mix, "equilibrium")) {
        flavor_mix_type = 0;
    }
    else if (CCTK_Equals(flavor_mix, "maximal")) {
        flavor_mix_type = 1;
    }
    else {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Unknown flavor mix type: \"%s\"", flavor_mix);
        CCTK_ERROR(msg);
    }

    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, psi_bssn);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_flavor_mix,
                k, THC_M1_NGHOST, cctk_lsh[2]-THC_M1_NGHOST,
                j, THC_M1_NGHOST, cctk_lsh[1]-THC_M1_NGHOST,
                i, THC_M1_NGHOST, cctk_lsh[0]-THC_M1_NGHOST) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            if (flavor_mix_rho >= 0.0 && rho[ijk] > CGS_GCC*flavor_mix_rho) {
                continue;
            }

            int const ijke = CCTK_VectGFIndex3D(cctkGH, i, j, k, 0);
            int const ijka = CCTK_VectGFIndex3D(cctkGH, i, j, k, 1);
            int const ijkx = CCTK_VectGFIndex3D(cctkGH, i, j, k, 2);
            int const ijky = CCTK_VectGFIndex3D(cctkGH, i, j, k, 3);

            CCTK_REAL const ne = max(rad_N_floor, rN[ijke]);
            CCTK_REAL const na = max(rad_N_floor, rN[ijka]);
            CCTK_REAL const nx = max(rad_N_floor, rN[ijkx]);
            CCTK_REAL const ny = max(rad_N_floor, rN[ijky]);

            CCTK_REAL const N  = ne + na + nx + ny;
            CCTK_REAL const Ne = ne - na;
            CCTK_REAL const Nx = nx - ny;

            // Maximally mixed states (possibly with negative occupation number)
            CCTK_REAL ne_mm, na_mm, nx_mm, ny_mm;
            if (0 == flavor_mix_type) {
                ne_mm = -N/6 + Ne/2 + sqrt(4*N*N + 12*Ne*Ne - 3*Nx*Nx)/6;
                na_mm = ne_mm - Ne;
                nx_mm = 0.5*(N + Ne + Nx) - ne_mm;
                ny_mm = nx_mm - Nx;
            }
            else {
                ne_mm = N/6 + Ne/2;
                na_mm = N/6 - Ne/2;
                nx_mm = N/3 + Nx/2;
                ny_mm = N/3 - Nx/2;
            }

            // Now enforce positivity
            // Basically we are reducing the number of neutrinos that can
            // participate in the flavor mixing
            CCTK_REAL alpha = 0;
            if (ne_mm < 0) {
                alpha = max(alpha, -ne_mm/(ne - ne_mm));
            }
            if (na_mm < 0) {
                alpha = max(alpha, -na_mm/(na - na_mm));
            }
            if (nx_mm < 0) {
                alpha = max(alpha, -nx_mm/(nx - nx_mm));
            }
            if (ny_mm < 0) {
                alpha = max(alpha, -ny_mm/(ny - ny_mm));
            }
            alpha = min(1.0, alpha);

            // Flavor mixed number densities
            CCTK_REAL const ne_mix = max(rad_N_floor, alpha*ne + (1 - alpha)*ne_mm);
            CCTK_REAL const na_mix = max(rad_N_floor, alpha*na + (1 - alpha)*na_mm);
            CCTK_REAL const nx_mix = max(rad_N_floor, alpha*nx + (1 - alpha)*nx_mm);
            CCTK_REAL const ny_mix = max(rad_N_floor, alpha*ny + (1 - alpha)*ny_mm);

            // Save the results back
            rN[ijke] = ne_mix;
            rN[ijka] = na_mix;
            rN[ijkx] = nx_mix;
            rN[ijky] = ny_mix;

            // Compute transition probabilities for neutrinos (needed to update rE and rF)
            CCTK_REAL const P_e_to_e = min(1., ne_mix/ne);
            CCTK_REAL const P_e_to_x = 1. - P_e_to_e;
            CCTK_REAL const P_x_to_x = min(1., nx_mix/nx);
            CCTK_REAL const P_x_to_e = 1. - P_x_to_x;

            // Transition probabilities for anti-neutrinos
            CCTK_REAL const P_a_to_a = min(1., na_mix/na);
            CCTK_REAL const P_a_to_y = 1. - P_a_to_a;
            CCTK_REAL const P_y_to_y = min(1., ny_mix/ny);
            CCTK_REAL const P_y_to_a = 1. - P_y_to_y;

            // Update energy densities
            rE[ijke] = rE[ijke]*(P_e_to_e - P_e_to_x) + rE[ijkx]*P_x_to_e;
            rE[ijkx] = rE[ijke]*P_e_to_x + rE[ijkx]*(P_x_to_x - P_x_to_e);
            rE[ijka] = rE[ijka]*(P_a_to_a - P_a_to_y) + rE[ijky]*P_y_to_a;
            rE[ijky] = rE[ijka]*P_a_to_y + rE[ijky]*(P_y_to_y - P_y_to_a);

            // Same for the fluxes
            rFx[ijke] = rFx[ijke]*(P_e_to_e - P_e_to_x) + rFx[ijkx]*P_x_to_e;
            rFx[ijkx] = rFx[ijke]*P_e_to_x + rFx[ijkx]*(P_x_to_x - P_x_to_e);
            rFx[ijka] = rFx[ijka]*(P_a_to_a - P_a_to_y) + rFx[ijky]*P_y_to_a;
            rFx[ijky] = rFx[ijka]*P_a_to_y + rFx[ijky]*(P_y_to_y - P_y_to_a);

            rFy[ijke] = rFy[ijke]*(P_e_to_e - P_e_to_x) + rFy[ijkx]*P_x_to_e;
            rFy[ijkx] = rFy[ijke]*P_e_to_x + rFy[ijkx]*(P_x_to_x - P_x_to_e);
            rFy[ijka] = rFy[ijka]*(P_a_to_a - P_a_to_y) + rFy[ijky]*P_y_to_a;
            rFy[ijky] = rFy[ijka]*P_a_to_y + rFy[ijky]*(P_y_to_y - P_y_to_a);

            rFz[ijke] = rFz[ijke]*(P_e_to_e - P_e_to_x) + rFz[ijkx]*P_x_to_e;
            rFz[ijkx] = rFz[ijke]*P_e_to_x + rFz[ijkx]*(P_x_to_x - P_x_to_e);
            rFz[ijka] = rFz[ijka]*(P_a_to_a - P_a_to_y) + rFz[ijky]*P_y_to_a;
            rFz[ijky] = rFz[ijka]*P_a_to_y + rFz[ijky]*(P_y_to_y - P_y_to_a);

            // Now we apply the floors
            tensor::metric<4> g_dd;
            tensor::inv_metric<4> g_uu;
            geom.get_metric(ijk, &g_dd);
            geom.get_inv_metric(ijk, &g_uu);

            tensor::generic<CCTK_REAL, 4, 1> F_d;

            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                         rFx[i4D], rFy[i4D], rFz[i4D], &F_d);
                apply_floor(g_uu, &rE[i4D], &F_d);
                unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
            }
        } UTILS_ENDLOOP3(thc_m1_flavor_mix);
    }
}
