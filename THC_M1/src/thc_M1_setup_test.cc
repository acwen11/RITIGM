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


#include <cmath>
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"
#include "utils.hh"

using namespace std;
using namespace thc::m1;
using namespace utils;

namespace {

// Estimate the fractional volume of the intersection between a cell and a
// sphere of radius R with center in the origin
CCTK_REAL volume(CCTK_REAL R, CCTK_REAL xp, CCTK_REAL yp, CCTK_REAL zp,
        CCTK_REAL dx, CCTK_REAL dy, CCTK_REAL dz) {
    int const NPOINTS = 10;

    int inside = 0;
    int count = 0;
    for (int i = 0; i < NPOINTS; ++i) {
        CCTK_REAL const myx = (xp - dx/2.) + (i + 0.5)*(dx/NPOINTS);
        for (int j = 0; j < NPOINTS; ++j) {
            CCTK_REAL const myy = (yp - dy/2.) + (j + 0.5)*(dy/NPOINTS);
            for (int k = 0; k < NPOINTS; ++k) {
                CCTK_REAL const myz = (zp - dz/2.) + (k + 0.5)*(dz/NPOINTS);
                count++;
                if (SQ(myx) + SQ(myy) + SQ(myz) < SQ(R)) {
                    inside++;
                }
            }
        }
    }
    return static_cast<CCTK_REAL>(inside)/static_cast<CCTK_REAL>(count);
}

} // namespace

extern "C" void THC_M1_SetupTest(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_SetupTest");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_setup_test,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if (CCTK_Equals(thc_m1_test, "beam")) {
                CCTK_REAL nx = beam_test_dir[0];
                CCTK_REAL ny = beam_test_dir[1];
                CCTK_REAL nz = beam_test_dir[2];
                CCTK_REAL n2 = SQ(nx) + SQ(ny) + SQ(nz);
                if (n2 > 0) {
                    CCTK_REAL nn = sqrt(n2);
                    nx /= nn;
                    ny /= nn;
                    nz /= nn;
                }
                else {
                    nx = 1.0;
                    nz = ny = 0.0;
                }

                CCTK_REAL proj = nx*x[ijk] + ny*y[ijk] + nz*z[ijk];
                CCTK_REAL offset2 = SQ(x[ijk] - nx*x[ijk]) +
                    SQ(y[ijk] - ny*y[ijk]) + SQ(z[ijk] - nz*z[ijk]);
                if (proj < beam_position && offset2 < SQ(beam_width)) {
                    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                        int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                        rE[i4D]   = 1.0;
                        rN[i4D]   = 1.0;
                        rFx[i4D]  = nx;
                        rFy[i4D]  = ny;
                        rFz[i4D]  = nz;
                    }
                }
                else {
                    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                        int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                        rE[i4D] = 0.0;
                        rN[i4D] = 0.0;
                        rFx[i4D] = 0.0;
                        rFy[i4D] = 0.0;
                        rFz[i4D] = 0.0;
                    }
                }
            }
            else if (CCTK_Equals(thc_m1_test, "diff")) {
                for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                    if (CCTK_Equals(diff_profile, "step")) {
                        if (x[ijk] > -0.5 && x[ijk] < 0.5) {
                            rE[i4D] = 1.0;
                        }
                        else {
                            rE[i4D] = 0.0;
                        }
                    }
                    else if (CCTK_Equals(diff_profile, "gaussian")) {
                        rE[i4D] = exp(-SQ(3*x[ijk]));
                    }
                    rN[i4D] = rE[i4D];
                    // Use thick closure to compute the fluxes
                    CCTK_REAL const W = fidu_w_lorentz[ijk];
                    CCTK_REAL const Jo3 = rE[i4D]/(4*SQ(W) - 1);
                    rFx[i4D] = 4*SQ(W)*fidu_velx[ijk]*Jo3;
                    rFy[i4D] = 4*SQ(W)*fidu_vely[ijk]*Jo3;
                    rFz[i4D] = 4*SQ(W)*fidu_velz[ijk]*Jo3;
                }
            }
            else if (CCTK_Equals(thc_m1_test, "equil")) {
                assert(ngroups == 1);
                assert(nspecies == 3);
                CCTK_REAL const W = fidu_w_lorentz[ijk];
                for (int is = 0; is < nspecies; ++is) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, is);
                    CCTK_REAL const Jnu = equil_nudens_1[is];
                    rE[i4D]  = (4.*W*W - 1.)/3.*Jnu;
                    rFx[i4D] = 4./3.*SQ(W)*fidu_velx[ijk]*Jnu;
                    rFy[i4D] = 4./3.*SQ(W)*fidu_vely[ijk]*Jnu;
                    rFz[i4D] = 4./3.*SQ(W)*fidu_velz[ijk]*Jnu;
                    rN[i4D]  = equil_nudens_0[is]*W;
                }
            }
            else if (CCTK_Equals(thc_m1_test, "kerrschild") ||
                     CCTK_Equals(thc_m1_test, "shadow") ||
                     CCTK_Equals(thc_m1_test, "sphere")) {
                for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                    rE[i4D] = 0.0;
                    rN[i4D] = 0.0;
                    rFx[i4D] = 0.0;
                    rFy[i4D] = 0.0;
                    rFz[i4D] = 0.0;
                }
            }
            else {
                char msg[BUFSIZ];
                snprintf(msg, BUFSIZ, "Unknown test problem \"%s\"",
                        thc_m1_test);
                CCTK_ERROR(msg);
            }
        } UTILS_ENDLOOP3(thc_m1_setup_test);
    }

}

extern "C" void THC_M1_KerrSchild_Mask(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_KerrSchild_Mask");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_kerrschild_mask,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            if (SQ(x[ijk]) + SQ(y[ijk]) + SQ(z[ijk]) < SQ(kerr_mask_radius)) {
                thc_m1_mask[ijk] = 1;
                for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                    rN[i4D] = 0;
                    rE[i4D] = 0;
                    rFx[i4D] = 0;
                    rFy[i4D] = 0;
                    rFz[i4D] = 0;
                }
            }
            else {
              thc_m1_mask[ijk] = 0;
            }
        } UTILS_ENDLOOP3(thc_m1_kerrschild_mask);
    }
}

extern "C" void THC_M1_SetupTest_Hydro(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_SetupTest_Hydro");
    }

    CCTK_REAL dx = CCTK_DELTA_SPACE(0);
    CCTK_REAL dy = CCTK_DELTA_SPACE(1);
    CCTK_REAL dz = CCTK_DELTA_SPACE(2);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_setup_test_hydro,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            if (CCTK_Equals(thc_m1_test, "shadow") || CCTK_Equals(thc_m1_test, "sphere")) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                rho[ijk] = volume(1.0, x[ijk], y[ijk], z[ijk], dx, dy, dz);
                vel[CCTK_VectGFIndex3D(cctkGH, i, j, k, 0)] = 0.0;
                vel[CCTK_VectGFIndex3D(cctkGH, i, j, k, 1)] = 0.0;
                vel[CCTK_VectGFIndex3D(cctkGH, i, j, k, 2)] = 0.0;
            }
            else {
                char msg[BUFSIZ];
                snprintf(msg, BUFSIZ, "Unknown test problem \"%s\"",
                        thc_m1_test);
                CCTK_ERROR(msg);
            }
        } UTILS_ENDLOOP3(thc_m1_setup_test_hydro);
    }
}
