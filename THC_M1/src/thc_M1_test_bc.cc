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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#define SQ(X) ((X)*(X))

using namespace utils;
using namespace std;

extern "C" void THC_M1_KerrBCs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_KerrBCs");
    }

    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);

    if (cctk_bbox[0]) {
        for (int k = 0; k < cctk_lsh[2]; ++k)
        for (int j = 0; j < cctk_lsh[1]; ++j)
        for (int i = 0; i < cctk_nghostzones[0]; ++i) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            tensor::metric<4> g_dd;
            tensor::generic<CCTK_REAL, 4, 1> beta_u;
            tensor::generic<CCTK_REAL, 4, 1> beta_d;
            geom.get_metric(ijk, &g_dd);
            geom.get_shift_vec(ijk, &beta_u);
            tensor::contract(g_dd, beta_u, &beta_d);

            // Here we want to satisfy two conditions
            //   F_i F^i = (1 - eps) E^2
            //   alp F^i = beta^i E + a \delta^i_x E
            // The coefficient a can be found as solution of
            //   g_xx a^2 + 2 a beta_x + [beta^2 - alp^2 (1 - eps)] = 0
            CCTK_REAL const eps = 0.01;
            CCTK_REAL const g_xx = g_dd(1,1);
            CCTK_REAL const beta_x = beta_d(1);
            CCTK_REAL const beta2 = tensor::dot(beta_u, beta_d);
            CCTK_REAL const a = (-beta_x + sqrt(SQ(beta_x) - beta2 +
                        SQ(alp[ijk])*(1 - eps)))/g_xx;

            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                if (abs(y[ijk]) <= kerr_beam_width &&
                        z[ijk] >= kerr_beam_position &&
                        z[ijk] <= kerr_beam_position + kerr_beam_width) {
                    CCTK_REAL const E = volform[ijk];

                    tensor::generic<CCTK_REAL, 4, 1> F_u;
                    F_u(0) = 0;
                    F_u(1) = a*E/alp[ijk] + beta_u(1)*E/alp[ijk];
                    F_u(2) = beta_u(2)*E/alp[ijk];
                    F_u(3) = beta_u(3)*E/alp[ijk];

                    tensor::generic<CCTK_REAL, 4, 1> F_d;
                    tensor::contract(g_dd, F_u, &F_d);

                    rE[i4D] = E;
                    rFx[i4D] = F_d(1);
                    rFy[i4D] = F_d(2);
                    rFz[i4D] = F_d(3);
                    rN[i4D] = 1.0;
                }
                else {
                    rE[i4D] = 0.0;
                    rFx[i4D] = 0.0;
                    rFy[i4D] = 0.0;
                    rFz[i4D] = 0.0;
                    rN[i4D] = 0.0;
                }
            }
        }
    }
    if (cctk_bbox[1]) {
        for (int k = 0; k < cctk_lsh[2]; ++k)
        for (int j = 0; j < cctk_lsh[1]; ++j)
        for (int i = cctk_lsh[0]-cctk_nghostzones[0]; i < cctk_lsh[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
        }
    }
    if (cctk_bbox[2]) {
        for (int k = 0; k < cctk_lsh[2]; ++k)
        for (int j = 0; j < cctk_nghostzones[1]; ++j)
        for (int i = 0; i < cctk_lsh[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4Db = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            int const i4Di = CCTK_VectGFIndex3D(cctkGH,
                    i, cctk_nghostzones[1], k, ig);
            rE[i4Db] = rE[i4Di];
            rFx[i4Db] = rFx[i4Di];
            rFy[i4Db] = rFy[i4Di];
            rFz[i4Db] = rFz[i4Di];
            rN[i4Db] = rN[i4Di];
        }
    }
    if (cctk_bbox[3]) {
        for (int k = 0; k < cctk_lsh[2]; ++k)
        for (int j = cctk_lsh[1]-cctk_nghostzones[1]; j < cctk_lsh[1]; ++j)
        for (int i = 0; i < cctk_lsh[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4Db = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            int const i4Di = CCTK_VectGFIndex3D(cctkGH,
                    i, cctk_lsh[1] - cctk_nghostzones[1] - 1, k, ig);
            rE[i4Db] = rE[i4Di];
            rFx[i4Db] = rFx[i4Di];
            rFy[i4Db] = rFy[i4Di];
            rFz[i4Db] = rFz[i4Di];
            rN[i4Db] = rN[i4Di];
        }
    }
    if (cctk_bbox[4]) {
        for (int k = 0; k < cctk_nghostzones[2]; ++k)
        for (int j = 0; j < cctk_lsh[1]; ++j)
        for (int i = 0; i < cctk_lsh[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
        }
    }
    if (cctk_bbox[5]) {
        for (int k = cctk_lsh[2]-cctk_nghostzones[2]; k < cctk_lsh[2]; ++k)
        for (int j = 0; j < cctk_lsh[1]; ++j)
        for (int i = 0; i < cctk_lsh[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
        }
    }
}

extern "C" void THC_M1_ShadowBCs(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_ShadowBCs");
    }

    if (cctk_bbox[0]) {
        for (int k = 0; k < cctk_lsh[2]; ++k)
        for (int j = 0; j < cctk_lsh[1]; ++j)
        for (int i = 0; i < cctk_nghostzones[0]; ++i)
        for (int ig = 0; ig < ngroups*nspecies; ++ig) {
            int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
            if (abs(y[CCTK_GFINDEX3D(cctkGH, i, j, k)]) < 1.5) {
                rE[i4D] = 1.0;
                rFx[i4D] = 1.0;
                rFy[i4D] = 0.0;
                rFz[i4D] = 0.0;
                rN[i4D] = 1.0;
            }
        }
    }
}
