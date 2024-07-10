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
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#define CGS_GCC (1.619100425158886e-18)

extern "C" void THC_M1_FiducialVelocity(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_FiducialVelocity");
    }

    size_t siz = UTILS_GFSIZE(cctkGH);

    CCTK_REAL const * velx = &vel[0*siz];
    CCTK_REAL const * vely = &vel[1*siz];
    CCTK_REAL const * velz = &vel[2*siz];

    if (CCTK_Equals(fiducial_velocity, "fluid")) {
        std::memcpy(fidu_velx, velx, siz*sizeof(CCTK_REAL));
        std::memcpy(fidu_vely, vely, siz*sizeof(CCTK_REAL));
        std::memcpy(fidu_velz, velz, siz*sizeof(CCTK_REAL));
        std::memcpy(fidu_w_lorentz, w_lorentz, siz*sizeof(CCTK_REAL));
    }
    else if (CCTK_Equals(fiducial_velocity, "mixed")) {
        UTILS_LOOP3(fidu_vel_mixed,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL fac = 1.0/std::max(dens[ijk],
                    fiducial_velocity_rho_fluid*CGS_GCC);

            fidu_velx[ijk] = velx[ijk]*dens[ijk]*fac;
            fidu_vely[ijk] = vely[ijk]*dens[ijk]*fac;
            fidu_velz[ijk] = velz[ijk]*dens[ijk]*fac;

            CCTK_REAL const fidu_vel_x =
                gxx[ijk]*fidu_velx[ijk] +
                gxy[ijk]*fidu_vely[ijk] +
                gxz[ijk]*fidu_velz[ijk];
            CCTK_REAL const fidu_vel_y =
                gxy[ijk]*fidu_velx[ijk] +
                gyy[ijk]*fidu_vely[ijk] +
                gyz[ijk]*fidu_velz[ijk];
            CCTK_REAL const fidu_vel_z =
                gxz[ijk]*fidu_velx[ijk] +
                gyz[ijk]*fidu_vely[ijk] +
                gzz[ijk]*fidu_velz[ijk];

            CCTK_REAL const v2 =
                fidu_vel_x*fidu_velx[ijk] +
                fidu_vel_y*fidu_vely[ijk] +
                fidu_vel_z*fidu_velz[ijk];
            fidu_w_lorentz[ijk] = 1.0/sqrt(1.0 - v2);
        } UTILS_ENDLOOP3(fidu_vel_mixed);
    }
    else {
        std::fill(fidu_velx, fidu_velx + siz, 0.0);
        std::fill(fidu_vely, fidu_vely + siz, 0.0);
        std::fill(fidu_velz, fidu_velz + siz, 0.0);
        std::fill(fidu_w_lorentz, fidu_w_lorentz + siz, 1.0);
    }

}

