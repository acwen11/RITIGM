//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2021, David Radice <david.radice@psu.edu>
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

#include "thc_M1_macro.h"

#include "utils.hh"

using namespace std;

extern "C" void THC_M1_Analysis(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_Analysis");
    }

    CCTK_REAL const mb = AverageBaryonMass();

    assert(nspecies == 3);
    assert(ngroups == 1);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_m1_analysis,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk  = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int const ijke = CCTK_VectGFIndex3D(cctkGH, i, j, k, 0);
            int const ijka = CCTK_VectGFIndex3D(cctkGH, i, j, k, 1);
            int const ijkx = CCTK_VectGFIndex3D(cctkGH, i, j, k, 2);

            CCTK_REAL const nb = rho[ijk]/mb;
            ynue[ijk] = rnnu[ijke]/volform[ijk]/nb;
            ynua[ijk] = rnnu[ijka]/volform[ijk]/nb;
            ynux[ijk] = rnnu[ijkx]/volform[ijk]/nb;

            CCTK_REAL const egas = rho[ijk]*(1 + eps[ijk]);
            CCTK_REAL const enue = rJ[ijke]/volform[ijk];
            CCTK_REAL const enua = rJ[ijka]/volform[ijk];
            CCTK_REAL const enux = rJ[ijkx]/volform[ijk];
            CCTK_REAL const etot = egas + enue + enua + enux;
            znue[ijk] = enue/etot;
            znua[ijk] = enua/etot;
            znux[ijk] = enux/etot;
        } UTILS_ENDLOOP3(thc_m1_analysis);
    }
}
