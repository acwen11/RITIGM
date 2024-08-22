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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void THC_M1_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (CCTK_Equals(thc_m1_test, "const")) {
        if (ngroups != 1 || nspecies != 3) {
            CCTK_PARAMWARN("If \"THC_M1::thc_m1_test\" is set to \"const\", then"
                    "\"THC_M1::ngroups\" must be 1 and \"THC_M1::nspecies\" must be 3");
        }
    }

    if (CCTK_Equals(thc_m1_test, "shadow") ||
        CCTK_Equals(thc_m1_test, "sphere")) {
        if (!CCTK_Equals(initial_hydro, "THC_M1")) {
            CCTK_PARAMWARN("This test requires HydroBase::initial_data "
                    "to be \"THC_M1\"");
        }
    }

    if (optimize_prolongation) {
        if (cctk_nghostzones[0] < 4 ||
            cctk_nghostzones[1] < 4 ||
            cctk_nghostzones[2] < 4) {
            CCTK_PARAMWARN("THC_M1::optimize_prolongation requires at least "
                    "four ghost points");
        }
    }
}
