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

extern "C" void THC_M1_InitTimeIntegrator(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_InitTimeIntegrator");
    }

    *TimeIntegratorStage = 2;
}

extern "C" void THC_M1_FinalizeTimeIntegrator(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (!QueryProlongating() || *TimeIntegratorStage != 0) {
        CCTK_VERROR("Unexpected prolongation state %d or time integrator stage %d. Expected 'true' and 0.",
                    (int)QueryProlongating(), (int)*TimeIntegratorStage);
    }

    if(verbose) {
        CCTK_INFO("THC_M1_FinalizeTimeIntegrator");
    }
}
