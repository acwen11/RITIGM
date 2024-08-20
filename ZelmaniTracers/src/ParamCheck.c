#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void ZelmaniTracers_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(seed_method, "on_grid") || manual_pos_tracers) {
        if(CCTK_Equals(TracerProblem, "BinaryNS") &&
                !CCTK_Equals(symm, "quadrant")) {
            if(ntracers != 2*nradialshells*ntheta*nphi) {
                CCTK_PARAMWARN("ntracers should be equal to 2*nradialshells*"
                        "ntheta*nphi");
            }
        }
        else if(ntracers != nradialshells*ntheta*nphi) {
            CCTK_PARAMWARN("ntracers should be equal to nradialshells*"
                    "ntheta*nphi!");
        }
    }
}
