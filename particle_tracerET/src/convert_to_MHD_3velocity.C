#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/* The 3-velocity as defined in IllinoisGRMHD: v^i = u^i/u^0 is the definition of velocity consistent with that used in the B-field induction equation.
   We adopt this definition here, since using other definitions will result in the particles not following the frozen-in MHD condition.
   This function converts the HydroBase/Valencia definition of 3-velocity to the IllinoisGRMHD one.
   It needs HydroBase vel[] gridfunctions and ADMBase's alp & betai gridfunctions. */
void convert_to_MHD_3velocity(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

	double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

	// IllinoisGRMHD defines v^i = u^i/u^0, consistent with the
	// magnetic induction equation 3-velocity.

	// Meanwhile, the ET/HydroBase formalism, called the Valencia
	// formalism, splits the 4 velocity into a purely spatial part
	// and a part that is normal to the spatial hypersurface:
	// u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
	// where n^a is the unit normal vector to the spatial hypersurface,
	// n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
	// is defined in HydroBase as the vel[] vector gridfunction.
	// Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
	// of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
	// the standard Lorentz factor.

	// Note that n^i = - \beta^i / \alpha, so
	// u^a = \Gamma (n^a + U^a)
	// -> u^i = \Gamma ( U^i - \beta^i / \alpha )
	// which implies
	// v^i = u^i/u^0
	//     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
	//     = \alpha ( U^i - \beta^i / \alpha )
	//     = \alpha U^i - \beta^i
	MHDvx[index] = alp[index]*ETvx - betax[index];
	MHDvy[index] = alp[index]*ETvy - betay[index];
	MHDvz[index] = alp[index]*ETvz - betaz[index];

      }
}
