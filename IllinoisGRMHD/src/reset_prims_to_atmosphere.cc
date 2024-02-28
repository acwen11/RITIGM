#include "cctk.h"

#include "IllinoisGRMHD_headers.h"

void reset_prims_to_atmosphere( const igm_eos_parameters eos,
																const CCTK_REAL rho_atm,
																const CCTK_REAL T_atm,
																const CCTK_REAL P_atm,
																const CCTK_REAL eps_atm,
																const CCTK_REAL S_atm,
                                CCTK_REAL *restrict PRIMS ) {

  // Just a simple reset to atmospheric values.
  // Velocities are set to zero. 
  // For a radial atmopheric falloff, these values are now position
  // dependent.

  PRIMS[RHOB         ] = rho_atm;
  PRIMS[PRESSURE     ] = P_atm;
  PRIMS[EPSILON      ] = eps_atm;
  PRIMS[ENTROPY      ] = S_atm;

	// Y_e is still assumed constant.
  if( eos.is_Tabulated ) {
    PRIMS[YEPRIM     ] = eos.Ye_atm;
    PRIMS[TEMPERATURE] = T_atm;
  }
  PRIMS[VX           ] = 0;
  PRIMS[VY           ] = 0;
  PRIMS[VZ           ] = 0;  
  
}
