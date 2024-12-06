#include "stdio.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
/*
 * Set `none` boundary conditions on IGM evolved variables, as these are set internally.
 *  
*/

void select_BCs_IllinoisGRMHD(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

	// Add calls to Boundary_SelectVarForBC to match what's done in Baikal, THC?
	int ierr;

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::rho_star", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::rho_star!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::tau", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::tau!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::mhd_st_x", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::mhd_st_x!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::mhd_st_y", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::mhd_st_y!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::mhd_st_z", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::mhd_st_z!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::Ye_star", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::Ye_star!");	

  ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, cctk_nghostzones[0], -1, "IllinoisGRMHD::S_star", "none");
  if (ierr < 0) CCTK_ERROR("Failed to register BC for IllinoisGRMHD::S_star!");	
}
