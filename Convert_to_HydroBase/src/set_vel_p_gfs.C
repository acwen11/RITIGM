#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include <iostream>

extern "C"
void init_vel_p_gfs(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;  // Declare all grid functions from interface.ccl
   DECLARE_CCTK_PARAMETERS; // Declare all parameters from param.ccl

   CCTK_INFO ("Setting vel_p fields to zero");
#pragma omp parallel for
   for (int k = 0; k < cctk_lsh[2]; k++) // loop over the z direction
   {
     for (int j = 0; j < cctk_lsh[1]; j++) // loop over the y direction
     {
       for (int i = 0; i < cctk_lsh[0]; i++) // loop over the x direction
       {
         const size_t index  = CCTK_GFINDEX3D(cctkGH, i, j, k);
					prev_vx[index] = 0.0;
					prev_vy[index] = 0.0;
					prev_vz[index] = 0.0;
				}
			}
		}
}

extern "C"
void set_vel_p_gfs(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;  // Declare all grid functions from interface.ccl
   DECLARE_CCTK_PARAMETERS; // Declare all parameters from param.ccl

	if(Convert_to_HydroBase_every==0 || cctk_iteration%Convert_to_HydroBase_every!=0 || !Need_vel_p) return;

   CCTK_INFO ("Setting vel_p to IGM vel.");
#pragma omp parallel for
   for (int k = 0; k < cctk_lsh[2]; k++) // loop over the z direction
   {
     for (int j = 0; j < cctk_lsh[1]; j++) // loop over the y direction
     {
       for (int i = 0; i < cctk_lsh[0]; i++) // loop over the x direction
       {
         const size_t index  = CCTK_GFINDEX3D(cctkGH, i, j, k);
					prev_vx[index] = vx[index];
					prev_vy[index] = vy[index];
					prev_vz[index] = vz[index];
				}
			}
		}
}
