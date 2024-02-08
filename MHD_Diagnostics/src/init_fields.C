#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include <iostream>

extern "C"
void init_fields(CCTK_ARGUMENTS) 
{
   DECLARE_CCTK_ARGUMENTS;  // Declare all grid functions from interface.ccl
   DECLARE_CCTK_PARAMETERS; // Declare all parameters from param.ccl

   CCTK_INFO ("Setting  diagnostic fields to zero");
#pragma omp parallel for
   for (int k = 0; k < cctk_lsh[2]; k++) // loop over the z direction
   {
     for (int j = 0; j < cctk_lsh[1]; j++) // loop over the y direction
     {
       for (int i = 0; i < cctk_lsh[0]; i++) // loop over the x direction
       {
         const size_t index  = CCTK_GFINDEX3D(cctkGH, i, j, k);

        smallbt[index] = 0.0;
	smallbx[index] = 0.0;
	smallby[index] = 0.0;
	smallbz[index] = 0.0;
	smallb2[index] = 0.0;
	Poynx[index]   = 0.0;
	Poyny[index]   = 0.0;
	Poynz[index]   = 0.0;
	Poyn2x[index]   = 0.0;
	Poyn2y[index]   = 0.0;
	Poyn2z[index]   = 0.0;
	sqrt_gamma[index] = 0.0;
	minus_one_minus_u_0[index] = 0.0;

	expansion_scalar[index]               = 0.0;
	expansion_no_time_der[index]               = 0.0;
	shear_spatial_tensor_xx[index]        = 0.0;
	shear_spatial_tensor_xy[index]        = 0.0;
	shear_spatial_tensor_xz[index]        = 0.0;
	shear_spatial_tensor_yy[index]        = 0.0;
	shear_spatial_tensor_yz[index]        = 0.0;
	shear_spatial_tensor_zz[index]        = 0.0;
	kin_vorticity_spatial_xy[index]       = 0.0;
	kin_vorticity_spatial_xz[index]       = 0.0;
	kin_vorticity_spatial_yz[index]       = 0.0;
	kin_acceleration_spatial_x[index]     = 0.0;
	kin_acceleration_spatial_y[index]     = 0.0;
        kin_acceleration_spatial_z[index]     = 0.0;

        sigma4bb[index] = 0.0;
	sigma4Ut[index] = 0.0;
	sigma4Ux[index] = 0.0;
	sigma4Uy[index] = 0.0;
	sigma4Uz[index] = 0.0;
	sigma4Trace[index] = 0.0;
	omega4Ut[index] = 0.0;
	omega4Ux[index] = 0.0;
	omega4Uy[index] = 0.0;
	omega4Uz[index] = 0.0;

	normB[index] = 0.0;
	normcurlB[index] = 0.0;
	a4sq[index] = 0.0;

	/*
	expansion_4[index]               = 0.0;
	shear_tt[index]                  = 0.0;
	shear_tx[index]                  = 0.0;
	shear_ty[index]                  = 0.0;
	shear_tz[index]                  = 0.0;
	shear_xx[index]                  = 0.0;
	shear_xy[index]                  = 0.0;
	shear_xz[index]                  = 0.0;
	shear_yy[index]                  = 0.0;
	shear_yz[index]                  = 0.0;
        shear_zz[index]                  = 0.0;
	omega_tx[index]                  = 0.0;
	omega_ty[index]                  = 0.0;
	omega_tz[index]                  = 0.0;
	omega_xy[index]                  = 0.0;
	omega_xz[index]                  = 0.0;
        omega_yz[index]                  = 0.0;
	acc_t[index]                     = 0.0;
	acc_x[index]                     = 0.0;
	acc_y[index]                     = 0.0;
	acc_z[index]                     = 0.0;
*/

       }
     }
   }

   CCTK_INFO ("Setting diagnostic fields to zero was successful");
}
