#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_Table.h"
#include "util_String.h"

//extern void Interpolate_density_many_pts(cGH *cctkGH,int interp_num_points,double *particle_x_temp,double *particle_y_temp,double *particle_z_temp, double *particle_density_temp);
//
void Interpolate_density_many_pts(cGH *cctkGH,int interp_num_points,double *particle_x_temp,double *particle_y_temp,double *particle_z_temp, double *particle_density_temp, double* particle_volform_temp) {
  DECLARE_CCTK_PARAMETERS;

  double *wl_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gxx_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gxy_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gxz_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gyy_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gyz_temp  = (double *)malloc(sizeof(double)*interp_num_points);
  double *gzz_temp  = (double *)malloc(sizeof(double)*interp_num_points);

  // Set up handles
  const int coord_system_handle = CCTK_CoordSystemHandle(coordinate_system);
  if(coord_system_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "can't get coordinate system handle for coordinate system \"%s\"!", coordinate_system);

  const int operator_handle = CCTK_InterpHandle(interpolator);
  if(operator_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "couldn't find interpolator \"%s\"!", interpolator);

  int param_table_handle = Util_TableCreateFromString(interpolator_options);
  if(param_table_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "bad interpolator parameter(s) \"%s\"!", interpolator_options);

  CCTK_INT operand_indices[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operand_indices, "operand_indices");

  CCTK_INT operation_codes[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operation_codes, "operation_codes");

  const void* interp_coords[3]
    = { (const void *) particle_x_temp,
        (const void *) particle_y_temp,
        (const void *) particle_z_temp };

  // 3d input arrays
  CCTK_STRING input_array_names[8] = { "HydroBase::rho",
																			"HydroBase::w_lorentz",
																			"ADMBase::gxx",
																			"ADMBase::gxy",
																			"ADMBase::gxz",
																			"ADMBase::gyy",
																			"ADMBase::gyz",
																			"ADMBase::gzz"};

  CCTK_INT input_array_indices[8];
  for(int i=0; i<8; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0)
      CCTK_VWARN(CCTK_WARN_ABORT, "COULD NOT FIND VARIABLE '%s'.", input_array_names[i]);
  }

  CCTK_INT output_array_types[8];
  for(int i=0; i<8; i++)
    output_array_types[i] = CCTK_VARIABLE_REAL;

  void * output_arrays[8] = { (void *) particle_density_temp,
															(void *) wl_temp,
															(void *) gxx_temp,
															(void *) gxy_temp,
															(void *) gxz_temp,
															(void *) gyy_temp,
															(void *) gyz_temp,
															(void *) gzz_temp};

  // actual interpolation call
  int ierr = CCTK_InterpGridArrays(cctkGH,
                                   3, // number of dimensions
                                   operator_handle,
                                   param_table_handle,
                                   coord_system_handle,
                                   interp_num_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                   8, // Number of input arrays
                                   input_array_indices,
                                   8, // Number of output arrays
                                   output_array_types,
                                   output_arrays);

  if(ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT,"Could not destroy table");
	for (int ii=0; ii<interp_num_points; ii++) {
		const CCTK_REAL detg = -gxz_temp[ii]*gxz_temp[ii]*gyy_temp[ii]
													 + 2*gxy_temp[ii]*gxz_temp[ii]*gyz_temp[ii]
													 - gxx_temp[ii]*gyz_temp[ii]*gyz_temp[ii]
													 - gxy_temp[ii]*gxy_temp[ii]*gzz_temp[ii] 
													 + gxx_temp[ii]*gyy_temp[ii]*gzz_temp[ii]; 
		particle_volform_temp[ii] = wl_temp[ii] * sqrt(fabs(detg));
	}
}
