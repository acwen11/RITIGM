#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_Table.h"
#include "util_String.h"

void Interpolate_density_many_pts(cGH *cctkGH,int interp_num_points,double *particle_x_temp,double *particle_y_temp,double *particle_z_temp, double *particle_density_temp) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_STRING coord_system = "cart3d";

  // Set up handles
  const int coord_system_handle = CCTK_CoordSystemHandle(coord_system);
  if(coord_system_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "can't get coordinate system handle for coordinate system \"%s\"!", coord_system);

  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if(operator_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "couldn't find interpolator \"%s\"!", interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if(param_table_handle < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "bad interpolator parameter(s) \"%s\"!", interpolator_pars);

  CCTK_INT operand_indices[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operand_indices, "operand_indices");

  CCTK_INT operation_codes[1] = {0};
  Util_TableSetIntArray(param_table_handle, 1, operation_codes, "operation_codes");

  const void* interp_coords[3]
    = { (const void *) particle_x_temp,
        (const void *) particle_y_temp,
        (const void *) particle_z_temp };

  // 3d input arrays
  CCTK_STRING input_array_names[1] = { "HydroBase::rho" };

  CCTK_INT input_array_indices[1] = {CCTK_VarIndex(input_array_names[0])};
  if(input_array_indices[0] < 0)
    CCTK_VWARN(CCTK_WARN_ABORT, "COULD NOT FIND VARIABLE '%s'.", input_array_names[0]);

  CCTK_INT output_array_types[1] = {CCTK_VARIABLE_REAL};

  void * output_arrays[1] = { (void *) particle_density_temp };

  // actual interpolation call
  int ierr = CCTK_InterpGridArrays(cctkGH,
                                   3, // number of dimensions
                                   operator_handle,
                                   param_table_handle,
                                   coord_system_handle,
                                   interp_num_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                   1, // Number of input arrays
                                   input_array_indices,
                                   1, // Number of output arrays
                                   output_array_types,
                                   output_arrays);
  if(ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT,"Could not destroy table");
}
