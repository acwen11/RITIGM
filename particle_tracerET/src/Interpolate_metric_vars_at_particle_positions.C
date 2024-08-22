#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_Table.h"
#include "util_String.h"

// Interpolate g_ij and lapse in the same manner as this thorn
// handles outputting 4-velocities.
//
// We will interpolate a total of 7 quantities:
#define NUM_INPUT_ARRAYS  2
#define NUM_OUTPUT_ARRAYS 2

void Interpolate_metric_vars_at_particle_positions(
      cGH *cctkGH,
      int interp_num_points,
      double *particle_x,
      double *particle_y,
      double *particle_z,
      double *particle_psi_bssn,
      double *particle_alp) {

  DECLARE_CCTK_PARAMETERS;

  const int coord_system_handle = CCTK_CoordSystemHandle("cart3d");
  if(coord_system_handle < 0) CCTK_WARN(CCTK_WARN_ABORT, "can't get coordinate system handle for coordinate system \"cart3d\"!");

  const int operator_handle = CCTK_InterpHandle(interpolator_name);
  if(operator_handle < 0) CCTK_VWARN(CCTK_WARN_ABORT, "couldn't find interpolator \"%s\"!", interpolator_name);

  int param_table_handle = Util_TableCreateFromString(interpolator_pars);
  if(param_table_handle < 0) CCTK_VWARN(CCTK_WARN_ABORT, "bad interpolator parameter(s) \"%s\"!", interpolator_pars);

  const void* interp_coords[3]
    = { (const void *) particle_x,
        (const void *) particle_y,
        (const void *) particle_z };

  CCTK_STRING input_array_names[NUM_INPUT_ARRAYS]
    = { "IllinoisGRMHD::psi_bssn",
        "ADMBase::alp" };

  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];
  for(int i=0; i<NUM_INPUT_ARRAYS; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0)
      CCTK_VWARN(CCTK_WARN_ABORT, "COULD NOT FIND VARIABLE '%s'.", input_array_names[i]);
  }

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];
  for(int i=0; i<NUM_OUTPUT_ARRAYS; i++)
    output_array_types[i] = CCTK_VARIABLE_REAL;

  void *output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) particle_psi_bssn,
        (void *) particle_alp };

  CCTK_INT num_arrays = NUM_INPUT_ARRAYS;

  int ierr = CCTK_InterpGridArrays(cctkGH,
                                   3,                   // number of dimensions
                                   operator_handle,
                                   param_table_handle,
                                   coord_system_handle,
                                   interp_num_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                   num_arrays,          // Number of input arrays
                                   (const CCTK_INT *)(&input_array_indices[0]),
                                   num_arrays,          // Number of output arrays
                                   (const CCTK_INT *)(&output_array_types[0]),
                                   &output_arrays[0]);

  if(ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT, "Could not destroy table");
}
