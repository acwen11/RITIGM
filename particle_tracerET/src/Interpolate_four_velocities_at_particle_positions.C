#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "util_Table.h"
#include "util_String.h"

// This function is used to interpolate the particles'
// four-velocities to the particles' positions.
//
// We will interpolate a total of 8 quantities:
// u^{mu} and u_{mu}.
#define NUM_INPUT_ARRAYS  8
#define NUM_OUTPUT_ARRAYS 8

void Interpolate_four_velocities_at_particle_positions(
      cGH *cctkGH,
      int interp_num_points,
      double *particle_x,
      double *particle_y,
      double *particle_z,
      double *particle_u4U0,
      double *particle_u4U1,
      double *particle_u4U2,
      double *particle_u4U3,
      double *particle_u4D0,
      double *particle_u4D1,
      double *particle_u4D2,
      double *particle_u4D3 ) {

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
    = { "particle_tracerET::u4U0GF",
	"particle_tracerET::u4U1GF",
	"particle_tracerET::u4U2GF",
        "particle_tracerET::u4U3GF",
        "particle_tracerET::u4D0GF",
	"particle_tracerET::u4D1GF",
	"particle_tracerET::u4D2GF",
        "particle_tracerET::u4D3GF" };

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
    = { (void *) particle_u4U0,
	(void *) particle_u4U1,
	(void *) particle_u4U2,
        (void *) particle_u4U3,
        (void *) particle_u4D0,
	(void *) particle_u4D1,
	(void *) particle_u4D2,
        (void *) particle_u4D3 };

  CCTK_INT shift      = 0;
  CCTK_INT num_arrays = NUM_INPUT_ARRAYS;
  if( !output_four_velocity_u4U ) {
    shift      = 4;
    num_arrays = 4;
  }
  if( !output_four_velocity_u4D )
    num_arrays = 4;

  int ierr = CCTK_InterpGridArrays(cctkGH,
                                   3,                   // number of dimensions
                                   operator_handle,
                                   param_table_handle,
                                   coord_system_handle,
                                   interp_num_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                   num_arrays,          // Number of input arrays
                                   (const CCTK_INT *)(&input_array_indices[0]+shift),
                                   num_arrays,          // Number of output arrays
                                   (const CCTK_INT *)(&output_array_types[0]+shift),
                                   &output_arrays[0]+shift);

  if(ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT, "Could not destroy table");
}
