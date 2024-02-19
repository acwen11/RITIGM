#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"
#include "util_String.h"


#define NUM_INPUT_ARRAYS  3
#define NUM_OUTPUT_ARRAYS 3


void Interpolate_velocities_at_particle_positions(cGH *cctkGH,int interp_num_points,double *particle_x,double *particle_y,double *particle_z, double *particle_velx,double *particle_vely,double *particle_velz) {
  DECLARE_CCTK_PARAMETERS;

  for(int i=0;i<interp_num_points;i++) {
    if(fabs(particle_x[i])>out_of_bounds_xyz ||
       fabs(particle_y[i])>out_of_bounds_xyz ||
       fabs(particle_z[i])>out_of_bounds_xyz) {

      printf("Interpolate_velocities_at_particle_positions: FOUND A PARTICLE %d WAS OUT OF BOUNDS.\n",i);

      particle_x[i] = 1e-100;
      particle_y[i] = 1e-100;
      particle_z[i] = 1e-100;
    }
  }


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

  CCTK_INT operand_indices[NUM_INPUT_ARRAYS]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i=0; i<NUM_INPUT_ARRAYS; i++) operand_indices[i] = i;
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operand_indices, "operand_indices");

  CCTK_INT operation_codes[NUM_INPUT_ARRAYS]; //NUM_OUTPUT_ARRAYS + MAX_NUMBER_EXTRAS];
  for(int i=0; i<NUM_INPUT_ARRAYS; i++) operation_codes[i] = 0;
  Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                        operation_codes, "operation_codes");

  const void* interp_coords[NUM_INPUT_ARRAYS]
    = { (const void *) particle_x,
        (const void *) particle_y,
        (const void *) particle_z };

  // 3d input arrays
  CCTK_STRING input_array_names[NUM_INPUT_ARRAYS]
    = { "particle_tracerET::MHDvx",
	"particle_tracerET::MHDvy",
	"particle_tracerET::MHDvz" };

  CCTK_INT input_array_indices[NUM_INPUT_ARRAYS];
  for(int i=0; i<NUM_INPUT_ARRAYS; i++) {
    input_array_indices[i] = CCTK_VarIndex(input_array_names[i]);
    if(input_array_indices[i] < 0)
      CCTK_VWARN(CCTK_WARN_ABORT, "COULD NOT FIND VARIABLE '%s'.", input_array_names[i]);
  }

  CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS];
  for(int i=0; i<NUM_OUTPUT_ARRAYS; i++)
    output_array_types[i] = CCTK_VARIABLE_REAL;

  void * output_arrays[NUM_OUTPUT_ARRAYS]
    = { (void *) particle_velx,
	(void *) particle_vely,
	(void *) particle_velz };

  // actual interpolation call
  int ierr = CCTK_InterpGridArrays(cctkGH,
                                   3, // number of dimensions
                                   operator_handle,
                                   param_table_handle,
                                   coord_system_handle,
                                   interp_num_points,
                                   CCTK_VARIABLE_REAL,
                                   interp_coords,
                                   NUM_INPUT_ARRAYS, // Number of input arrays
                                   input_array_indices,
                                   NUM_OUTPUT_ARRAYS, // Number of output arrays
                                   output_array_types,
                                   output_arrays);
  if(ierr < 0) {
    CCTK_WARN(CCTK_WARN_ALERT, "interpolation screwed up");
    Util_TableDestroy(param_table_handle);
    exit(1);
  }

  ierr = Util_TableDestroy(param_table_handle);
  if(ierr != 0) CCTK_WARN(CCTK_WARN_ABORT, "Could not destroy table");
}
