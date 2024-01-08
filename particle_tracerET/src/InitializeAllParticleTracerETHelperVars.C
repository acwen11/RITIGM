#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// Split a comma-separated string of integers
// into an array of integers.
static int split_comma_separated_string_convert_to_ints(char* string, int** arr) {
  if( string == NULL ) return 0;

  // Count number of fields
  int nf = 1;
  for(char* ptr=string;*ptr!='\0';ptr++) if(*ptr==',') nf++;

  int i = 0;
  int* arr_l = (int*)malloc(sizeof(int)*nf);
  char* token  = strtok(string, ",");
  while( token ) {
    arr_l[i++] = atoi(token);
    token = strtok(NULL, ",");
  }
  *arr = arr_l;
  return nf;
}

static int count_duplicates_in_str(const char* string) {

  // If empty string return 0
  if( string[0] == '\0' ) return 0;

  // Sanitize the input by removing the '[' at the
  // beginning and the ']' at the end of the string
  char str[strlen(string)-2];
  memcpy(str, string+1, sizeof(char)*(strlen(string)-2));

  int* arr = NULL;
  int n = split_comma_separated_string_convert_to_ints(str, &arr);

  // Count the duplicates
  int dup = 0;
  for(int i=0;i<n;i++)
    for(int j=i+1;j<n;j++)
      if(arr[i] == arr[j]) {
        dup++;
        break;
      }

  free(arr);
  return dup;
}

/* This routine should only be called at cctk_iteration=start_tracing_particles_iteration */
void InitializeAllParticleTracerETHelperVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *RK4IterationCounter = 1;

  // Get number of maximum and active refinement levels
  const int mrl = GetMaxRefinementLevels(cctkGH);
  const int arl = GetRefinementLevels(cctkGH);
  *initial_number_of_active_refinement_levels = arl;

  // Check if user input is sane
  const int update_freq_min_granularity = 1<<(mrl-arl);
  const int output_freq_min_granularity = 2*update_freq_min_granularity;
  if( update_RK4_freq%update_freq_min_granularity != 0 )
    CCTK_VERROR("update_RK4_freq must be a multiple of %d, but got %d", update_freq_min_granularity, update_RK4_freq);

  if( output_freq == -3 ) {
    // Default value: 2*update_RK4_Freq
    *file_output_freq = 2*update_RK4_freq;
  }
  else if( output_freq == -1 || output_freq == 0 ) {
    // Disabled (only makes sense if update_RK4_freq=0
    if( update_RK4_freq != 0 )
      CCTK_ERROR("Particles' positions will be updated but no output requested. Please adjust the particle_tracerET::update_RK4_freq or the particle_tracerET::output_freq parameters.");
    *file_output_freq = 0;
  }
  else {
    if( output_freq == -2 ) {
      // Use the output frequency set by IO::out_every
      *file_output_freq = out_every;
      if( (*file_output_freq)%update_freq_min_granularity != 0 )
        CCTK_VERROR("output_freq must be a multiple of %d, but got %d (from IO::out_every)", output_freq_min_granularity, *file_output_freq);
    }
    else {
      // Custom user input
      *file_output_freq = output_freq;
      if( ((*file_output_freq)%update_freq_min_granularity != 0) )
        CCTK_VERROR("output_freq must be a multiple of %d, but got %d (from user input)", output_freq_min_granularity, *file_output_freq);
    }
  }

  // Further sanitize the output frequency
  if( (*file_output_freq) < 2*update_RK4_freq ) {
    CCTK_VWARN(CCTK_WARN_ALERT,"The requested output frequency (%d) is too small and will result in duplicated outputs.",*file_output_freq);
    *file_output_freq = 2*update_RK4_freq;
    CCTK_VWARN(CCTK_WARN_ALERT,"Adjusting output frequency to smallest value that does not result in duplicates (%d).",*file_output_freq);
  }

  // Set the time step based on dt at the *finest possible*
  // refinement level, even if it is not currently active.
  // Notice the factor of two is needed even when the update
  // frequency is set to one because the positions of the
  // particles are only updated every other time step.
  const char** time_refinement_factors = (const char**)CCTK_ParameterGet("time_refinement_factors", "Carpet", NULL);
  const int duplicates = count_duplicates_in_str(*time_refinement_factors);
  const CCTK_REAL dtmin = cctk_delta_time/(1<<(mrl-1-duplicates));
  *RK4_delta_time = 2*update_RK4_freq*dtmin;

  CCTK_VINFO("Initialized RK4IterationCounter to %d",*RK4IterationCounter);
  CCTK_VINFO("Number of requested refinement levels          : %d", mrl);
  CCTK_VINFO("Number of *active*  refinement levels          : %d", arl);
  CCTK_VINFO("Number of duplicates in time_refinement_factors: %d", duplicates);
  CCTK_VINFO("Requested particle position update every %d iterations of the finest refinement level (time step %e)", update_RK4_freq, update_RK4_freq*dtmin);
  CCTK_VINFO("Particle position is only updated every other iteration, so time step set to: %e", *RK4_delta_time);
  CCTK_VINFO("Particle positions will be written to file every %d iterations", *file_output_freq);

}
