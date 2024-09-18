#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <math.h>

extern "C" void GridMover_Move_Grids(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

	position_x[index_region] = path_radius * cos(path_freq * cctk_time);
	position_y[index_region] = path_radius * sin(path_freq * cctk_time);
	position_z[index_region] = 0;
}
