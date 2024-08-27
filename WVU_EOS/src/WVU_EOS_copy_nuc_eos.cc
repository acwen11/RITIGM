// Thorn      : WVU_EOS
// File       : WVU_EOS_Tabulated_known_T.cc
// Author(s)  : Allen Wen
// Description: This file copies EOS table constants

#include "cctk.h"
#include "WVU_EOS_Tabulated_headers.hh"

extern "C"
void WVU_EOS_share_nuc_eos_impl( CCTK_INT* out_nrho,
																 CCTK_INT* out_nT,
																 CCTK_INT* out_nYe,
																 CCTK_REAL* out_rhomin,
																 CCTK_REAL* out_rhomax,
																 CCTK_REAL* out_Tmin,
																 CCTK_REAL* out_Tmax,
																 CCTK_REAL* out_Yemin,
																 CCTK_REAL* out_Yemax,
																 CCTK_REAL* out_drho,
																 CCTK_REAL* out_dT,
																 CCTK_REAL* out_dYe) {

	using namespace nuc_eos;
	using namespace nuc_eos_private;

	*out_nrho = nrho;
	*out_nT = ntemp;
	*out_nYe = nye;
	*out_rhomin = eos_rhomin;
	*out_rhomax = eos_rhomax;
	*out_Tmin = eos_tempmin;
	*out_Tmax = eos_tempmax;
	*out_Yemin = eos_yemin;
	*out_Yemax = eos_yemax;
	*out_drho = drho;
	*out_dT = dtemp;
	*out_dYe = dye;
}
