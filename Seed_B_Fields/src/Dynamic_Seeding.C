#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void Evaluate_SeedingVars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration%seed_every!=0) {

	 *EvalMyScalar = 0;

  } else {
	  
	  if (use_separation && (*bns_sep_tot < seeding_separation)) {

		  *EvalMyScalar = 1;

		  CCTK_INFO ("##################");
			CCTK_INFO ("Seeding magnetic field");
			CCTK_INFO ("BNS separation is below threshold");
			CCTK_INFO ("##################");
	  
	  } else if (use_time && (cctk_time > seeding_time)) {

		  *EvalMyScalar = 1;

		  CCTK_INFO ("##################");
			CCTK_INFO ("Seeding magnetic field");
			CCTK_INFO ("Evolution time is below threshold");
			CCTK_INFO ("##################");

	  } else {

		  *EvalMyScalar = 0;

		  CCTK_INFO ("##################");
			CCTK_INFO ("Thresholds not met, skipping");
			CCTK_INFO ("##################");

	  }

  }

*EvalMyScalar = *EvalMyScalar && !(*DoneSeeding);

}

void Init_Seeding_Vars(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *EvalMyScalar = 0;
  *DoneSeeding = 0;

}

void Done_Seeding_Var(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *DoneSeeding = 1;
  CCTK_INFO ("##################");
  CCTK_INFO ("Seeding magnetic field finished");
  CCTK_INFO ("done_seeding was set to 1");
  CCTK_INFO ("##################");

}




