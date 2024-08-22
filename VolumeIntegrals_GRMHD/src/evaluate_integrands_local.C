#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <string.h>

#include "integrands.C"

void VI_GRMHD_ComputeIntegrand(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // rhs: m_neutron_MeV * MeV_to_erg * PRESSGF * LENGTHGF^3
  // TODO: Get from Margherita, make also Margherita baryon_mass consistent with 
  //       value used to generate .h5 table -> Put in param.ccl!
  constexpr double my_baryon_mass = 939.565379 * 1.60217733e-6 * 1.80123683248503e-39 * (6.77269222552442e-06*6.77269222552442e-06*6.77269222552442e-06);

  int which_integral = NumIntegrals - *IntegralCounter + 1;

  /* Note: Must extend this if/else statement if adding a new integrand! */
  if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"centerofmass")) { 
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	  
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  CoM_integrand(VolIntegrand1,VolIntegrand2,VolIntegrand3,VolIntegrand4, index,w_lorentz, rho, gxx,gxy,gxz, gyy,gyz,gzz,x,y,z); 
    
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"coordvolume")) { 
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	  
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  CoordVol_integrand(VolIntegrand1, index, rho, dens_atmo); 
    
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"restmass")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	    
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  M0_integrand(VolIntegrand1, index,w_lorentz, rho, gxx,gxy,gxz, gyy,gyz,gzz);
   
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"usepreviousintegrands")) {
   
	  /* Do Nothing; the action for this is below. */

  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"one")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	    
	    int index = CCTK_GFINDEX3D(cctkGH,i,j,k); 
	    VolIntegrand1[index]=1.0;
   
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"density_weighted_norm_B_field")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k); 
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    
	    
	mean_density_weighted_B(VolIntegrand1,VolIntegrand2,index,
		          indexX,indexY,indexZ,
                          w_lorentz,rho,
			  Bvec,
                          gxx,gxy,gxz,gyy,gyz,gzz);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"kinetic_energy")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);

	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;
	    
	kinetic(VolIntegrand1,VolIntegrand2,VolIntegrand3,
		 index,
                 indexX,indexY,indexZ,
                 vel,w_lorentz,rho,eps,press,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"kinetic_energy_palenzuela")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);

	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	kinetic_palenzuela(VolIntegrand1,VolIntegrand2,VolIntegrand3,
		 index,
                 indexX,indexY,indexZ,
                 vel,w_lorentz,rho,eps,press,
		 alp,betax,betay,betaz,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"kinetic_energy_shibata")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);

	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);
	    
	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	kinetic_shibata(VolIntegrand1,VolIntegrand2,VolIntegrand3,
		 index,
                 indexX,indexY,indexZ,
                 vel,w_lorentz,rho,eps,press,
		 alp,betax,betay,betaz,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"kinetic_energy_total")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);

	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);
	    
	kinetic_tot(VolIntegrand1,index,
                 indexX,indexY,indexZ,
                 vel,w_lorentz,rho,eps,press,
                 gxx,gxy,gxz,gyy,gyz,gzz);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"thermal_energy")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	    
	thermal(VolIntegrand1,VolIntegrand2,VolIntegrand3,
			  index,
                          w_lorentz,rho,
			  eps,entropy,
                          gxx,gxy,gxz,gyy,gyz,gzz,
			  my_baryon_mass);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"magnetic_energy_total")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	magnetic_tot(VolIntegrand1,VolIntegrand2,index,
                 indexX,indexY,indexZ,
                 vel,Bvec,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"em_energy_ab")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	magnetic_tot_12(VolIntegrand1,VolIntegrand2,index,
                 indexX,indexY,indexZ,
                 rho,dens_a,dens_b,
		 vel,Bvec,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"em_energy_cd")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	magnetic_tot_12(VolIntegrand1,VolIntegrand2,index,
                 indexX,indexY,indexZ,
                 rho,dens_c,dens_d,
		 vel,Bvec,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);;
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"volume_average_norm_B_ab")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	volume_norm_B_12(VolIntegrand1,VolIntegrand2,VolIntegrand3,VolIntegrand4,
		 index,
                 indexX,indexY,indexZ,
		 rho,dens_a,dens_b,
                 Bvec,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"volume_average_norm_B_cd")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    
	int indexX = CCTK_GFINDEX4D(cctkGH,i,j,k,0);
        int indexY = CCTK_GFINDEX4D(cctkGH,i,j,k,1);
        int indexZ = CCTK_GFINDEX4D(cctkGH,i,j,k,2);    

	// Get BNS COM from BNSTrackerGen
	// Make sure we have updated values 
	// before this integral is called!

        double cms_x = *bns_cms_x;
        double cms_y = *bns_cms_y;

	volume_norm_B_12(VolIntegrand1,VolIntegrand2,VolIntegrand3,VolIntegrand4,
		 index,
                 indexX,indexY,indexZ,
		 rho,dens_c,dens_d,
                 Bvec,
                 gxx,gxy,gxz,gyy,gyz,gzz,
		 x,y,z,cms_x,cms_y);
    }
  } else if(CCTK_EQUALS(Integration_quantity_keyword[which_integral],"magnetic_energy_comov")) {
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { 
	
	int index  = CCTK_GFINDEX3D(cctkGH,i,j,k);    

	magnetic_co(VolIntegrand1,index,
                 smallb2,w_lorentz,
                 gxx,gxy,gxz,gyy,gyz,gzz);
    }
  } else {
    /* Print a warning if no integrand is computed because Integration_quantity_keyword unrecognized. */
    printf("VolumeIntegrals: WARNING: Integrand not computed. Did not understand Integration_quantity_keyword[%d] = %s\n",which_integral,Integration_quantity_keyword[which_integral]);
  }

  if(cctk_iteration==0) {
    volintegral_inside_sphere__center_x[which_integral]=volintegral_sphere__center_x_initial[which_integral];
    volintegral_inside_sphere__center_y[which_integral]=volintegral_sphere__center_y_initial[which_integral];
    volintegral_inside_sphere__center_z[which_integral]=volintegral_sphere__center_z_initial[which_integral];

    volintegral_outside_sphere__center_x[which_integral]=volintegral_sphere__center_x_initial[which_integral];
    volintegral_outside_sphere__center_y[which_integral]=volintegral_sphere__center_y_initial[which_integral];
    volintegral_outside_sphere__center_z[which_integral]=volintegral_sphere__center_z_initial[which_integral];
  }

  if(volintegral_sphere__tracks__amr_centre[which_integral]!=-1) {
    int which_centre = volintegral_sphere__tracks__amr_centre[which_integral];

    volintegral_inside_sphere__center_x[which_integral] = position_x[which_centre];
    volintegral_inside_sphere__center_y[which_integral] = position_y[which_centre];
    volintegral_inside_sphere__center_z[which_integral] = position_z[which_centre];

    volintegral_outside_sphere__center_x[which_integral] = position_x[which_centre];
    volintegral_outside_sphere__center_y[which_integral] = position_y[which_centre];
    volintegral_outside_sphere__center_z[which_integral] = position_z[which_centre];
  }

  /* ZERO OUT INTEGRATION REGIONS */

  /* The below code supports zeroing out of arbitrary spherical shells. 
     In the case of integration INSIDE a full sphere, this code also supports
     moving spheres if AMR centre tracking is enabled. This can be used to 
     track compact objects. 

     Here's one way:
     Generally the lapse at the center of these objects is minimized. 
     If the object has a length scale of R and is known to be centered at x,y,z, compute 
     X^i = Integral [ (1-lapse) x^i dV] / Integral [(1-lapse) dV]
     over the spherical volume centered at x,y,z with radius R. 
     This should yield X,Y,Z = x,y,z to good approximation. 
     If you enable AMR centre tracking, it should work just as well
     as any other method for tracking compact objects, if not better. */

  /* Set integrands to zero outside a sphere centered at x,y,z. 
     I.e., this results in the integral being restricted INSIDE sphere */
  if(volintegral_inside_sphere__radius[which_integral]>0.0) {
    double radius=volintegral_inside_sphere__radius[which_integral];
    double xprime=volintegral_sphere__center_x_initial[which_integral];
    double yprime=volintegral_sphere__center_y_initial[which_integral];
    double zprime=volintegral_sphere__center_z_initial[which_integral];
    if(cctk_iteration>0 && (amr_centre__tracks__volintegral_inside_sphere[which_integral]!=-1 || volintegral_sphere__tracks__amr_centre[which_integral]!=-1) ) {
      xprime=volintegral_inside_sphere__center_x[which_integral];
      yprime=volintegral_inside_sphere__center_y[which_integral];
      zprime=volintegral_inside_sphere__center_z[which_integral];
    }
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double x_minus_xprime = x[index] - xprime;
	  double y_minus_xprime = y[index] - yprime;
	  double z_minus_xprime = z[index] - zprime;
	  if(sqrt(x_minus_xprime*x_minus_xprime + y_minus_xprime*y_minus_xprime + z_minus_xprime*z_minus_xprime)>radius) 
	    VolIntegrand1[index]=VolIntegrand2[index]=VolIntegrand3[index]=VolIntegrand4[index]=0.0;
	}
  }

  /* Set integrands to zero inside a sphere centered at x,y,z. 
     I.e., this results in the integral being restricted OUTSIDE sphere.
     Combine this with above to get spherical shell. 
  
     Note that volume integrals outside a sphere are fixed at the 
     original sphere center position for all time, unlike volume 
     integrals inside a sphere. We do this because the latter are 
     generally used for tracking moving compact objects or other things. */
  if(volintegral_outside_sphere__radius[which_integral]>0.0) {
    double radius=volintegral_outside_sphere__radius[which_integral];
    double xprime=volintegral_outside_sphere__center_x[which_integral];
    double yprime=volintegral_outside_sphere__center_y[which_integral];
    double zprime=volintegral_outside_sphere__center_z[which_integral];
#pragma omp parallel for
    for (int k=0;k<cctk_lsh[2];k++) for (int j=0;j<cctk_lsh[1];j++) for (int i=0;i<cctk_lsh[0];i++) { int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double x_minus_xprime = x[index] - xprime;
	  double y_minus_xprime = y[index] - yprime;
	  double z_minus_xprime = z[index] - zprime;
	  if(sqrt(x_minus_xprime*x_minus_xprime + y_minus_xprime*y_minus_xprime + z_minus_xprime*z_minus_xprime)<=radius) 
	    VolIntegrand1[index]=VolIntegrand2[index]=VolIntegrand3[index]=VolIntegrand4[index]=0.0;
	}
  }
}
