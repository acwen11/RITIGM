#include <cstdio>
#include <cassert>
#include <vector>
#include <ios>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <algorithm>

using std::min;
using std::max;

extern "C"
void Seed_B_Fields_Func (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO ("Seeding B fields");

  double dX = CCTK_DELTA_SPACE(0);
  double dY = CCTK_DELTA_SPACE(1);
  double dZ = CCTK_DELTA_SPACE(2);
  double Ax_yshift_staggering = 0.;
  double Ay_xshift_staggering = 0.;

  if(enable_IllinoisGRMHD_staggered_A_fields) {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          int indexip1=CCTK_GFINDEX3D(cctkGH,min(i+1,cctk_lsh[0]-1),j,k);
          int indexjp1=CCTK_GFINDEX3D(cctkGH,i,min(j+1,cctk_lsh[1]-1),k);
          int indexkp1=CCTK_GFINDEX3D(cctkGH,i,j,min(k+1,cctk_lsh[2]-1));
          int indexip1kp1=CCTK_GFINDEX3D(cctkGH,min(i+1,cctk_lsh[0]-1),j,min(k+1,cctk_lsh[2]-1));
          int indexjp1kp1=CCTK_GFINDEX3D(cctkGH,i,min(j+1,cctk_lsh[1]-1),min(k+1,cctk_lsh[2]-1));

	  // Get BNS locations from BNSTrackerGen
	  // Difference to Seed_Magnetic_Fields.C:
	  //   In case of BNS we calculate the locations
	  //   by using BNSTrackerGen 

          double loc_x1 = *bns_x_1;
          double loc_x2 = *bns_x_2;
	  double loc_y1 = *bns_y_1;
	  double loc_y2 = *bns_y_2;

          double xL = x[index];
          double yL = y[index];
          double zL = z[index];
          double x1 = xL - loc_x1;
          double x2 = xL - loc_x2;
	  double y1 = yL - loc_y1;
	  double y2 = yL - loc_y2;

	  double r0 = sqrt(xL*xL + yL*yL + zL*zL);

	  double r_at_Ax_stagger = sqrt(xL*xL + (yL+0.5*dY)*(yL+0.5*dY) + (zL+0.5*dZ)*(zL+0.5*dZ));
          double r_at_Ay_stagger = sqrt((xL+0.5*dX)*(xL+0.5*dX) + yL*yL + (zL+0.5*dZ)*(zL+0.5*dZ));

          double r1 = sqrt(x1*x1 + y1*y1 + zL*zL);
          double r2 = sqrt(x2*x2 + y2*y2 + zL*zL);

	  double r1_at_Ax_stagger = sqrt( x1*x1 + (y1+0.5*dY)*(y1+0.5*dY) + (zL+0.5*dZ)*(zL+0.5*dZ) );
	  double r1_at_Ay_stagger = sqrt( (x1+0.5*dX)*(x1+0.5*dX) + y1*y1 + (zL+0.5*dZ)*(zL+0.5*dZ) );

	  double r2_at_Ax_stagger = sqrt( x2*x2 + (y2+0.5*dY)*(y2+0.5*dY) + (zL+0.5*dZ)*(zL+0.5*dZ) );
	  double r2_at_Ay_stagger = sqrt( (x2+0.5*dX)*(x2+0.5*dX) + y2*y2 + (zL+0.5*dZ)*(zL+0.5*dZ) );

	  if(CCTK_EQUALS(A_field_type, "poloidal_A_interior"))
	    {
	      double PL = press[index]; // Assumes HydroBase pressure is set!

	      double PLip1 = press[indexip1];
	      double PLjp1 = press[indexjp1];
	      double PLkp1 = press[indexkp1];
	      double PLip1kp1 = press[indexip1kp1];
	      double PLjp1kp1 = press[indexjp1kp1];
          
	      double Pressure_at_Ax_stagger = 0.25*(PL + PLjp1 + PLkp1 + PLjp1kp1); 

	      double Pressure_at_Ay_stagger = 0.25*(PL + PLip1 + PLkp1 + PLip1kp1);

	      double A_b_x  = A_b * exp( - g_width * pow(r_at_Ax_stagger-g_r0,2.0) );
	      double A_b_y  = A_b * exp( - g_width * pow(r_at_Ay_stagger-g_r0,2.0) );

	      double A_b_x1 = A_b * exp( - g_width * pow(r1_at_Ax_stagger-g_r0,2.0) );
	      double A_b_y1 = A_b * exp( - g_width * pow(r1_at_Ay_stagger-g_r0,2.0) );
	      double A_b_x2 = A_b * exp( - g_width * pow(r2_at_Ax_stagger-g_r0,2.0) );
	      double A_b_y2 = A_b * exp( - g_width * pow(r2_at_Ay_stagger-g_r0,2.0) );
 
	      if(!have_two_NSs) {       	
		if( r0<=r_NS1 && r0>=r_NS1_low ) 
			  {
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -(yL + 0.5*dY)*A_b_x*pow(max(Pressure_at_Ax_stagger-P_cut,0.0),n_s);
		    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (xL + 0.5*dX)*A_b_y*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
		  }
		if( r0>r_NS1 ) 
		  {	
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No external B-field.
                    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
		  }
		if( r0<r_NS1_low ) 
		  {	
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No internal field if desired
                    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
		  }
	      }
	      else {
		if(r1<=r_NS1 && r1>=r_NS1_low ) 
		  {
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (x1 + 0.5*dX)*A_b_y1*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
		    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -(y1 + 0.5*dY)*A_b_x1*pow(max(Pressure_at_Ax_stagger-P_cut,0.0),n_s);
		  }
		if(r2<=r_NS2 && r2>=r_NS2_low ) 
	          {
		    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  (x2 + 0.5*dX)*A_b_y2*pow(max(Pressure_at_Ay_stagger-P_cut,0.0),n_s);
		    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -(y2 + 0.5*dY)*A_b_x2*pow(max(Pressure_at_Ax_stagger-P_cut,0.0),n_s);
		  }
		if(r1>r_NS1 && r2>r_NS2) 
		  {	
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No external B-field.
                    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
		  }
		if(r1<r_NS1_low || r2<r_NS2_low) 
		  {	
	            Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No internal field if desired
                    Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
		  }
	      }
	    }
	  // Modified to include poloidal vector potential option
	  else if(CCTK_EQUALS(A_field_type, "dipolar_A_everywhere"))
	    {

              CCTK_WARN(1, "dipolar_A_everywhere is currently not available for dynamic seeding!");

	      /*

	      if(!enable_IllinoisGRMHD_staggered_A_fields)
		{
		  CCTK_WARN(1, "dipolar_A_everywhere requires a staggered grid");
		}

	      double pi = 3.141592653589793;

	      // \varpi^2
	      double varpi1_2 = x1 * x1 + yL * yL;
	      double varpi2_2 = x2 * x2 + yL * yL;

	      // Current loop radius
	      double r01 = r_zero_NS1 * r_NS1;
	      double r02 = r_zero_NS2 * r_NS2;

	      // phi-component of vector potential in spherical basis
	      // Eq (2) in Paschalidis et al PRD 88 021504(R) (2013)
	      double Ap1 = pi * r01 * r01 * I_zero_NS1 * pow(r01 * r01 + r1 * r1, -1.5) * (1.0 + 1.875 * r01 * r01 * (r01 * r01 + varpi1_2) * pow(r01 * r01 + r1 * r1, -2.0));
	      double Ap2 = pi * r02 * r02 * I_zero_NS2 * pow(r02 * r02 + r2 * r2, -1.5) * (1.0 + 1.875 * r02 * r02 * (r02 * r02 + varpi2_2) * pow(r02 * r02 + r2 * r2, -2.0));

	      // x-component of vector potential in Cartesian basis
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = - (yL + 0.5 * dY) * (Ap1 + Ap2);

	      // y-component of vector potential in Cartesian basis
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = (x1 + 0.5 * dX) * Ap1 + (x2 + 0.5 * dX) * Ap2;

	      */
	    }

	  // z-component of vector potential in Cartesian basis
          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] =  0.0;
          Aphi[index]  = 0.0;

          // Copy to IllinoisGRMHD internal gfs
	  Ax[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	  Ay[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	  Az[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];
	  psi6phi[index] = Aphi[index];

	  // Now, copy to other timelevels, too
	  // We also do _p_p but because of subsequent timelevel
	  // rotation we probably just need _p

          psi6phi_p[index] = psi6phi[index];
          Ax_p[index] = Ax[index]; 
          Ay_p[index] = Ay[index];
          Az_p[index] = Az[index];

          psi6phi_p_p[index] = psi6phi[index];
          Ax_p_p[index] = Ax[index]; 
          Ay_p_p[index] = Ay[index];
          Az_p_p[index] = Az[index];
          
        }

  } else {
#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {

          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);
          double PL = press[index]; // Assumes HydroBase pressure is set!

          double r_at_A = sqrt(x[index]*x[index] + y[index]*y[index] + z[index]*z[index]);

          double A_b_r  = A_b * exp( - g_width * pow(r_at_A-g_r0,2.0) );

	  if( r_at_A<=r_NS1 && r_at_A>=r_NS1_low ) 
	    {
              Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = -y[index]*A_b_r*pow(max(PL-P_cut,0.0),n_s);
              Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] =  x[index]*A_b_r*pow(max(PL-P_cut,0.0),n_s);
	    }
	  if( r_at_A>r_NS1 ) 
            {	
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No external B-field.
              Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
            }
          if( r_at_A<r_NS1_low ) 
	    {	
	      Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)] = 0.0; // No internal field if desired
              Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)] = 0.0;
	    }

          Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] =  0.0;
          Aphi[index]  = 0.0;

          // Copy to IllinoisGRMHD internal gfs
	  Ax[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	  Ay[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	  Az[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];
	  psi6phi[index] = Aphi[index];

	  // Now, copy to other timelevels, too
	  // We also do _p_p but because of subsequent timelevel
	  // rotation we probably just need _p

          psi6phi_p[index] = psi6phi[index];
          Ax_p[index] = Ax[index]; 
          Ay_p[index] = Ay[index];
          Az_p[index] = Az[index];

          psi6phi_p_p[index] = psi6phi[index];
          Ax_p_p[index] = Ax[index]; 
          Ay_p_p[index] = Ay[index];
          Az_p_p[index] = Az[index];
        }
  }
}

