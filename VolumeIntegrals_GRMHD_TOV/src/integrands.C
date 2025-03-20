#define _USE_MATH_DEFINES
#include <cmath>

#ifndef WVU_INTEGRANDS_C
#define WVU_INTEGRANDS_C

/*
  NRPy+ Python code that generates below C code:
  (Note: used commit 6a86d39aa3e4c863caa51c2fc7a04f3d06727167 in
         https://github.com/zachetienne/nrpytutorial)

import sympy as sp
import grid as gri
import indexedexp as ixp
from outputC import *

DIM=3

# Need to read from memory & write to memory manually due to speed limiter if() statement.
with open("compute_rhostar.h","w") as file:
    file.write("""
// This function computes
//  rho_* = alpha*u^0*sqrt(gamma)*rho_0
//        = w_lorentz*sqrt(gamma)*rho_0 ,
// where gamma = determinant of 3-metric,
// rho_0 is the rest-mass density, and
// w_lorentz is the Lorentz factor\n\n""")

    for i in range(DIM):
        for j in range(i,DIM):
            file.write("const CCTK_REAL gammaDD"+str(i)+str(j)+" = g"+chr(ord('x')+i)+chr(ord('x')+j)+"[index];\n")
    file.write("""
const CCTK_REAL w_lorentz = w_lorentz[index];
const CCTK_REAL rho0      = rhozero[  index];

CCTK_REAL w_lorentz_limited = w_lorentz;
if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;
""")

    rho0,w_lorentz_limited = sp.symbols("rho0 w_lorentz_limited")
    gammaDD                = ixp.declarerank2("gammaDD", "sym01",DIM=3)

    dummy, detgamma = ixp.symm_matrix_inverter3x3(gammaDD)

    rhostar = w_lorentz_limited*sp.sqrt(detgamma)*rho0

    rho_star_str_ugly = outputC(rhostar, "const CCTK_REAL rhostar", filename="returnstring",params="includebraces=False")
    # Beautify rho_star_str_ugly
    rho_star_str = rho_star_str_ugly.replace(
        "pow(gammaDD02, 2)","(gammaDD02*gammaDD02)").replace(
        "pow(gammaDD12, 2)","(gammaDD12*gammaDD12)").replace(
        "pow(gammaDD01, 2)","(gammaDD01*gammaDD01)")

    file.write(rho_star_str)
*/

inline CCTK_REAL compute_rho_star(const int index, const CCTK_REAL *restrict w_lorentzGF, const CCTK_REAL *restrict rho0GF,
                                  const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                                  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  DECLARE_CCTK_PARAMETERS;

  // This function computes
  //  rho_* = alpha*u^0*sqrt(gamma)*rho_0
  //        = w_lorentz*sqrt(gamma)*rho_0 ,
  // where gamma = determinant of 3-metric,
  // rho_0 is the rest-mass density, and
  // w_lorentz is the Lorentz factor

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL w_lorentz = w_lorentzGF[index];
  const CCTK_REAL rho0      = rho0GF[     index];

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;
  /*
   *  Original SymPy expression:
   *  "rhostar = rho0*w_lorentz_limited*sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*gammaDD12**2 - gammaDD01**2*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - gammaDD02**2*gammaDD11)"
   */
  const CCTK_REAL rhostar = rho0*w_lorentz_limited*sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*(gammaDD12*gammaDD12) - (gammaDD01*gammaDD01)*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - (gammaDD02*gammaDD02)*gammaDD11);

  return rhostar;
}

inline CCTK_REAL compute_sqrtgamma(const int index,
                                   const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                                   const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL sqrtgamma = sqrt(gammaDD00*gammaDD11*gammaDD22 - gammaDD00*(gammaDD12*gammaDD12) - (gammaDD01*gammaDD01)*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - (gammaDD02*gammaDD02)*gammaDD11);

  return sqrtgamma;
}


/* Center of Mass: */
inline void CoM_integrand(double *VolIntegrand1,double *VolIntegrand2,double *VolIntegrand3,double *VolIntegrand4, const int index,
                          const CCTK_REAL *restrict w_lorentz, const CCTK_REAL *restrict rho0,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
                          const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z) {
  double rho_starL = compute_rho_star(index, w_lorentz, rho0, gxx,gxy,gxz, gyy,gyz,gzz);
  VolIntegrand1[index] = rho_starL*x[index];
  VolIntegrand2[index] = rho_starL*y[index];
  VolIntegrand3[index] = rho_starL*z[index];
  VolIntegrand4[index] = rho_starL;
}

/* Rest Mass: */
inline void M0_integrand(double *VolIntegrand1, const int index,
                          const CCTK_REAL *restrict w_lorentz, const CCTK_REAL *restrict rho0,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  double rho_starL = compute_rho_star(index, w_lorentz, rho0, gxx,gxy,gxz, gyy,gyz,gzz);
  VolIntegrand1[index] = rho_starL;
}

/* Density weighted norm of magnetic field strength: */
inline void mean_density_weighted_B(double *VolIntegrand1,double *VolIntegrand2,const int index,
		          const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict w_lorentz, const CCTK_REAL *restrict rho0,
			  const CCTK_REAL *restrict Bvec,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL Bx        = Bvec[index_vecX];
  const CCTK_REAL By        = Bvec[index_vecY];
  const CCTK_REAL Bz        = Bvec[index_vecZ];
	
  double rho_starL = compute_rho_star(index, w_lorentz, rho0,gxx,gxy,gxz, gyy,gyz,gzz);

  double norm_B_sq =    gammaDD00*Bx*Bx+
	             2.*gammaDD01*Bx*By+
		     2.*gammaDD02*Bx*Bz+
                        gammaDD11*By*By+
	             2.*gammaDD12*By*Bz+
		        gammaDD22*Bz*Bz;
  VolIntegrand1[index] = rho_starL*sqrt(norm_B_sq);
  VolIntegrand2[index] = rho_starL;
}

/* Volume averaged norm of B in region (dens_1,dens_2]: */
inline void volume_norm_B_12(double *VolIntegrand1,double *VolIntegrand2, double *VolIntegrand3, double *VolIntegrand4,
		          const int index,
		          const int index_vecX,const int index_vecY, const int index_vecZ,
			  const CCTK_REAL *restrict rho0, const CCTK_REAL dens_1, const CCTK_REAL dens_2,
			  const CCTK_REAL *restrict Bvec,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {

  const CCTK_REAL my_rho = rho0[index];

  if ( my_rho > dens_1 && my_rho <= dens_2 ) {

	  const CCTK_REAL gammaDD00 = gxx[index];
	  const CCTK_REAL gammaDD01 = gxy[index];
	  const CCTK_REAL gammaDD02 = gxz[index];
	  const CCTK_REAL gammaDD11 = gyy[index];
	  const CCTK_REAL gammaDD12 = gyz[index];
	  const CCTK_REAL gammaDD22 = gzz[index];

	  const CCTK_REAL B_contra [3] {Bvec[index_vecX],Bvec[index_vecY],Bvec[index_vecZ]};

	  CCTK_REAL posx      = x[index]-cms_x;
	  CCTK_REAL posy      = y[index]-cms_y;
	  CCTK_REAL posz      = z[index];

	// Cylindrical coordinates and some 
	// quantities for transformation of vectors. 

	  const double posr2      = posx*posx+posy*posy;
	  const double posr       = sqrt(posr2);
 
	  double xy_over_r2   = posx*posy/posr2;
	  double x2_over_r2   = posx*posx/posr2;
	  double y2_over_r2   = posy*posy/posr2;

	  if ( posr <=  1e-15 ) {
  
		 xy_over_r2   = 0.0; 
		 x2_over_r2   = 0.5; 
		 y2_over_r2   = 0.5; 
	  }

	  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz, gyy,gyz,gzz);

	// Covariant B

	  const CCTK_REAL B_cov [3] {gammaDD00*B_contra[0]+gammaDD01*B_contra[1]+gammaDD02*B_contra[2],
                                    gammaDD01*B_contra[0]+gammaDD11*B_contra[1]+gammaDD12*B_contra[2],
                                    gammaDD02*B_contra[0]+gammaDD12*B_contra[1]+gammaDD22*B_contra[2]};

	  const CCTK_REAL B_sq =  B_cov[0]*B_contra[0] + B_cov[1]*B_contra[1] + B_cov[2]*B_contra[2];

	// Finally, compute B_phi*B^phi
 
	  const CCTK_REAL B_2_tor = B_cov[0]*B_contra[0]*y2_over_r2
		- B_cov[0]*B_contra[1]*xy_over_r2
		- B_cov[1]*B_contra[0]*xy_over_r2
		+ B_cov[1]*B_contra[1]*x2_over_r2;

	// Finally, B_r*B^r + B_z*B^z

	  const CCTK_REAL B_2_pol = B_cov[0]*B_contra[0]*x2_over_r2
		+ B_cov[0]*B_contra[1]*xy_over_r2
		+ B_cov[1]*B_contra[0]*xy_over_r2
		+ B_cov[1]*B_contra[1]*y2_over_r2
		+ B_cov[2]*B_contra[2];

	  VolIntegrand1[index] = sqrtgamma*sqrt(abs(B_sq));
	  VolIntegrand2[index] = sqrtgamma*sqrt(abs(B_2_tor));
	  VolIntegrand3[index] = sqrtgamma*sqrt(abs(B_2_pol));
	  VolIntegrand4[index] = sqrtgamma;

  } else {

	  VolIntegrand1[index] = 0.0;
	  VolIntegrand2[index] = 0.0;
	  VolIntegrand3[index] = 0.0;
	  VolIntegrand4[index] = 0.0;
  }

}

/* Kinetic energy as given by https://arxiv.org/pdf/1206.5911.pdf */
/* See also https://arxiv.org/pdf/2102.01346.pdf for rotational energy */
/* See also https://arxiv.org/abs/astro-ph/0605331 */
inline void kinetic_shibata(double *VolIntegrand1, double *VolIntegrand2, double *VolIntegrand3, const int index,
		    const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict velGF,
                          const CCTK_REAL *restrict w_lorentzGF, 
			  const CCTK_REAL *restrict rho0, const CCTK_REAL *restrict epsGF, 
			  const CCTK_REAL *restrict pressGF,
			  const CCTK_REAL *restrict alpGF, const CCTK_REAL *restrict betaxGF, const CCTK_REAL *restrict betayGF, const CCTK_REAL *restrict betazGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
  			  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {
  
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL my_rho    = rho0[index];
  const CCTK_REAL lfac      = w_lorentzGF[index];
  const CCTK_REAL my_press  = pressGF[index];
  const CCTK_REAL my_eps    = epsGF[index];

  const CCTK_REAL vx        = velGF[index_vecX];
  const CCTK_REAL vy        = velGF[index_vecY];
  const CCTK_REAL vz        = velGF[index_vecZ];

  const CCTK_REAL my_lapse  = alpGF[index];
  const CCTK_REAL my_shiftx = betaxGF[index];
  const CCTK_REAL my_shifty = betayGF[index];
  const CCTK_REAL my_shiftz = betazGF[index];

  const CCTK_REAL posx      = x[index]-cms_x;
  const CCTK_REAL posy      = y[index]-cms_y;
  const CCTK_REAL posz      = z[index];

  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz,gyy,gyz,gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if(lfac > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

// Covariant 3 velocity

  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz};

  const CCTK_REAL v_2 = v_cov[0]*vx + v_cov[1]*vy + v_cov[2]*vz;

  const CCTK_REAL vz_2 = v_cov[2]*vz;

// Scalar product with shift

  const CCTK_REAL shift_v = v_cov[0]*my_shiftx + v_cov[1]*my_shifty + v_cov[2]*my_shiftz;

  const CCTK_REAL shiftz_vz = v_cov[2]*my_shiftz;

// This is just the total kinetic energy

  VolIntegrand1[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(my_lapse*v_2-shift_v);

// Cylindrical coordinates and some 
// quantities for transformation of vectors. 

  const double posr2      = posx*posx+posy*posy;
  const double posr       = sqrt(posr2);
 
  double xy_over_r2   = posx*posy/posr2;
  double x2_over_r2   = posx*posx/posr2;
  double y2_over_r2   = posy*posy/posr2;

  if ( posr <=  1e-15 ) {
  
	 xy_over_r2   = 0.0; 
	 x2_over_r2   = 0.5; 
	 y2_over_r2   = 0.5; 
  }

// Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2 {0.0};

  v_phi2 = v_cov[0]*vx*y2_over_r2
	  - v_cov[0]*vy*xy_over_r2
	  - v_cov[1]*vx*xy_over_r2
	  + v_cov[1]*vy*x2_over_r2;

// Finally, compute v_phi*shift^phi

  CCTK_REAL v_phi_shift_phi {0.0};

  v_phi_shift_phi = v_cov[0]*my_shiftx*y2_over_r2
	  - v_cov[0]*my_shifty*xy_over_r2
	  - v_cov[1]*my_shiftx*xy_over_r2
	  + v_cov[1]*my_shifty*x2_over_r2;

// This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(my_lapse*v_phi2-v_phi_shift_phi);

// Kinetic energy in z-direction

  VolIntegrand3[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(my_lapse*vz_2-shiftz_vz);

}

/* Kinetic energy as given by https://arxiv.org/pdf/2112.08413.pdf */
inline void kinetic_palenzuela(double *VolIntegrand1, double *VolIntegrand2, double *VolIntegrand3, const int index,
		    const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict velGF,
                          const CCTK_REAL *restrict w_lorentzGF, 
			  const CCTK_REAL *restrict rho0, const CCTK_REAL *restrict epsGF, 
			  const CCTK_REAL *restrict pressGF,
			  const CCTK_REAL *restrict alpGF, const CCTK_REAL *restrict betaxGF, const CCTK_REAL *restrict betayGF, const CCTK_REAL *restrict betazGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
  			  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {
  
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL my_rho    = rho0[index];
  const CCTK_REAL lfac      = w_lorentzGF[index];
  const CCTK_REAL my_press  = pressGF[index];
  const CCTK_REAL my_eps    = epsGF[index];

  const CCTK_REAL vx        = velGF[index_vecX];
  const CCTK_REAL vy        = velGF[index_vecY];
  const CCTK_REAL vz        = velGF[index_vecZ];

  const CCTK_REAL my_lapse  = alpGF[index];
  const CCTK_REAL my_shiftx = betaxGF[index];
  const CCTK_REAL my_shifty = betayGF[index];
  const CCTK_REAL my_shiftz = betazGF[index];

  const CCTK_REAL posx      = x[index]-cms_x;
  const CCTK_REAL posy      = y[index]-cms_y;
  const CCTK_REAL posz      = z[index];

  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz,gyy,gyz,gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if(lfac > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

// Covariant 3 velocity

  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz};

  const CCTK_REAL v_2 = v_cov[0]*vx + v_cov[1]*vy + v_cov[2]*vz;

  const CCTK_REAL vz_2 = v_cov[2]*vz;

// Scalar product with shift

  const CCTK_REAL shift_v = v_cov[0]*my_shiftx + v_cov[1]*my_shifty + v_cov[2]*my_shiftz;

  const CCTK_REAL shiftz_vz = v_cov[2]*my_shiftz;

// This is just the total kinetic energy

  VolIntegrand1[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(v_2-shift_v/my_lapse);

// Cylindrical coordinates and some 
// quantities for transformation of vectors. 

  const double posr2      = posx*posx+posy*posy;
  const double posr       = sqrt(posr2);
 
  double xy_over_r2   = posx*posy/posr2;
  double x2_over_r2   = posx*posx/posr2;
  double y2_over_r2   = posy*posy/posr2;

  if ( posr <=  1e-15 ) {
  
	 xy_over_r2   = 0.0; 
	 x2_over_r2   = 0.5; 
	 y2_over_r2   = 0.5; 
  }

// Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2 {0.0};

  v_phi2 = v_cov[0]*vx*y2_over_r2
	  - v_cov[0]*vy*xy_over_r2
	  - v_cov[1]*vx*xy_over_r2
	  + v_cov[1]*vy*x2_over_r2;

// Finally, compute v_phi*shift^phi

  CCTK_REAL v_phi_shift_phi {0.0};

  v_phi_shift_phi = v_cov[0]*my_shiftx*y2_over_r2
	  - v_cov[0]*my_shifty*xy_over_r2
	  - v_cov[1]*my_shiftx*xy_over_r2
	  + v_cov[1]*my_shifty*x2_over_r2;

// This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(v_phi2-v_phi_shift_phi/my_lapse);

// Kinetic energy in z-direction

  VolIntegrand3[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*(vz_2-shiftz_vz/my_lapse);

}

/* Kinetic energy density given by 0.5*rhoh*W^2v^2 */
inline void kinetic(double *VolIntegrand1, double *VolIntegrand2, double *VolIntegrand3, const int index,
		    const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict velGF,
                          const CCTK_REAL *restrict w_lorentzGF, 
			  const CCTK_REAL *restrict rho0, const CCTK_REAL *restrict epsGF, 
			  const CCTK_REAL *restrict pressGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
  			  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {
  
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL my_rho    = rho0[index];
  const CCTK_REAL lfac      = w_lorentzGF[index];
  const CCTK_REAL my_press  = pressGF[index];
  const CCTK_REAL my_eps    = epsGF[index];

  const CCTK_REAL vx        = velGF[index_vecX];
  const CCTK_REAL vy        = velGF[index_vecY];
  const CCTK_REAL vz        = velGF[index_vecZ];

  const CCTK_REAL posx      = x[index]-cms_x;
  const CCTK_REAL posy      = y[index]-cms_y;
  const CCTK_REAL posz      = z[index];

  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz,gyy,gyz,gzz);

  CCTK_REAL w_lorentz_limited = lfac;

  if(lfac > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

// Covariant 3 velocity

  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz};

  const CCTK_REAL v_2 = v_cov[0]*vx + v_cov[1]*vy + v_cov[2]*vz;

  const CCTK_REAL vz_2 = v_cov[2]*vz;

// This is just the total kinetic energy

  VolIntegrand1[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*(w_lorentz_limited*w_lorentz_limited*v_2);

// Cylindrical coordinates and some 
// quantities for transformation of vectors. 

  const double posr2      = posx*posx+posy*posy;
  const double posr       = sqrt(posr2);
 
  double xy_over_r2   = posx*posy/posr2;
  double x2_over_r2   = posx*posx/posr2;
  double y2_over_r2   = posy*posy/posr2;

  if ( posr <=  1e-15 ) {
  
	 xy_over_r2   = 0.0; 
	 x2_over_r2   = 0.5; 
	 y2_over_r2   = 0.5; 
  }

// Finally, compute v_phi*v^phi

  CCTK_REAL v_phi2 {0.0};

  v_phi2 = v_cov[0]*vx*y2_over_r2
	  - v_cov[0]*vy*xy_over_r2
	  - v_cov[1]*vx*xy_over_r2
	  + v_cov[1]*vy*x2_over_r2;

// This is the rotational kinetic energy wrt to the z-axis

  VolIntegrand2[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*w_lorentz_limited*w_lorentz_limited*v_phi2;

// Kinetic energy in z-direction

  VolIntegrand3[index] = 0.5*sqrtgamma*(my_rho*(1.0+my_eps)+my_press)*(w_lorentz_limited*w_lorentz_limited*vz_2);

}

/* Total kinetic energy as defined by total energy in hydro sector - rest mass - thermal energy */
inline void kinetic_tot(double *VolIntegrand1, const int index,
		    const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict velGF,
                          const CCTK_REAL *restrict w_lorentzGF,
			  const CCTK_REAL *restrict rho0, const CCTK_REAL *restrict epsGF, 
			  const CCTK_REAL *restrict pressGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
  			  const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {
  
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL my_rho    = rho0[index];
  const CCTK_REAL my_press  = pressGF[index];
  const CCTK_REAL my_eps    = epsGF[index];
  
  CCTK_REAL w_lorentz = w_lorentzGF[index];

  const CCTK_REAL vx        = velGF[index_vecX];
  const CCTK_REAL vy        = velGF[index_vecY];
  const CCTK_REAL vz        = velGF[index_vecZ];

  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz,gyy,gyz,gzz);

// Covariant 3 velocity

  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz};

  const CCTK_REAL v_2 = v_cov[0]*vx + v_cov[1]*vy + v_cov[2]*vz;

  const CCTK_REAL w_lorentz_2_inv = 1.0 - v_2;

  if(w_lorentz_2_inv > 0.0){
	  w_lorentz = sqrt(1./w_lorentz_2_inv);
  }

  CCTK_REAL w_lorentz_limited = w_lorentz;

  if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

  CCTK_REAL wm1 = w_lorentz_limited -1.0;

  if(wm1< 0.0) wm1 = 0.0;  

  VolIntegrand1[index] = sqrtgamma*( my_press*(w_lorentz_limited*w_lorentz_limited*v_2) + 
		                     (my_rho*(1.0+my_eps))*w_lorentz_limited*wm1 );
}


/* Thermal energy:
 * First integrand is just sqrtgamma*W*rho*eps
 * Second integrand is adiabatic heat transfer defined by the following formula 
 *                  sqrtgamma*W*rho* ( eps - eps(rho,Ye,s_adiabat) )
 * s_adiabat is the entropy given by pure advection from previous to current state, hence no heat transfer
 * Important: Second integrand is not working for now !!!
 * Third integrand is sqrtgamma*W*entropy_density
*/
inline void thermal(double *VolIntegrand1, double *VolIntegrand2, double *VolIntegrand3,
		          const int index,
                          const CCTK_REAL *restrict w_lorentzGF, const CCTK_REAL *restrict rho0GF,
			  const CCTK_REAL *restrict epsGF,
			  const CCTK_REAL *restrict entropyGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const double my_baryon_mass) {

  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL w_lorentz = w_lorentzGF[index];
  const CCTK_REAL rho0      = rho0GF[index];
  const CCTK_REAL eps       = epsGF[index];
  const CCTK_REAL entropy = entropyGF[index]; // this is supposed to be the entropy per baryon
                                              // hence, the entropy density is given by entropy*rho0/baryon_mass
					      // for now we use neutron mass

  CCTK_REAL w_lorentz_limited = w_lorentz;
  if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

  double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz, gyy,gyz,gzz);
  
  VolIntegrand1[index] = sqrtgamma*rho0*eps*w_lorentz_limited;

  VolIntegrand2[index] = sqrtgamma*rho0*w_lorentz_limited*(eps);

  // Total entropy
  VolIntegrand3[index] = sqrtgamma*rho0*entropy*w_lorentz_limited/my_baryon_mass;
}

inline void magnetic_co(double *VolIntegrand1, const int index,
                          const CCTK_REAL *restrict smallb2GF,
			  const CCTK_REAL *restrict w_lorentzGF,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz) {

  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL smallb2 = smallb2GF[index];
  const CCTK_REAL w_lorentz = w_lorentzGF[index];

  CCTK_REAL w_lorentz_limited = w_lorentz;

  if(w_lorentz > CoM_integrand_GAMMA_SPEED_LIMIT) w_lorentz_limited = CoM_integrand_GAMMA_SPEED_LIMIT;

  double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz, gyy,gyz,gzz);
  
  VolIntegrand1[index] = sqrtgamma*w_lorentz_limited*smallb2/2.0;
}

/* Magnetic energy split into azimuthal energy component and rest: */
inline void magnetic_tot(double *VolIntegrand1, double *VolIntegrand2, const int index,
		          const int index_vecX,const int index_vecY, const int index_vecZ,
                          const CCTK_REAL *restrict vel,
			  const CCTK_REAL *restrict Bvec,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {

// Notice that CoM_integrand_GAMMA_SPEED_LIMIT is applied to the integrands above
// involving the Lorentz factor w_lorentz. The be consistent we probably 
// would have to limit the velocity vel, too. 
// This is not done so far.

  const CCTK_REAL gammaDD00 = gxx[index];
  const CCTK_REAL gammaDD01 = gxy[index];
  const CCTK_REAL gammaDD02 = gxz[index];
  const CCTK_REAL gammaDD11 = gyy[index];
  const CCTK_REAL gammaDD12 = gyz[index];
  const CCTK_REAL gammaDD22 = gzz[index];

  const CCTK_REAL B_contra [3] {Bvec[index_vecX],Bvec[index_vecY],Bvec[index_vecZ]};

  const CCTK_REAL vx        = vel[index_vecX];
  const CCTK_REAL vy        = vel[index_vecY];
  const CCTK_REAL vz        = vel[index_vecZ];

  CCTK_REAL posx      = x[index]-cms_x;
  CCTK_REAL posy      = y[index]-cms_y;
  CCTK_REAL posz      = z[index];

// Cylindrical coordinates and some 
// quantities for transformation of vectors. 

  const double posr2      = posx*posx+posy*posy;
  const double posr       = sqrt(posr2);
 
  double xy_over_r2   = posx*posy/posr2;
  double x2_over_r2   = posx*posx/posr2;
  double y2_over_r2   = posy*posy/posr2;

  if ( posr <=  1e-15 ) {
  
	 xy_over_r2   = 0.0; 
	 x2_over_r2   = 0.5; 
	 y2_over_r2   = 0.5; 
  }

  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz, gyy,gyz,gzz);
  const double gamma = sqrtgamma*sqrtgamma;

// Covariant quantities

  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz}; 

  const CCTK_REAL B_cov [3] {gammaDD00*B_contra[0]+gammaDD01*B_contra[1]+gammaDD02*B_contra[2],
                                    gammaDD01*B_contra[0]+gammaDD11*B_contra[1]+gammaDD12*B_contra[2],
                                    gammaDD02*B_contra[0]+gammaDD12*B_contra[1]+gammaDD22*B_contra[2]};

// Contravariant electric field * sqrtgamma

  const CCTK_REAL sqrtgammaE_contra [3] {B_cov[1]*v_cov[2]-B_cov[2]*v_cov[1],
                                                B_cov[2]*v_cov[0]-B_cov[0]*v_cov[2],
                                                B_cov[0]*v_cov[1]-B_cov[1]*v_cov[0]}; 

// Covariant electric field * sqrtgamma

  const CCTK_REAL sqrtgammaE_cov [3] {gammaDD00*sqrtgammaE_contra[0]+gammaDD01*sqrtgammaE_contra[1]+gammaDD02*sqrtgammaE_contra[2],
                                             gammaDD01*sqrtgammaE_contra[0]+gammaDD11*sqrtgammaE_contra[1]+gammaDD12*sqrtgammaE_contra[2],
                                             gammaDD02*sqrtgammaE_contra[0]+gammaDD12*sqrtgammaE_contra[1]+gammaDD22*sqrtgammaE_contra[2]};


// Finally, compute E_phi*E^phi*sqrtgamma and B_phi*B_phi*sqrtgamma

  CCTK_REAL E_2 {0.0};
  CCTK_REAL B_2 {0.0};

  E_2 = sqrtgammaE_cov[0]*sqrtgammaE_contra[0]*y2_over_r2
	  - sqrtgammaE_cov[0]*sqrtgammaE_contra[1]*xy_over_r2
	  - sqrtgammaE_cov[1]*sqrtgammaE_contra[0]*xy_over_r2
	  + sqrtgammaE_cov[1]*sqrtgammaE_contra[1]*x2_over_r2;
  E_2 /= sqrtgamma;

  B_2 = B_cov[0]*B_contra[0]*y2_over_r2
	- B_cov[0]*B_contra[1]*xy_over_r2
	- B_cov[1]*B_contra[0]*xy_over_r2
	+ B_cov[1]*B_contra[1]*x2_over_r2;
  B_2 *= sqrtgamma;
   
  VolIntegrand1[index] = (E_2+B_2)/(8.0*M_PI);

// Finally, compute E_r*E^r*sqrtgamma and B_r*B_r*sqrtgamma

  E_2 = 0.0;
  B_2 = 0.0;

  E_2 = sqrtgammaE_cov[0]*sqrtgammaE_contra[0]*x2_over_r2
	  + sqrtgammaE_cov[0]*sqrtgammaE_contra[1]*xy_over_r2
	  + sqrtgammaE_cov[1]*sqrtgammaE_contra[0]*xy_over_r2
	  + sqrtgammaE_cov[1]*sqrtgammaE_contra[1]*y2_over_r2;
  E_2 /= sqrtgamma;

  B_2 = B_cov[0]*B_contra[0]*x2_over_r2
	+ B_cov[0]*B_contra[1]*xy_over_r2
	+ B_cov[1]*B_contra[0]*xy_over_r2
	+ B_cov[1]*B_contra[1]*y2_over_r2;
  B_2 *= sqrtgamma;
   
  VolIntegrand2[index] = (E_2+sqrtgammaE_cov[2]*sqrtgammaE_contra[2]/sqrtgamma
		  +B_2+B_cov[2]*B_contra[2]*sqrtgamma)/(8.0*M_PI);
}

/* Magnetic energy split into azimuthal energy component and rest: */
/* Only for specific density region: from dens_1 to dens_2 */
inline void magnetic_tot_12(double *VolIntegrand1, double *VolIntegrand2, const int index,
		          const int index_vecX,const int index_vecY, const int index_vecZ,
			  const CCTK_REAL *restrict rho0, const CCTK_REAL dens_1, const CCTK_REAL dens_2,
                          const CCTK_REAL *restrict vel,
			  const CCTK_REAL *restrict Bvec,
                          const CCTK_REAL *restrict gxx,const CCTK_REAL *restrict gxy,const CCTK_REAL *restrict gxz,
                          const CCTK_REAL *restrict gyy,const CCTK_REAL *restrict gyz,const CCTK_REAL *restrict gzz,
			  const CCTK_REAL *restrict x,const CCTK_REAL *restrict y,const CCTK_REAL *restrict z,
			  const double cms_x, const double cms_y) {

// Notice that CoM_integrand_GAMMA_SPEED_LIMIT is applied to the integrands above
// involving the Lorentz factor w_lorentz. The be consistent we probably 
// would have to limit the velocity vel, too. 
// This is not done so far.

  const CCTK_REAL my_rho = rho0[index];

  if ( my_rho > dens_1 && my_rho <= dens_2 ) {

	  const CCTK_REAL gammaDD00 = gxx[index];
	  const CCTK_REAL gammaDD01 = gxy[index];
	  const CCTK_REAL gammaDD02 = gxz[index];
	  const CCTK_REAL gammaDD11 = gyy[index];
	  const CCTK_REAL gammaDD12 = gyz[index];
	  const CCTK_REAL gammaDD22 = gzz[index];

	  const CCTK_REAL B_contra [3] {Bvec[index_vecX],Bvec[index_vecY],Bvec[index_vecZ]};

	  const CCTK_REAL vx        = vel[index_vecX];
	  const CCTK_REAL vy        = vel[index_vecY];
	  const CCTK_REAL vz        = vel[index_vecZ];

	  CCTK_REAL posx      = x[index]-cms_x;
	  CCTK_REAL posy      = y[index]-cms_y;
	  CCTK_REAL posz      = z[index];

	// Cylindrical coordinates and some 
	// quantities for transformation of vectors. 

	  const double posr2      = posx*posx+posy*posy;
	  const double posr       = sqrt(posr2);
 
	  double xy_over_r2   = posx*posy/posr2;
	  double x2_over_r2   = posx*posx/posr2;
	  double y2_over_r2   = posy*posy/posr2;

	  if ( posr <=  1e-15 ) {
  
		 xy_over_r2   = 0.0; 
		 x2_over_r2   = 0.5; 
		 y2_over_r2   = 0.5; 
	  }

	  const double sqrtgamma = compute_sqrtgamma(index,gxx,gxy,gxz, gyy,gyz,gzz);
	  const double gamma = sqrtgamma*sqrtgamma;

	// Covariant quantities

	  const CCTK_REAL v_cov [3] {gammaDD00*vx+gammaDD01*vy+gammaDD02*vz,
                                    gammaDD01*vx+gammaDD11*vy+gammaDD12*vz,
                                    gammaDD02*vx+gammaDD12*vy+gammaDD22*vz}; 

	  const CCTK_REAL B_cov [3] {gammaDD00*B_contra[0]+gammaDD01*B_contra[1]+gammaDD02*B_contra[2],
                                    gammaDD01*B_contra[0]+gammaDD11*B_contra[1]+gammaDD12*B_contra[2],
                                    gammaDD02*B_contra[0]+gammaDD12*B_contra[1]+gammaDD22*B_contra[2]};

	// Contravariant electric field * sqrtgamma

	  const CCTK_REAL sqrtgammaE_contra [3] {B_cov[1]*v_cov[2]-B_cov[2]*v_cov[1],
                                                B_cov[2]*v_cov[0]-B_cov[0]*v_cov[2],
                                                B_cov[0]*v_cov[1]-B_cov[1]*v_cov[0]}; 

	// Covariant electric field * sqrtgamma

	  const CCTK_REAL sqrtgammaE_cov [3] {gammaDD00*sqrtgammaE_contra[0]+gammaDD01*sqrtgammaE_contra[1]+gammaDD02*sqrtgammaE_contra[2],
                                             gammaDD01*sqrtgammaE_contra[0]+gammaDD11*sqrtgammaE_contra[1]+gammaDD12*sqrtgammaE_contra[2],
                                             gammaDD02*sqrtgammaE_contra[0]+gammaDD12*sqrtgammaE_contra[1]+gammaDD22*sqrtgammaE_contra[2]};


	// Finally, compute E_phi*E^phi*sqrtgamma and B_phi*B_phi*sqrtgamma

	  CCTK_REAL E_2 {0.0};
	  CCTK_REAL B_2 {0.0};

	  E_2 = sqrtgammaE_cov[0]*sqrtgammaE_contra[0]*y2_over_r2
		  - sqrtgammaE_cov[0]*sqrtgammaE_contra[1]*xy_over_r2
		  - sqrtgammaE_cov[1]*sqrtgammaE_contra[0]*xy_over_r2
		  + sqrtgammaE_cov[1]*sqrtgammaE_contra[1]*x2_over_r2;
	  E_2 /= sqrtgamma;

	  B_2 = B_cov[0]*B_contra[0]*y2_over_r2
		- B_cov[0]*B_contra[1]*xy_over_r2
		- B_cov[1]*B_contra[0]*xy_over_r2
		+ B_cov[1]*B_contra[1]*x2_over_r2;
	  B_2 *= sqrtgamma;
   
	  VolIntegrand1[index] = (E_2+B_2)/(8.0*M_PI);

	// Finally, compute E_r*E^r*sqrtgamma and B_r*B_r*sqrtgamma

	  E_2 = 0.0;
	  B_2 = 0.0;

	  E_2 = sqrtgammaE_cov[0]*sqrtgammaE_contra[0]*x2_over_r2
		  + sqrtgammaE_cov[0]*sqrtgammaE_contra[1]*xy_over_r2
		  + sqrtgammaE_cov[1]*sqrtgammaE_contra[0]*xy_over_r2
		  + sqrtgammaE_cov[1]*sqrtgammaE_contra[1]*y2_over_r2;
	  E_2 /= sqrtgamma;

	  B_2 = B_cov[0]*B_contra[0]*x2_over_r2
		+ B_cov[0]*B_contra[1]*xy_over_r2
		+ B_cov[1]*B_contra[0]*xy_over_r2
		+ B_cov[1]*B_contra[1]*y2_over_r2;
	  B_2 *= sqrtgamma;
   
	  VolIntegrand2[index] = (E_2+sqrtgammaE_cov[2]*sqrtgammaE_contra[2]/sqrtgamma
			  +B_2+B_cov[2]*B_contra[2]*sqrtgamma)/(8.0*M_PI);

  } else {

	  VolIntegrand1[index] = 0.0;
	  VolIntegrand2[index] = 0.0;
  }

}
/* Neutrino Luminosity: */
//inline void Neutrino_lum_integrand(double *VolIntegrand1, double *VolIntegrand2,double *VolIntegrand3,const int index,
//                          const CCTK_REAL *restrict lum_nue,const CCTK_REAL *restrict lum_nua,const CCTK_REAL *restrict lum_nux) {
//  VolIntegrand1[index] = lum_nue[index];
//  VolIntegrand2[index] = lum_nua[index];
//  VolIntegrand3[index] = lum_nux[index];
//}
#endif
