#include "cctk.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

const double ONE_OVER_SQRT_4PI = 1/sqrt(4*M_PI);
template<typename T>
inline T SQR(T t) { return t*t; }

void compute_bi_b2_Poyn_fluxET(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(smallbPoynET_compute_every<=0 || cctk_iteration%smallbPoynET_compute_every!=0) return;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	double gxxL = gxx[index];
	double gxyL = gxy[index];
	double gxzL = gxz[index];
	double gyyL = gyy[index];
	double gyzL = gyz[index];
	double gzzL = gzz[index];

	double det = -gxzL*gxzL*gyyL + 2*gxyL*gxzL*gyzL - gxxL*gyzL*gyzL - gxyL*gxyL*gzzL + gxxL*gyyL*gzzL;
	double invdet = 1.0 / det;
	double gupxxL=(-gyzL*gyzL + gyyL*gzzL)*invdet;
	double gupxyL=( gxzL*gyzL - gxyL*gzzL)*invdet;
	double gupyyL=(-gxzL*gxzL + gxxL*gzzL)*invdet;
	double gupxzL=(-gxzL*gyyL + gxyL*gyzL)*invdet;
	double gupyzL=( gxyL*gxzL - gxxL*gyzL)*invdet;
	double gupzzL=(-gxyL*gxyL + gxxL*gyyL)*invdet;

	double lapse = alp[index];
	double shiftx = betax[index];
	double shifty = betay[index];
	double shiftz = betaz[index];

	double ONE_OVER_LAPSE = 1.0/alp[index];
	double ONE_OVER_LAPSE_SQRT_4PI = ONE_OVER_LAPSE*ONE_OVER_SQRT_4PI;

	double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

	// IllinoisGRMHD defines v^i = u^i/u^0.
        
	// Meanwhile, the ET/HydroBase formalism, called the Valencia 
	// formalism, splits the 4 velocity into a purely spatial part
	// and a part that is normal to the spatial hypersurface:
	// u^a = G (n^a + U^a), (Eq. 14 of arXiv:1304.5544; G=W, U^a=v^a)
	// where n^a is the unit normal vector to the spatial hypersurface,
	// n_a = {-\alpha,0,0,0}, and U^a is the purely spatial part, which
	// is defined in HydroBase as the vel[] vector gridfunction.
	// Then u^a n_a = - \alpha u^0 = G n^a n_a = -G, and
	// of course \alpha u^0 = 1/sqrt(1+γ^ij u_i u_j) = \Gamma,
	// the standard Lorentz factor.

	// Note that n^i = - \beta^i / \alpha, so 
	// u^a = \Gamma (n^a + U^a) 
	// -> u^i = \Gamma ( U^i - \beta^i / \alpha )
	// which implies
	// v^i = u^i/u^0
	//     = \Gamma/u^0 ( U^i - \beta^i / \alpha ) <- \Gamma = \alpha u^0
	//     = \alpha ( U^i - \beta^i / \alpha )
	//     = \alpha U^i - \beta^i

	double vxL = lapse*ETvx - shiftx;
	double vyL = lapse*ETvy - shifty;
	double vzL = lapse*ETvz - shiftz;

	// Derivation of first equation:
	// \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2 
	//   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
	//   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
	//   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
	//   = 1 - 1/(u^0 \alpha)^2 <= 1
	double one_minus_one_over_alpha_u0_squared = (gxxL* SQR(vxL + shiftx) +
						      2.0*gxyL*(vxL + shiftx)*(vyL + shifty) +
						      2.0*gxzL*(vxL + shiftx)*(vzL + shiftz) +
						      gyyL* SQR(vyL + shifty) +
						      2.0*gyzL*(vyL + shifty)*(vzL + shiftz) +
						      gzzL* SQR(vzL + shiftz) )*SQR(ONE_OVER_LAPSE);
	/*** Check for superluminal velocity ***/
	/* Don't bother with speed limits in this analysis thorn.
	//FIXME: Instead of >1.0, should be one_minus_one_over_alpha_u0_squared > ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED, for consistency with conserv_to_prims routines
	if(one_minus_one_over_alpha_u0_squared > 1.0) {
	double ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED = 1.0-1.0/SQR(GAMMA_SPEED_LIMIT);
	double correction_fac = sqrt(ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED/one_minus_one_over_alpha_u0_squared);
	vxL = (vxL + shiftx)*correction_fac-shiftx;
	vyL = (vyL + shifty)*correction_fac-shifty;
	vzL = (vzL + shiftz)*correction_fac-shiftz;
	one_minus_one_over_alpha_u0_squared=ONE_MINUS_ONE_OVER_GAMMA_SPEED_LIMIT_SQUARED;
	stats.failure_checker+=1000;
	}
	*/
	// A = 1.0-one_minus_one_over_alpha_u0_squared = 1-(1-1/(al u0)^2) = 1/(al u0)^2
	// 1/sqrt(A) = al u0
	//double alpha_u0_minus_one = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared)-1.0;
	//u0_out          = (alpha_u0_minus_one + 1.0)*ONE_OVER_LAPSE;
	double alpha_u0 = 1.0/sqrt(1.0-one_minus_one_over_alpha_u0_squared);
	double u0L = alpha_u0*ONE_OVER_LAPSE;
 
	double Bx_center = Bvec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
	double By_center = Bvec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
	double Bz_center = Bvec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

	// NOW COMPUTE b^{\mu} and b^2 = b^{\mu} b^{\nu} g_{\mu \nu}
	double ONE_OVER_U0 = 1.0/u0L;
	double shiftx_plus_vx = (shiftx+vxL);
	double shifty_plus_vy = (shifty+vyL);
	double shiftz_plus_vz = (shiftz+vzL);

	// Eq. 56 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	//  u_i = gamma_{ij} u^0 (v^j + beta^j), gamma_{ij} is the physical metric
	double u_x_over_u0 =  gxxL*shiftx_plus_vx + gxyL*shifty_plus_vy + gxzL*shiftz_plus_vz;
	double u_y_over_u0 =  gxyL*shiftx_plus_vx + gyyL*shifty_plus_vy + gyzL*shiftz_plus_vz;
	double u_z_over_u0 =  gxzL*shiftx_plus_vx + gyzL*shifty_plus_vy + gzzL*shiftz_plus_vz;

	// Eqs. 23 and 31 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	//   Compute alpha sqrt(4 pi) b^t = u_i B^i
	double alpha_sqrt_4pi_bt = ( u_x_over_u0*Bx_center + u_y_over_u0*By_center + u_z_over_u0*Bz_center ) * u0L;
	// Eq. 24 in http://arxiv.org/pdf/astro-ph/0503420.pdf:
	// b^i = B^i_u / sqrt(4 pi)
	// b^i = ( B^i/alpha + B^0_u u^i ) / ( u^0 sqrt(4 pi) )
	// b^i = ( B^i/alpha +  sqrt(4 pi) b^t u^i ) / ( u^0 sqrt(4 pi) )
	// b^i = ( B^i +  alpha sqrt(4 pi) b^t u^i ) / ( alpha u^0 sqrt(4 pi) )
	// b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t u^i/u^0 ) / ( alpha sqrt(4 pi) )
	// b^i = ( B^i/u^0 +  alpha sqrt(4 pi) b^t v^i ) / ( alpha sqrt(4 pi) )
	double smallbxL = (Bx_center*ONE_OVER_U0 + vxL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	double smallbyL = (By_center*ONE_OVER_U0 + vyL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	double smallbzL = (Bz_center*ONE_OVER_U0 + vzL*alpha_sqrt_4pi_bt)*ONE_OVER_LAPSE_SQRT_4PI;
	// Eq. 23 in http://arxiv.org/pdf/astro-ph/0503420.pdf, with alpha sqrt (4 pi) b^2 = u_i B^i already computed above
	double smallbtL = alpha_sqrt_4pi_bt * ONE_OVER_LAPSE_SQRT_4PI;

	/* POYNTING FLUX: */
	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) 
	double LAPSE_SQUARED = lapse*lapse;
	double ONE_OVER_LAPSE_SQUARED = ONE_OVER_LAPSE*ONE_OVER_LAPSE;
	double g4dn[4][4],g4up[4][4];
    
	// g^{\mu \nu} = upper four-metric.
	g4up[0][0]              = -ONE_OVER_LAPSE_SQUARED;
	g4up[0][1] = g4up[1][0] = ONE_OVER_LAPSE_SQUARED*shiftx;
	g4up[0][2] = g4up[2][0] = ONE_OVER_LAPSE_SQUARED*shifty;
	g4up[0][3] = g4up[3][0] = ONE_OVER_LAPSE_SQUARED*shiftz;
	g4up[1][1]              = gupxxL - ONE_OVER_LAPSE_SQUARED*shiftx*shiftx;
	g4up[1][2] = g4up[2][1] = gupxyL - ONE_OVER_LAPSE_SQUARED*shiftx*shifty;
	g4up[1][3] = g4up[3][1] = gupxzL - ONE_OVER_LAPSE_SQUARED*shiftx*shiftz;
	g4up[2][2]              = gupyyL - ONE_OVER_LAPSE_SQUARED*shifty*shifty;
	g4up[2][3] = g4up[3][2] = gupyzL - ONE_OVER_LAPSE_SQUARED*shifty*shiftz;
	g4up[3][3]              = gupzzL - ONE_OVER_LAPSE_SQUARED*shiftz*shiftz;

	/* Compute beta_i */
	double shift_x = ( gxxL*(shiftx) + gxyL*(shifty) +
			   gxzL*(shiftz) );
	double shift_y = ( gxyL*(shiftx) + gyyL*(shifty) +
			   gyzL*(shiftz) );
	double shift_z = ( gxzL*(shiftx) + gyzL*(shifty) +
			   gzzL*(shiftz) );
	
	// g_{00}               = - alpha^2 + gamma_{ij} beta^i beta^j = - alpha^2 beta_i beta^i
	g4dn[0][0]              = -LAPSE_SQUARED + (shiftx*shift_x + shifty*shift_y + shiftz*shift_z);
	// g_{0i} =  gamma_{ij} beta^j = beta_i
	g4dn[0][1] = g4dn[1][0] = shift_x;
	g4dn[0][2] = g4dn[2][0] = shift_y;
	g4dn[0][3] = g4dn[3][0] = shift_z;
	// g_{ij} =  gamma_{ij} <- 3 metric
	g4dn[1][1] =              gxxL;
	g4dn[1][2] = g4dn[2][1] = gxyL;
	g4dn[1][3] = g4dn[3][1] = gxzL;
	g4dn[2][2] =              gyyL;
	g4dn[2][3] = g4dn[3][2] = gyzL;
	g4dn[3][3] =              gzzL;

	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) 

	// First compute u_0.
	double uup[4] = {u0L,u0L*vxL,u0L*vyL,u0L*vzL};
	double u_0=0.0; for(int ii=0;ii<4;ii++) u_0 += g4dn[0][ii]*uup[ii];

	// Next compute b_0.
	double sbup[4] = {smallbtL,smallbxL,smallbyL,smallbzL};
	double smallb_0L=0.0; for(int ii=0;ii<4;ii++) smallb_0L += g4dn[0][ii]*sbup[ii];

	// Next compute g^i_0:
	double g4upx_0=0.0; for(int ii=0;ii<4;ii++) g4upx_0 += g4up[1][ii]*g4dn[0][ii];
	double g4upy_0=0.0; for(int ii=0;ii<4;ii++) g4upy_0 += g4up[2][ii]*g4dn[0][ii];
	double g4upz_0=0.0; for(int ii=0;ii<4;ii++) g4upz_0 += g4up[3][ii]*g4dn[0][ii];

	// Next compute b^2:
	double smallb[4] = {smallbtL,smallbxL,smallbyL,smallbzL};
	double smallb2L=0.0; for(int ii=0;ii<4;ii++) for(int jj=0;jj<4;jj++) smallb2L += g4dn[ii][jj]*smallb[ii]*smallb[jj];

	// S^i = -alp*Tem^i_0 = -alp*(b^2 u^i u_0 + b^2/2 g^i_0 - b^i b_0) 
	double PoynxL = -lapse*(smallb2L*(u0L*vxL)*u_0 + 0.5*smallb2L*g4upx_0 - smallbxL*smallb_0L);
	double PoynyL = -lapse*(smallb2L*(u0L*vyL)*u_0 + 0.5*smallb2L*g4upy_0 - smallbyL*smallb_0L);
	double PoynzL = -lapse*(smallb2L*(u0L*vzL)*u_0 + 0.5*smallb2L*g4upz_0 - smallbzL*smallb_0L);

	// Magnetic Energy Calculation
	// double B_x = gxxL*(Bx_center) + gxyL*(By_center) + gxzL*(Bz_center);
	// double B_y = gxyL*(Bx_center) + gyyL*(By_center) + gyzL*(Bz_center);
	// double B_z = gxzL*(Bx_center) + gyzL*(By_center) + gzzL*(Bz_center);
	// double B2big = Bx_center * B_x + By_center * B_y + Bz_center * B_z;

	// TODO: Check this. Multiply by 4pi to convert back to Gaussian units.
	// double EMenergyL = 4 * M_PI * (B2big - 0.5 * smallb2L) * std::sqrt(det); 
	
	// From Baumgarte + Shapiro eq. 5.152
	double EMenergyL = (smallb2L * (std::pow(alpha_u0, 2) - 0.5) - std::pow(lapse * smallbtL, 2)) * std::sqrt(det);

	minus_one_minus_u_0[index] = -1.0 - u_0;

	smallbt[index] = smallbtL;
	smallbx[index] = smallbxL;
	smallby[index] = smallbyL;
	smallbz[index] = smallbzL;
	smallb2[index] = smallb2L;

	Poynx[index] = PoynxL;
	Poyny[index] = PoynyL;
	Poynz[index] = PoynzL;

	EMenergy[index] = EMenergyL;
	Bmag[index] = std::sqrt(B2big);

	double grid_spacing_L = 1.0;
	for(int ii=0; ii<cctk_dim; ii++){
		grid_spacing_L *= cctk_delta_space[ii];
	}
	*grid_spacing_product_mhd = grid_spacing_L;
      } // Loop over grid
}

void compute_total_energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(smallbPoynET_compute_every<=0 || cctk_iteration%smallbPoynET_compute_every!=0) return;

  const int reduction_handle = CCTK_ReductionHandle("sum");
  const int var_index = CCTK_VarIndex("smallbPoynET::EMenergy");

	double tot_energy_reduce;
  const int ierr = CCTK_Reduce(cctkGH, -1, reduction_handle,
       1, CCTK_VARIABLE_REAL, &tot_energy_reduce, 1, var_index);

  if (ierr != 0){ 
     CCTK_WARN(0, "Failed to integrate EM energy");
	}

  *Total_Energy = tot_energy_reduce * (*grid_spacing_product_mhd);
}
