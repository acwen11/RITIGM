#include <math.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "GlobalDerivative.h"
#include "loopcontrol.h"
#include "util.h"

void
CTGBase_Convert_CTG_to_ADM_detg(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT istart[3], iend[3];
  CCTK_REAL n = *detg_exponent;

  if (verbose) CCTK_INFO("Converting CTGamma to ADM variables (detg^n version).");

  Util_GetGridRanges(cctkGH, istart, iend);

#pragma omp parallel
  LC_LOOP3 (CTGBase_Convert_CTG_to_ADM_detg,
            i, j, k,
            0, 0, 0,
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      CCTK_REAL e4phi = pow(phi[ijk],1.0/(3.0*n));

      gxx[ijk] = e4phi * gamma11[ijk];
      gxy[ijk] = e4phi * gamma12[ijk];
      gxz[ijk] = e4phi * gamma13[ijk];
      gyy[ijk] = e4phi * gamma22[ijk];
      gyz[ijk] = e4phi * gamma23[ijk];
      gzz[ijk] = e4phi * gamma33[ijk];

			if (alp[ijk] < 0 || gxx[ijk] < 0 || gyy[ijk] < 0 || gzz[ijk] < 0
					|| isnan(alp[ijk] * gxx[ijk] * gyy[ijk] * gzz[ijk])) {
				CCTK_VINFO("PERFORMING SPACETIME RESET HACK at ijk = %d %d %d, xyz = %e %e %e!", 
					i, j, k, x[ijk], y[ijk], z[ijk]);
				CCTK_VINFO("g_ij = %e %e %e %e %e %e", gxx[ijk], gxy[ijk], gxz[ijk], gyy[ijk], gyz[ijk], gzz[ijk]);
				CCTK_VINFO("Z4c gamma = %e %e %e %e %e %e", gamma11[ijk], gamma12[ijk], gamma13[ijk], gamma22[ijk], gamma23[ijk], gamma33[ijk]);
				CCTK_VINFO("Z4c conformal fac = %e", phi[ijk]);

				if( GetRefinementLevel(cctkGH) > 6 ) {
					CCTK_ERROR("Spacetime error detected in Z4c!");
				}
				else {
					CCTK_VWARN(CCTK_WARN_ALERT,"Found spacetime error in Z4c, but not at finest level. Proceeding with caution...");
          CCTK_VWARN(CCTK_WARN_ALERT, "PERFORMING SPACETIME RESET HACK at ijk = %d %d %d, xyz = %e %e %e!",
            i, j, k, x[ijk], y[ijk], z[ijk]);
				}

				// Bad hack to reset spacetime gfs
				alp[ijk] = alp_p[ijk];
				betax[ijk] = betax_p[ijk];
				betay[ijk] = betay_p[ijk];
				betaz[ijk] = betaz_p[ijk];
				gxx[ijk] = gxx_p[ijk];
				gxy[ijk] = gxy_p[ijk];
				gxz[ijk] = gxz_p[ijk];
				gyy[ijk] = gyy_p[ijk];
				gyz[ijk] = gyz_p[ijk];
				gzz[ijk] = gzz_p[ijk];

				phi[ijk] = phi_p[ijk];
				gamma11[ijk] = gamma11_p[ijk];
				gamma12[ijk] = gamma12_p[ijk];
				gamma13[ijk] = gamma13_p[ijk];
				gamma22[ijk] = gamma22_p[ijk];
				gamma23[ijk] = gamma23_p[ijk];
				gamma33[ijk] = gamma33_p[ijk];
				
				K[ijk] = K_p[ijk];
				Khat[ijk] = Khat_p[ijk];
				Theta[ijk] = Theta_p[ijk];

				A11[ijk] = A11_p[ijk];
				A12[ijk] = A12_p[ijk];
				A13[ijk] = A13_p[ijk];
				A22[ijk] = A22_p[ijk];
				A23[ijk] = A23_p[ijk];
				A33[ijk] = A33_p[ijk];

				Gamma1[ijk] = Gamma1_p[ijk];
				Gamma2[ijk] = Gamma2_p[ijk];
				Gamma3[ijk] = Gamma3_p[ijk];

				CCTK_VINFO("After fix at ijk = %d %d %d, xyz = %e %e %e!", i, j, k, x[ijk], y[ijk], z[ijk]);
				CCTK_VINFO("g_ij = %e %e %e %e %e %e", gxx[ijk], gxy[ijk], gxz[ijk], gyy[ijk], gyz[ijk], gzz[ijk]);
				CCTK_VINFO("Z4c gamma = %e %e %e %e %e %e", gamma11[ijk], gamma12[ijk], gamma13[ijk], gamma22[ijk], gamma23[ijk], gamma33[ijk]);
				CCTK_VINFO("Z4c conformal fac = %e", phi[ijk]);
    	}  

      CCTK_REAL trKfac = K[ijk] / 3.0;

      kxx[ijk] = e4phi * (A11[ijk] + trKfac*gamma11[ijk]);
      kxy[ijk] = e4phi * (A12[ijk] + trKfac*gamma12[ijk]);
      kxz[ijk] = e4phi * (A13[ijk] + trKfac*gamma13[ijk]);
      kyy[ijk] = e4phi * (A22[ijk] + trKfac*gamma22[ijk]);
      kyz[ijk] = e4phi * (A23[ijk] + trKfac*gamma23[ijk]);
      kzz[ijk] = e4phi * (A33[ijk] + trKfac*gamma33[ijk]);
	}
  LC_ENDLOOP3 (CTGBase_Convert_CTG_to_ADM_detg);

  return;
}


void
CTGBase_Convert_ADM_to_CTG_detg(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ni, nj, nk; 
  CCTK_INT istart[3], iend[3];
  CCTK_INT * restrict imin[3], * restrict imax[3];
  CCTK_REAL * restrict q[3];
  CCTK_REAL ihx, ihy, ihz;
  CCTK_REAL n = *detg_exponent;

  DECLARE_EVOL_TYPE;

  CCTK_INFO("Converting ADM to CTGamma variables (detg^n version).");

  ni = cctk_lsh[0];
  nj = cctk_lsh[1];
  nk = cctk_lsh[2];

  /*
   * Grid spacings required by the finite difference operators.
   */
  ihx = 1.0 / CCTK_DELTA_SPACE(0);
  ihy = 1.0 / CCTK_DELTA_SPACE(1);
  ihz = 1.0 / CCTK_DELTA_SPACE(2);

  Util_SetStencil(cctkGH, (CCTK_POINTER_TO_CONST *) imin, (CCTK_POINTER_TO_CONST *) imax, (CCTK_POINTER_TO_CONST *) q, 1, 0);
  Util_GetGridRanges(cctkGH, istart, iend);

#pragma omp parallel
  LC_LOOP3 (CTGBase_Convert_ADM_to_CTG_detg0,
            i, j, k,
            0, 0, 0,
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      CCTK_REAL ig11, ig12, ig13, ig22, ig23, ig33;
      CCTK_REAL local_detg;
      CCTK_REAL psi4;
      CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

      CTGBase_invert_metric(gxx[ijk], gxy[ijk], gxz[ijk], gyy[ijk],
                            gyz[ijk], gzz[ijk],
                            local_detg,
                            ig11, ig12, ig13, ig22, ig23, ig33);

      detg[ijk] = local_detg;

      psi4 = pow(local_detg, 0.333333333333333333);

      igamma11[ijk] = ig11 * psi4;
      igamma12[ijk] = ig12 * psi4;
      igamma13[ijk] = ig13 * psi4;
      igamma22[ijk] = ig22 * psi4;
      igamma23[ijk] = ig23 * psi4;
      igamma33[ijk] = ig33 * psi4;
    }
  LC_ENDLOOP3 (CTGBase_Convert_ADM_to_CTG_detg0);

#pragma omp parallel
  LC_LOOP3 (CTGBase_Convert_ADM_to_CTG_detg1,
            i, j, k,
            0, 0, 0,
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2],
            cctk_lsh[0], cctk_lsh[1], cctk_lsh[2])
    {
      CCTK_REAL dadx=1, dbdx=0, dcdx=0, dady=0, dbdy=1, dcdy=0, dadz=0, dbdz=0,
        dcdz=1;

      CCTK_REAL digamma111=-424242.0, digamma121=-424242.0,
	digamma122=-424242.0;
      CCTK_REAL digamma131=-424242.0, digamma133=-424242.0,
	digamma222=-424242.0;
      CCTK_REAL digamma232=-424242.0, digamma233=-424242.0,
	digamma333=-424242.0;

      CCTK_INT ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
      CCTK_REAL pm4 = 1.0 / pow(detg[ijk], 0.333333333333333333);
      CCTK_REAL ig11 = igamma11[ijk] * pm4;
      CCTK_REAL ig12 = igamma12[ijk] * pm4;
      CCTK_REAL ig13 = igamma13[ijk] * pm4;
      CCTK_REAL ig22 = igamma22[ijk] * pm4;
      CCTK_REAL ig23 = igamma23[ijk] * pm4;
      CCTK_REAL ig33 = igamma33[ijk] * pm4;

      CCTK_REAL g11 = gxx[ijk];
      CCTK_REAL g12 = gxy[ijk];
      CCTK_REAL g13 = gxz[ijk];
      CCTK_REAL g22 = gyy[ijk];
      CCTK_REAL g23 = gyz[ijk];
      CCTK_REAL g33 = gzz[ijk];
          
      CCTK_REAL k11 = kxx[ijk];
      CCTK_REAL k12 = kxy[ijk];
      CCTK_REAL k13 = kxz[ijk];
      CCTK_REAL k22 = kyy[ijk];
      CCTK_REAL k23 = kyz[ijk];
      CCTK_REAL k33 = kzz[ijk];

      CCTK_REAL trk = ig11*k11 + ig22*k22 + ig33*k33 + 
        2.0*(ig12*k12 + ig13*k13 + ig23*k23);

      if (*general_coordinates)
        {
          dadx = J11[ijk];
          dady = J12[ijk];
          dadz = J13[ijk];
          dbdx = J21[ijk];
          dbdy = J22[ijk];
          dbdz = J23[ijk];
          dcdx = J31[ijk];
          dcdy = J32[ijk];
          dcdz = J33[ijk];
        }

      digamma111 = g_diff_dx(cctkGH, igamma11,
                             *general_coordinates,
                             dadx, dbdx, dcdx,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma121 = g_diff_dx(cctkGH, igamma12,
                             *general_coordinates,
                             dadx, dbdx, dcdx,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma131 = g_diff_dx(cctkGH, igamma13,
                             *general_coordinates,
                             dadx, dbdx, dcdx,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);


      digamma122 = g_diff_dy(cctkGH, igamma12,
                             *general_coordinates,
                             dady, dbdy, dcdy,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma222 = g_diff_dy(cctkGH, igamma22,
                             *general_coordinates,
                             dady, dbdy, dcdy,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma232 = g_diff_dy(cctkGH, igamma23,
                             *general_coordinates,
                             dady, dbdy, dcdy,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);

      digamma133 = g_diff_dz(cctkGH, igamma13,
                             *general_coordinates,
                             dadz, dbdz, dcdz,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma233 = g_diff_dz(cctkGH, igamma23,
                             *general_coordinates,
                             dadz, dbdz, dcdz,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
      digamma333 = g_diff_dz(cctkGH, igamma33,
                             *general_coordinates,
                             dadz, dbdz, dcdz,
                             i, j, k, ni, nj, nk,
                             imin[0], imax[0], imin[1], imax[1],
                             imin[2], imax[2], q[0], q[1], q[2],
                             ihx, ihy, ihz);
	  
      phi[ijk] = pow(detg[ijk], n);
          
      gamma11[ijk] = pm4 * gxx[ijk];
      gamma12[ijk] = pm4 * gxy[ijk];
      gamma13[ijk] = pm4 * gxz[ijk];
      gamma22[ijk] = pm4 * gyy[ijk];
      gamma23[ijk] = pm4 * gyz[ijk];
      gamma33[ijk] = pm4 * gzz[ijk];

      K[ijk] = trk;

      A11[ijk] = pm4 * (k11 - g11*trk/3.0);
      A12[ijk] = pm4 * (k12 - g12*trk/3.0);
      A13[ijk] = pm4 * (k13 - g13*trk/3.0);
      A22[ijk] = pm4 * (k22 - g22*trk/3.0);
      A23[ijk] = pm4 * (k23 - g23*trk/3.0);
      A33[ijk] = pm4 * (k33 - g33*trk/3.0);

      Gamma1[ijk] = -(digamma111 + digamma122 + digamma133);
      Gamma2[ijk] = -(digamma121 + digamma222 + digamma233);
      Gamma3[ijk] = -(digamma131 + digamma232 + digamma333);
      
      if (EVOL_TYPE == Z4c_TYPE) {
        Theta[ijk] = 0;
        Khat[ijk] = K[ijk];
      }
    }
  LC_ENDLOOP3 (CTGBase_Convert_ADM_to_CTG_detg1);

  CTGBase_free_stencil(imin, imax, q);

  return;
}
