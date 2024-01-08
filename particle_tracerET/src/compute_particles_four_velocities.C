#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

// This function computes u^{mu} and u_{mu}.
//
// We know that
//
// u^{0} = W / \alpha ,
//
// where the Lorentz factor is given by
//
// W = 1/sqrt(1-v^2) ,
//
// and v^2 = gamma_{ij}v^{i}v^{j}, where
// gamma_{ij} is the spatial metric and
// v^{i} is the valencia 3-velocity.
//
// Finally, we know that
//
// u^{i} = u^{0} ( alpha v^{i} - beta^{i} ) .
//
// To compute u_{mu} we need the 4-metric:
//
// g_{00} = -alpha^2 + beta^2 ,
// g_{0i} = g_{i0} = beta_{i} ,
// g_{ij} = gamma_{ij} .
void compute_particles_four_velocities(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(update_RK4_freq==0) return;

  if(cctk_iteration==start_tracing_particles_iteration || cctk_iteration%(2*update_RK4_freq)==0) {

    if(verbose>1) CCTK_VINFO("Computing four velocity at Ref. Lev. %d...", GetRefinementLevel(cctkGH));

#pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) {
      for(int j=0;j<cctk_lsh[1];j++) {
        for(int i=0;i<cctk_lsh[0];i++) {

          // Set gridfunction index
          int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

          // Read from gridfunctions
          const double alpha    = alp[index];
          const double betaU[3] = { betax[index], betay[index], betaz[index] };
          const double gammaDD[3][3] = { {gxx[index], gxy[index], gxz[index]},
                                         {gxy[index], gyy[index], gyz[index]},
                                         {gxz[index], gyz[index], gzz[index]} };
          const double vU[3] = { vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)],
                                 vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)],
                                 vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)] };

          // Compute v^2 = gamma_{ij} v^{i} v^{j}
          double v2 = 0.0;
          for(int ii=0;ii<3;ii++)
            for(int jj=0;jj<3;jj++)
              v2 += gammaDD[ii][jj] * vU[ii] * vU[jj];

          // Compute W = 1/sqrt(1-v2)
          const double W = 1.0/sqrt(1.0-v2);

          // Compute u^{mu}
          double u4U[4] = {W/alpha, 0.0, 0.0, 0.0};
          for(int ii=0;ii<3;ii++)
            u4U[ii+1] = u4U[0] * ( alpha * vU[ii] - betaU[ii] );

          if( output_four_velocity_u4U ) {
            // Write to gridfunctions
            u4U0GF[index] = u4U[0];
            u4U1GF[index] = u4U[1];
            u4U2GF[index] = u4U[2];
            u4U3GF[index] = u4U[3];
          }

          if( output_four_velocity_u4D ) {
            // Compute beta_{i} = gamma_{ij} beta^{j}
            double betaD[3] = {0.0, 0.0, 0.0};
            for(int ii=0;ii<3;ii++)
              for(int jj=0;jj<3;jj++)
                betaD[ii] += gammaDD[ii][jj] * betaU[jj];

            // Compute beta^2 = beta^{i} beta_{i}
            double beta2 = 0.0;
            for(int ii=0;ii<3;ii++)
              beta2 += betaU[ii] * betaD[ii];

            // Now set g_{mu nu}
            double g4DD[4][4];
            g4DD[0][0] = -alpha*alpha + beta2;
            for(int ii=0;ii<3;ii++) {
              g4DD[0][ii+1] = g4DD[ii+1][0] = betaD[ii];
              for(int jj=0;jj<3;jj++)
                g4DD[ii+1][jj+1] = gammaDD[ii][jj];
            }

            // Compute u_{mu} = g_{mu nu} u^{nu}
            double u4D[4] = {0.0, 0.0, 0.0, 0.0};
            for(int mu=0;mu<4;mu++)
              for(int nu=0;nu<4;nu++)
                u4D[mu] += g4DD[mu][nu] * u4U[nu];

            // Check constrant u_{mu} u^{mu} = -1
            double u4sqr = 0.0;
            for(int mu=0;mu<4;mu++)
              u4sqr += u4D[mu] * u4U[mu];

            if( fabs(u4sqr+1) > 1e-12 )
              CCTK_VWARN(CCTK_WARN_ALERT,
                         "Found u_{mu} u^{mu} = %.15e != -1. "
                         "alp = %e, beta^{i} = (%e, %e, %e), gamma_{ij} = (%e, %e, %e, %e, %e, %e), vel^{i} = (%e, %e, %e)",
                         u4sqr,
                         alpha, betaU[0], betaU[1], betaU[2],
                         gammaDD[0][0], gammaDD[0][1], gammaDD[0][2], gammaDD[1][1], gammaDD[1][2], gammaDD[2][2],
                         vU[0], vU[1], vU[2]);

            // Write to gridfunctions
            u4D0GF[index] = u4D[0];
            u4D1GF[index] = u4D[1];
            u4D2GF[index] = u4D[2];
            u4D3GF[index] = u4D[3];
          }
        }
      }
    }
    if(verbose>1) CCTK_VINFO("Finished computing four velocity at Ref. Lev. %d", GetRefinementLevel(cctkGH));
  }
}
