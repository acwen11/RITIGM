/********************************
 * CONVERT ET ID TO IllinoisGRMHD
 *
 * Written in 2014 by Zachariah B. Etienne
 *
 * Sets metric & MHD variables needed
 * by IllinoisGRMHD, converting from
 * HydroBase and ADMBase.
 ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "IllinoisGRMHD_headers.h"

extern "C" void set_IllinoisGRMHD_GRHD_variables_only(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  igm_eos_parameters eos;
  initialize_igm_eos_parameters_from_input(igm_eos_key,cctk_time,eos);

	CCTK_VINFO("Converting ADM -> BSSN");

  // Convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                              gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                              gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                              phi_bssn,psi_bssn,lapm1);

	CCTK_VINFO("Converting RESET HB -> IGM");

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index=CCTK_GFINDEX3D(cctkGH,i,j,k);

        if( eos.is_Hybrid ) {
					CCTK_VERROR("This feature is not supported for hybrid EOS!");
        }

        rho_b  [index]           = rho[index];
        P      [index]           = press[index];
        igm_eps[index]           = eps[index];
        if( eos.is_Tabulated ) {
          igm_Ye[index]          = Y_e[index];
          igm_temperature[index] = temperature[index];
					//CCTK_VINFO("Copied T = %e from HB to %e in IGM", temperature[index], igm_temperature[index]);
        }
        if( eos.evolve_entropy ) {
          // In this case we expect another thorn,
          // such as ID_TabEOS_HydroQuantities,
          // to ahve already taken care of the
          // entropy initialization.
          igm_entropy[index]   = entropy[index];
        }

				// IGM's A, B, and v should have been recovered and not changed
        // Ax[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
        // Ay[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
        // Az[index] = Avec[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];
        // psi6phi[index] = Aphi[index];

        // double ETvx = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,0)];
        // double ETvy = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,1)];
        // double ETvz = vel[CCTK_GFINDEX4D(cctkGH,i,j,k,2)];

        // vx[index] = alp[index]*ETvx - betax[index];
        // vy[index] = alp[index]*ETvy - betay[index];
        // vz[index] = alp[index]*ETvz - betaz[index];

      }

	CCTK_VINFO("Enforcing limits on prims in ID_converter.");
  // Finally, enforce limits on primitives & compute conservative variables.
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
        static const int zero_int=0;
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        int ww;

				// Set floor, find scaled atmospheric density and temperature.
				const CCTK_REAL r_atmo        = MAX(r_atmo_min, r[index]);
				const CCTK_REAL r_pow         = atmo_falloff ? r_power : 0.;
				const CCTK_REAL r_pow_T       = atmo_falloff_T ? r_power_T : 0.;
				const CCTK_REAL rho_b_atm     = MAX(rho_b_atm_max*std::pow(r_atmo / r_atmo_min, r_pow), eos.rho_min);
				const CCTK_REAL T_atm					= MAX(igm_T_atm*std::pow(r_atmo / r_atmo_min, r_pow_T), eos.T_min);

        CCTK_REAL PRIMS[MAXNUMVARS];
        PRIMS[RHOB         ] = rho_b[index];
        PRIMS[PRESSURE     ] = P[index];
        PRIMS[VX           ] = vx[index];
        PRIMS[VY           ] = vy[index];
        PRIMS[VZ           ] = vz[index];
        PRIMS[BX_CENTER    ] = Bx[index];
        PRIMS[BY_CENTER    ] = By[index];
        PRIMS[BZ_CENTER    ] = Bz[index];
        PRIMS[EPSILON      ] = igm_eps[index];
        if( eos.evolve_entropy ) {
          PRIMS[ENTROPY      ] = igm_entropy[index];
        }
        if( eos.is_Tabulated ) {
          PRIMS[YEPRIM     ] = igm_Ye[index];
          PRIMS[TEMPERATURE] = igm_temperature[index];
        }

        double METRIC[NUMVARS_FOR_METRIC],dummy=0;
        ww=0;
        // FIXME: NECESSARY?
        //psi_bssn[index] = exp(phi[index]);
        METRIC[ww] = phi_bssn[index];ww++;
        METRIC[ww] = dummy;          ww++; // Don't need to set psi.
        METRIC[ww] = gtxx[index];    ww++;
        METRIC[ww] = gtxy[index];    ww++;
        METRIC[ww] = gtxz[index];    ww++;
        METRIC[ww] = gtyy[index];    ww++;
        METRIC[ww] = gtyz[index];    ww++;
        METRIC[ww] = gtzz[index];    ww++;
        METRIC[ww] = lapm1[index];   ww++;
        METRIC[ww] = betax[index];   ww++;
        METRIC[ww] = betay[index];   ww++;
        METRIC[ww] = betaz[index];   ww++;
        METRIC[ww] = gtupxx[index];  ww++;
        METRIC[ww] = gtupyy[index];  ww++;
        METRIC[ww] = gtupzz[index];  ww++;
        METRIC[ww] = gtupxy[index];  ww++;
        METRIC[ww] = gtupxz[index];  ww++;
        METRIC[ww] = gtupyz[index];  ww++;

        double CONSERVS[NUM_CONSERVS] = {0,0,0,0,0};
        double g4dn[4][4];
        double g4up[4][4];
        double TUPMUNU[10],TDNMUNU[10];

        struct output_stats stats; stats.failure_checker=0;
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(zero_int,PRIMS,stats,eos,
                                                                          METRIC,g4dn,g4up,TUPMUNU,TDNMUNU,CONSERVS,r[index],rho_b_atm, T_atm);

        rho_b      [index] = PRIMS[RHOB        ];
        P          [index] = PRIMS[PRESSURE    ];
        vx         [index] = PRIMS[VX          ];
        vy         [index] = PRIMS[VY          ];
        vz         [index] = PRIMS[VZ          ];
        igm_eps    [index] = PRIMS[EPSILON     ];

        rho_star   [index] = CONSERVS[RHOSTAR  ];
        mhd_st_x   [index] = CONSERVS[STILDEX  ];
        mhd_st_y   [index] = CONSERVS[STILDEY  ];
        mhd_st_z   [index] = CONSERVS[STILDEZ  ];
        tau        [index] = CONSERVS[TAUENERGY];

        if( eos.evolve_entropy ) {
          igm_entropy[index] = PRIMS[ENTROPY     ];
          S_star     [index] = CONSERVS[ENTSTAR  ];
        }

        // Tabulated EOS
        if( eos.is_Tabulated ) {
          // Primitives
          igm_Ye[index]          = PRIMS[YEPRIM      ];
          igm_temperature[index] = PRIMS[TEMPERATURE ];
          // Conservatives
          Ye_star[index]         = CONSERVS[YESTAR   ];
        }

        if(update_Tmunu) {
          ww=0;
          eTtt[index] = TDNMUNU[ww]; ww++;
          eTtx[index] = TDNMUNU[ww]; ww++;
          eTty[index] = TDNMUNU[ww]; ww++;
          eTtz[index] = TDNMUNU[ww]; ww++;
          eTxx[index] = TDNMUNU[ww]; ww++;
          eTxy[index] = TDNMUNU[ww]; ww++;
          eTxz[index] = TDNMUNU[ww]; ww++;
          eTyy[index] = TDNMUNU[ww]; ww++;
          eTyz[index] = TDNMUNU[ww]; ww++;
          eTzz[index] = TDNMUNU[ww];
        }
      }
	CCTK_VINFO("Done.");
}
