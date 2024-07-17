//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include <cassert>
#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"
#include "utils.hh"

using namespace thc::m1;
using namespace std;
using namespace utils;

#define NDIM 3

#define GFINDEX1D(__k, ig, iv) \
    ((iv) + (ig)*5 + (__k)*(5*ngroups*nspecies))

#define PINDEX1D(ig, iv) \
    ((iv) + (ig)*5)

namespace {

CCTK_REAL minmod2(CCTK_REAL rl, CCTK_REAL rp, CCTK_REAL th) {
    return min(1.0, min(th*rl, th*rp));
}

}

// Compute the numerical fluxes using a simple 2nd order flux-limited method
// with high-Peclet limit fix. The fluxes are then added to the RHSs.
extern "C" void THC_M1_CalcFluxes(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_CalcFluxes");
    }

    // Aliases for the RHS pointers
    CCTK_REAL ** rhs = new CCTK_REAL * [ngroups*nspecies*5];
    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
        int const i4D = CCTK_VectGFIndex3D(cctkGH, 0, 0, 0, ig);
        rhs[PINDEX1D(ig, 0)] = &rN_rhs[i4D];
        rhs[PINDEX1D(ig, 1)] = &rFx_rhs[i4D];
        rhs[PINDEX1D(ig, 2)] = &rFy_rhs[i4D];
        rhs[PINDEX1D(ig, 3)] = &rFz_rhs[i4D];
        rhs[PINDEX1D(ig, 4)] = &rE_rhs[i4D];
    }

    // Grid data
    CCTK_REAL const delta[3] = {
        CCTK_DELTA_SPACE(0),
        CCTK_DELTA_SPACE(1),
        CCTK_DELTA_SPACE(2),
    };
    CCTK_REAL const idelta[3] = {
        1.0/CCTK_DELTA_SPACE(0),
        1.0/CCTK_DELTA_SPACE(1),
        1.0/CCTK_DELTA_SPACE(2),
    };

    // Tensor fields
    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, psi_bssn);
    tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz, fidu_w_lorentz,
            fidu_velx, fidu_vely, fidu_velz);

#pragma omp parallel
    {
        for (int dir = 0; dir < NDIM; ++dir) {
            int index[3];
            int lsh[3];

            // We have to leave as the most internal loop the one on the
            // direction. For this reason we will remap the usual indices
            // i,j,k into different points of index[:].
            int ii, ij, ik;
            switch(dir) {
                case 0:
                    ii = 2;
                    ij = 1;
                    ik = 0;

                    lsh[0] = cctk_lsh[2];
                    lsh[1] = cctk_lsh[1];
                    lsh[2] = cctk_lsh[0];

                    break;
                case 1:
                    ii = 1;
                    ij = 2;
                    ik = 0;

                    lsh[0] = cctk_lsh[2];
                    lsh[1] = cctk_lsh[0];
                    lsh[2] = cctk_lsh[1];

                    break;
                case 2:
                    ii = 1;
                    ij = 0;
                    ik = 2;

                    lsh[0] = cctk_lsh[1];
                    lsh[1] = cctk_lsh[0];
                    lsh[2] = cctk_lsh[2];

                    break;
            }

            // Indices aliases
            int & i = index[ii];
            int & j = index[ij];
            int & k = index[ik];

            // Actual indices
            int __i, __j, __k;

            // Scratch space
            CCTK_REAL * cons;
            CCTK_REAL * flux;
            CCTK_REAL * cmax    = NULL;
            CCTK_REAL * flux_jm = NULL;
            CCTK_REAL * flux_jp = NULL;
            CCTK_REAL * d_ptr   = NULL;
            try {
                cons = new CCTK_REAL[5*ngroups*nspecies*lsh[2]];
                flux = new CCTK_REAL[5*ngroups*nspecies*lsh[2]];
                cmax = new CCTK_REAL[5*ngroups*nspecies*lsh[2]];
                flux_jm = new CCTK_REAL[5*ngroups*nspecies];
                flux_jp = new CCTK_REAL[5*ngroups*nspecies];
            }
            catch (std::bad_alloc &e) {
#pragma omp critical
                CCTK_ERROR("Out of memory!");
            }

#pragma omp for collapse(2)
            for (__i = THC_M1_NGHOST; __i < lsh[0] - THC_M1_NGHOST; ++__i)
            for (__j = THC_M1_NGHOST; __j < lsh[1] - THC_M1_NGHOST; ++__j) {
                // ----------------------------------------------
                // 1st pass compute the fluxes
                // . this could be made into a separate Cactus scheduled
                //   function to avoid some duplicated calculations
                for (__k = 0; __k < lsh[2]; ++__k) {
                    index[0] = __i;
                    index[1] = __j;
                    index[2] = __k;
                    int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                    tensor::inv_metric<3> gamma_uu;
                    tensor::inv_metric<4> g_uu;
                    tensor::generic<CCTK_REAL, 4, 1> beta_u;

                    geom.get_inv_metric(ijk, &gamma_uu);
                    geom.get_inv_metric(ijk, &g_uu);
                    geom.get_shift_vec(ijk, &beta_u);

                    tensor::generic<CCTK_REAL, 4, 1> u_u;
                    fidu.get(ijk, &u_u);

                    tensor::generic<CCTK_REAL, 4, 1> v_u;
                    pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);

                    tensor::generic<CCTK_REAL, 4, 1> H_d;
                    tensor::generic<CCTK_REAL, 4, 1> H_u;
                    tensor::generic<CCTK_REAL, 4, 1> F_d;
                    tensor::generic<CCTK_REAL, 4, 1> F_u;
                    tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
                    tensor::generic<CCTK_REAL, 4, 2> P_ud;

                    tensor::generic<CCTK_REAL, 4, 1> fnu_u;

                    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                        int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                        pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                                 rFx[i4D], rFy[i4D], rFz[i4D], &F_d);
                        pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);
                        pack_P_dd(betax[ijk], betay[ijk], betaz[ijk],
                                  rPxx[i4D], rPxy[i4D], rPxz[i4D],
                                  rPyy[i4D], rPyz[i4D], rPzz[i4D], &P_dd);
                        tensor::contract(g_uu, H_d, &H_u);
                        tensor::contract(g_uu, F_d, &F_u);
                        tensor::contract(g_uu, P_dd, &P_ud);

                        assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u);
                        CCTK_REAL const Gamma = compute_Gamma(
                                fidu_w_lorentz[ijk], v_u, rJ[i4D], rE[i4D], F_d);
                        // Note that nnu is densitized here
                        CCTK_REAL const nnu = rN[i4D]/Gamma;

                        cons[GFINDEX1D(__k, ig, 0)] = rN[i4D];
                        cons[GFINDEX1D(__k, ig, 1)] = rFx[i4D];
                        cons[GFINDEX1D(__k, ig, 2)] = rFy[i4D];
                        cons[GFINDEX1D(__k, ig, 3)] = rFz[i4D];
                        cons[GFINDEX1D(__k, ig, 4)] = rE[i4D];

                        assert(isfinite(cons[GFINDEX1D(__k, ig, 0)]));
                        assert(isfinite(cons[GFINDEX1D(__k, ig, 1)]));
                        assert(isfinite(cons[GFINDEX1D(__k, ig, 2)]));
                        assert(isfinite(cons[GFINDEX1D(__k, ig, 3)]));
                        assert(isfinite(cons[GFINDEX1D(__k, ig, 4)]));

                        flux[GFINDEX1D(__k, ig, 0)] =
                            alp[ijk] * nnu * fnu_u(dir+1);
                        flux[GFINDEX1D(__k, ig, 1)] =
                            calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir+1, 1);
                        flux[GFINDEX1D(__k, ig, 2)] =
                            calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir+1, 2);
                        flux[GFINDEX1D(__k, ig, 3)] =
                            calc_F_flux(alp[ijk], beta_u, F_d, P_ud, dir+1, 3);
                        flux[GFINDEX1D(__k, ig, 4)] =
                            calc_E_flux(alp[ijk], beta_u, rE[i4D], F_u, dir+1);

                        assert(isfinite(flux[GFINDEX1D(__k, ig, 0)]));
                        assert(isfinite(flux[GFINDEX1D(__k, ig, 1)]));
                        assert(isfinite(flux[GFINDEX1D(__k, ig, 2)]));
                        assert(isfinite(flux[GFINDEX1D(__k, ig, 3)]));
                        assert(isfinite(flux[GFINDEX1D(__k, ig, 4)]));

                        // Eigenvalues in the optically thin limit
                        //
                        // Using the optically thin eigenvalues seems to cause
                        // problems in some situations, possibly because the
                        // optically thin closure is acausal in certain
                        // conditions, so this is commented for now.
#if 0
                        CCTK_REAL const F2 = tensor::dot(F_u, F_d);
                        CCTK_REAL const F = std::sqrt(F2);
                        CCTK_REAL const fx = F_u(dir+1)*(F > 0 ? 1/F : 0);
                        CCTK_REAL const ffx = F_u(dir+1)*(F2 > 0 ? 1/F2 : 0);
                        CCTK_REAL const lam[3] = {
                              alp[ijk]*fx - beta_u(dir+1),
                            - alp[ijk]*fx - beta_u(dir+1),
                              alp[ijk]*rE[i4D]*ffx - beta_u(dir+1),
                        };
                        CCTK_REAL const cM1 = std::max(std::abs(lam[0]),
                                std::max(std::abs(lam[1]), std::abs(lam[2])));

                        // TODO optically thick characteristic speeds and combination
#endif
                        // Speed of light -- note that gamma_uu has NDIM=3
                        CCTK_REAL const clam[2] = {
                              alp[ijk]*std::sqrt(gamma_uu(dir,dir)) - beta_u(dir+1),
                            - alp[ijk]*std::sqrt(gamma_uu(dir,dir)) - beta_u(dir+1),
                        };
                        CCTK_REAL const clight = std::max(
                                std::abs(clam[0]), std::abs(clam[1]));

                        // In some cases the eigenvalues in the thin limit
                        // overestimate the actual characteristic speed and can
                        // become larger than c
                        cmax[GFINDEX1D(__k, ig, 0)] = clight;
                                                 // = std::min(clight, cM1);
                    }
                }

                // Cleanup the flux buffer
                memset(flux_jm, 0, 5*ngroups*nspecies*sizeof(CCTK_REAL));
                memset(flux_jp, 0, 5*ngroups*nspecies*sizeof(CCTK_REAL));

                // ----------------------------------------------
                // 2nd pass update the RHS
                for (__k = THC_M1_NGHOST-1; __k < lsh[2]-THC_M1_NGHOST; ++__k) {
                    index[0] = __i;
                    index[1] = __j;
                    index[2] = __k;
                    int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                    for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                        int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                        index[2]++;
                        int const i4Dp = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
                        index[2]--;

                        // Remove dissipation at high Peclet numbers
                        CCTK_REAL kapa = 0.5*(abs_1[i4D] + abs_1[i4Dp] +
                                scat_1[i4D] + scat_1[i4Dp]);
                        CCTK_REAL A = 1.0;
                        if (kapa*delta[dir] > 1.0) {
                            A = min(1.0, 1.0/(delta[dir]*kapa));
                            A = max(A, mindiss);
                        }

                        for (int iv = 0; iv < 5; ++iv) {
                            CCTK_REAL const ujm = cons[GFINDEX1D(__k-1, ig, iv)];
                            CCTK_REAL const uj = cons[GFINDEX1D(__k, ig, iv)];
                            CCTK_REAL const ujp = cons[GFINDEX1D(__k+1, ig, iv)];
                            CCTK_REAL const ujpp = cons[GFINDEX1D(__k+2, ig, iv)];

                            CCTK_REAL const fj = flux[GFINDEX1D(__k, ig, iv)];
                            CCTK_REAL const fjp = flux[GFINDEX1D(__k+1, ig, iv)];

                            CCTK_REAL const cc = cmax[GFINDEX1D(__k, ig, 0)];
                            CCTK_REAL const ccp = cmax[GFINDEX1D(__k+1, ig, 0)];
                            CCTK_REAL const cmx = std::max(cc, ccp);

                            CCTK_REAL const dup = ujpp - ujp;
                            CCTK_REAL const duc = ujp - uj;
                            CCTK_REAL const dum = uj - ujm;

                            bool sawtooth = false;
                            CCTK_REAL phi = 0;
                            if (dup*duc > 0 && dum*duc > 0) {
                                phi = minmod2(dum/duc, dup/duc, minmod_theta);
                            }
                            else if (dup*duc < 0 && dum*duc < 0) {
                                sawtooth = true;
                            }
                            assert(isfinite(phi));

                            CCTK_REAL const flux_low =
                                0.5*(fj + fjp - cmx*(ujp - uj));
                            CCTK_REAL const flux_high = 0.5*(fj + fjp);

                            CCTK_REAL flux_num = flux_high -
                                (sawtooth ? 1 : A)*(1 - phi)*(flux_high - flux_low);
                            flux_jp[PINDEX1D(ig, iv)] = flux_num;

                            if (!thc_m1_mask[ijk]) {
                                rhs[PINDEX1D(ig, iv)][ijk] += idelta[dir]*(
                                   flux_jm[PINDEX1D(ig, iv)] -
                                   flux_jp[PINDEX1D(ig, iv)])*
                                       static_cast<CCTK_REAL>(
                                          i >= THC_M1_NGHOST
                                       && i <  cctk_lsh[0] - THC_M1_NGHOST
                                       && j >= THC_M1_NGHOST
                                       && j <  cctk_lsh[1] - THC_M1_NGHOST
                                       && k >= THC_M1_NGHOST
                                       && k <  cctk_lsh[2] - THC_M1_NGHOST);
                                assert(isfinite(rhs[PINDEX1D(ig, iv)][ijk]));
                            }
                        }
                    }
                    // Rotate flux pointer
                    d_ptr = flux_jm;
                    flux_jm = flux_jp;
                    flux_jp = d_ptr;
                }
            }

            delete[] cons;
            delete[] flux;
            delete[] cmax;
            delete[] flux_jm;
            delete[] flux_jp;
        }
    }

    delete[] rhs;
}
