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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_printer.hh"
#include "utils_macro.h"
#include "utils_tensor.hh"

#include "thc_M1_closure.hh"
#include "thc_M1_macro.h"
#include "thc_M1_sources.hh"

using namespace utils;
using namespace thc::m1;
using namespace std;

extern "C" void THC_M1_CalcUpdate(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_CalcUpdate");
    }

    // Disable GSL error handler
    gsl_error_handler_t * gsl_err = gsl_set_error_handler_off();

    closure_t closure_fun;
    if (CCTK_Equals(closure, "Eddington")) {
        closure_fun = eddington;
    }
    else if (CCTK_Equals(closure, "Kershaw")) {
        closure_fun = kershaw;
    }
    else if (CCTK_Equals(closure, "Minerbo")) {
        closure_fun = minerbo;
    }
    else if (CCTK_Equals(closure, "thin")) {
        closure_fun = thin;
    }
    else {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Unknown closure \"%s\"", closure);
        CCTK_ERROR(msg);
    }

    // Setup Printer
    thc::Printer::start(
            "[INFO|THC|THC_M1_CalcUpdate]: ",
            "[WARN|THC|THC_M1_CalcUpdate]: ",
            "[ERR|THC|THC_M1_CalcUpdate]: ",
            m1_max_num_msg, m1_max_num_msg);

    // Steps
    // 1. F^m   = F^k + dt/2 [ A[F^k] + S[F^m]   ]
    // 2. F^k+1 = F^k + dt   [ A[F^m] + S[F^k+1] ]
    // At each step we solve an implicit problem in the form
    //    F = F^* + cdt S[F]
    // Where F^* = F^k + cdt A
    CCTK_REAL const dt = CCTK_DELTA_TIME / static_cast<CCTK_REAL>(
            *TimeIntegratorStage);
    --(*TimeIntegratorStage);

    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);
    tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz, fidu_w_lorentz,
            fidu_velx, fidu_vely, fidu_velz);

    int const siz = UTILS_GFSIZE(cctkGH);
    CCTK_REAL * sconx = &scon[0*siz];
    CCTK_REAL * scony = &scon[1*siz];
    CCTK_REAL * sconz = &scon[2*siz];

    CCTK_REAL const mb = AverageBaryonMass();

#pragma omp parallel
    {
        gsl_root_fsolver * gsl_solver_1d =
            gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        gsl_multiroot_fdfsolver * gsl_solver_nd =
            gsl_multiroot_fdfsolver_alloc(gsl_multiroot_fdfsolver_hybridsj, 4);
        UTILS_LOOP3_DYN(thc_m1_calc_update,
                k, THC_M1_NGHOST, cctk_lsh[2]-THC_M1_NGHOST,
                j, THC_M1_NGHOST, cctk_lsh[1]-THC_M1_NGHOST,
                i, THC_M1_NGHOST, cctk_lsh[0]-THC_M1_NGHOST) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            netabs[ijk] = 0;
            netheat[ijk] = 0;
            if (thc_m1_mask[ijk]) {
                continue;
            }

            tensor::metric<4> g_dd;
            tensor::inv_metric<4> g_uu;
            tensor::generic<CCTK_REAL, 4, 1> n_u;
            tensor::generic<CCTK_REAL, 4, 1> n_d;
            tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
            geom.get_metric(ijk, &g_dd);
            geom.get_inv_metric(ijk, &g_uu);
            geom.get_normal(ijk, &n_u);
            geom.get_normal_form(ijk, &n_d);
            geom.get_space_proj(ijk, &gamma_ud);

            tensor::generic<CCTK_REAL, 4, 1> u_u;
            tensor::generic<CCTK_REAL, 4, 1> u_d;
            tensor::generic<CCTK_REAL, 4, 2> proj_ud;
            fidu.get(ijk, &u_u);
            tensor::contract(g_dd, u_u, &u_d);
            calc_proj(u_d, u_u, &proj_ud);

            tensor::generic<CCTK_REAL, 4, 1> v_u;
            tensor::generic<CCTK_REAL, 4, 1> v_d;
            pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
            tensor::contract(g_dd, v_u, &v_d);

            tensor::generic<CCTK_REAL, 4, 1> F_d;
            tensor::generic<CCTK_REAL, 4, 1> Fstar_d;
            tensor::generic<CCTK_REAL, 4, 1> Fnew_d;

            //
            // Source RHS are stored here
            CCTK_REAL DrE[ngroups*nspecies];
            CCTK_REAL DrFx[ngroups*nspecies];
            CCTK_REAL DrFy[ngroups*nspecies];
            CCTK_REAL DrFz[ngroups*nspecies];
            CCTK_REAL DrN[ngroups*nspecies];
            CCTK_REAL DDxp[ngroups*nspecies];

            //
            // Step 1 -- compute the sources
            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

#if (THC_M1_SRC_METHOD == THC_M1_SRC_EXPL)
#warning "THC_M1: using explicit collisional solver!"
                //
                // Radiation fields
                CCTK_REAL E = rE[i4D];
                pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                         rFx[i4D], rFy[i4D], rFz[i4D], &F_d);
                tensor::generic<CCTK_REAL, 4, 1> F_u;
                tensor::generic<CCTK_REAL, 4, 1> S_d;
                tensor::generic<CCTK_REAL, 4, 1> tS_d;
                tensor::contract(g_uu, F_d, &F_u);

                //
                // Compute radiation quantities in the fluid frame
                CCTK_REAL J = rJ[i4D];
                CCTK_REAL const Gamma = compute_Gamma(
                        fidu_w_lorentz[ijk], v_u, J, E, F_d);
                tensor::generic<CCTK_REAL, 4, 1> H_d;
                pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

                //
                // Compute radiation sources
                calc_rad_sources(eta_1[i4D]*volform[ijk],
                        abs_1[i4D], scat_1[i4D], u_d, J, H_d, &S_d);
                DrE[ig] = dt*calc_rE_source(alp[ijk], n_u, S_d);

                calc_rF_source(alp[ijk], gamma_ud, S_d, &tS_d);
                DrFx[ig] = dt*tS_d(1);
                DrFy[ig] = dt*tS_d(2);
                DrFz[ig] = dt*tS_d(3);

                DrN[ig] = dt*alp[ijk]*(volform[ijk]*eta_0[i4D] - abs_0[i4D]*rN[i4D]/Gamma);

#else // (THC_M1_SRC_METHOD == THC_M1_SRC_EXPL)

                //
                // Here we boost to the fluid frame, compute fluid matter
                // interaction, and boost back. These values are used as
                // initial guess for the implicit solve.

                //
                // Advect radiation
                CCTK_REAL Estar = rE_p[i4D] + dt*rE_rhs[i4D];
                pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                        rFx_p[i4D] + dt*rFx_rhs[i4D],
                        rFy_p[i4D] + dt*rFy_rhs[i4D],
                        rFz_p[i4D] + dt*rFz_rhs[i4D],
                        &Fstar_d);
                apply_floor(g_uu, &Estar, &Fstar_d);
                CCTK_REAL Nstar = max(rN_p[i4D] + dt*rN_rhs[i4D], rad_N_floor);
                CCTK_REAL Enew;

                //
                // Compute quantities in the fluid frame
                tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
                calc_closure(cctkGH, i, j, k, ig,
                        closure_fun, gsl_solver_1d, g_dd, g_uu, n_d,
                        fidu_w_lorentz[ijk], u_u, v_d, proj_ud, Estar, Fstar_d,
                        &chi[i4D], &P_dd);

                tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;
                assemble_rT(n_d, Estar, Fstar_d, P_dd, &rT_dd);

                CCTK_REAL const Jstar = calc_J_from_rT(rT_dd, u_u);

                tensor::generic<CCTK_REAL, 4, 1> Hstar_d;
                calc_H_from_rT(rT_dd, u_u, proj_ud, &Hstar_d);

                //
                // Estimate interaction with matter
                CCTK_REAL const dtau = dt/fidu_w_lorentz[ijk];
                CCTK_REAL Jnew = (Jstar + dtau*eta_1[i4D]*volform[ijk])/(1 + dtau*abs_1[i4D]);

                // Only three components of H^a are independent H^0 is found by
                // requiring H^a u_a = 0
                CCTK_REAL const khat = (abs_1[i4D] + scat_1[i4D]);
                tensor::generic<CCTK_REAL, 4, 1> Hnew_d;
                for (int a = 1; a < 4; ++a) {
                    Hnew_d(a) = Hstar_d(a)/(1 + dtau*khat);
                }
                Hnew_d(0) = 0.0;
                for (int a = 1; a < 4; ++a) {
                    Hnew_d(0) -= Hnew_d(a)*(u_u(a)/u_u(0));
                }

                //
                // Update Tmunu
                CCTK_REAL const H2 = tensor::dot(g_uu, Hnew_d, Hnew_d);
#if (THC_M1_SRC_METHOD == THC_M1_SRC_BOOST)
                CCTK_REAL const xi = sqrt(H2)*(Jnew > rad_E_floor ? 1/Jnew : 0);
                chi[i4D] = closure_fun(xi);
#else // this is really only important in the thick limit, so take chi = 1/3
                chi[i4D] = 1.0/3.0;
#endif
                CCTK_REAL const dthick = 3.*(1. - chi[i4D])/2.;
                CCTK_REAL const dthin = 1. - dthick;

                for(int a = 0; a < 4; ++a)
                for(int b = a; b < 4; ++b) {
                    rT_dd(a,b) = Jnew*u_d(a)*u_d(b) + Hnew_d(a)*u_d(b) + Hnew_d(b)*u_d(a) +
                        dthin*Jnew*(Hnew_d(a)*Hnew_d(b)*(H2 > 0 ? 1/H2 : 0)) +
                        dthick*Jnew*(g_dd(a,b) + u_d(a)*u_d(b))/3;
                }

                //
                // Boost back to the lab frame
                Enew = calc_J_from_rT(rT_dd, n_u);
                calc_H_from_rT(rT_dd, n_u, gamma_ud, &Fnew_d);
                apply_floor(g_uu, &Enew, &Fnew_d);
#if (THC_M1_SRC_METHOD == THC_M1_SRC_IMPL)
                //
                // Compute interaction with matter
                (void)source_update(cctkGH, i, j, k, ig,
                        closure_fun, gsl_solver_1d, gsl_solver_nd, dt,
                        alp[ijk], g_dd, g_uu, n_d, n_u, gamma_ud, u_d, u_u,
                        v_d, v_u, proj_ud, fidu_w_lorentz[ijk], Estar, Fstar_d,
                        Estar, Fstar_d, volform[ijk]*eta_1[i4D],
                        abs_1[i4D], scat_1[i4D], &chi[i4D], &Enew, &Fnew_d);
                apply_floor(g_uu, &Enew, &Fnew_d);

                //
                // Update closure
                apply_closure(g_dd, g_uu, n_d, fidu_w_lorentz[ijk],
                        u_u, v_d, proj_ud, Enew, Fnew_d, chi[i4D], &P_dd);

                //
                // Compute new radiation energy density in the fluid frame
                tensor::symmetric2<CCTK_REAL, 4, 2> T_dd;
                assemble_rT(n_d, Enew, Fnew_d, P_dd, &T_dd);
                Jnew = calc_J_from_rT(T_dd, u_u);
#endif // (THC_M1_SRC_METHOD == THC_M1_SRC_IMPL)

                //
                // Compute changes in radiation energy and momentum
                DrE[ig]  = Enew - Estar;
                DrFx[ig] = Fnew_d(1) - Fstar_d(1);
                DrFy[ig] = Fnew_d(2) - Fstar_d(2);
                DrFz[ig] = Fnew_d(3) - Fstar_d(3);

                //
                // Compute updated Gamma
                CCTK_REAL const Gamma = compute_Gamma(
                        fidu_w_lorentz[ijk], v_u, Jnew, Enew, Fnew_d);

                //
                // N^k+1 = N^* + dt ( eta - abs N^k+1 )
                if (source_therm_limit < 0 || dt*abs_0[i4D] < source_therm_limit) {
                    DrN[ig] = (Nstar + dt*alp[ijk]*volform[ijk]*eta_0[i4D])/
                                (1 + dt*alp[ijk]*abs_0[i4D]/Gamma) - Nstar;
                }
                //
                // The neutrino number density is updated assuming the neutrino
                // average energies are those of the equilibrium
                else {
                    DrN[ig] = (nueave[i4D] > 0 ? Gamma*Jnew/nueave[i4D] - Nstar : 0.0);
                }
#endif // (THC_M1_SRC_METHOD == THC_M1_SRC_EXPL)

                //
                // Fluid lepton sources
                DDxp[ig] = -mb*(DrN[ig]*(ig == 0) - DrN[ig]*(ig == 1));
            }

            //
            // Step 2 -- limit the sources
            CCTK_REAL theta = 1.0;
            if (source_limiter >= 0) {
                theta = 1.0;
                CCTK_REAL DTau_sum = 0.0;
                CCTK_REAL DDxp_sum = 0.0;
                for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                    int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                    CCTK_REAL Estar = rE_p[i4D] + dt*rE_rhs[i4D];
                    if (DrE[ig] < 0) {
                        theta = min(-source_limiter*max(Estar, 0.0)/DrE[ig], theta);
                    }
                    DTau_sum -= DrE[ig];

                    CCTK_REAL Nstar = rN_p[i4D] + dt*rN_rhs[i4D];
                    if (DrN[ig] < 0) {
                        theta = min(-source_limiter*max(Nstar, 0.0)/DrN[ig], theta);
                    }
                    DDxp_sum += DDxp[ig];
                }
                CCTK_REAL const DYe = DDxp_sum/dens[ijk];
                if (DTau_sum < 0) {
                    theta = min(-source_limiter*max(tau[ijk], 0.0)/DTau_sum, theta);
                }
                if (DYe > 0) {
                    theta = min(source_limiter*max(source_Ye_max - Y_e[ijk], 0.0)/DYe, theta);
                }
                else if (DYe < 0) {
                    theta = min(source_limiter*min(source_Ye_min - Y_e[ijk], 0.0)/DYe, theta);
                }
                theta = max(0.0, theta);
            }

            //
            // Step 3 -- update fields
            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                //
                // Update radiation quantities
                CCTK_REAL E =  rE_p[i4D] + dt*rE_rhs[i4D]  + theta*DrE[ig];
                F_d(1)      = rFx_p[i4D] + dt*rFx_rhs[i4D] + theta*DrFx[ig];
                F_d(2)      = rFy_p[i4D] + dt*rFy_rhs[i4D] + theta*DrFy[ig];
                F_d(3)      = rFz_p[i4D] + dt*rFz_rhs[i4D] + theta*DrFz[ig];
                apply_floor(g_uu, &E, &F_d);

                CCTK_REAL N =  rN_p[i4D] + dt*rN_rhs[i4D]  + theta*DrN[ig];
                N = max(N, rad_N_floor);

                //
                // Compute back reaction on the fluid
                // NOTE: fluid backreaction is only needed at the last substep
                if (backreact && 0 == *TimeIntegratorStage) {
                    assert (ngroups == 1);
                    assert (nspecies == 3);

                    sconx[ijk]  -= theta*DrFx[ig];
                    scony[ijk]  -= theta*DrFy[ig];
                    sconz[ijk]  -= theta*DrFz[ig];
                    tau[ijk]    -= theta*DrE[ig];
                    densxp[ijk] += theta*DDxp[ig];
                    densxn[ijk] -= theta*DDxp[ig];

                    netabs[ijk]  += theta*DDxp[ig];
                    netheat[ijk] -= theta*DrE[ig];
                }

                //
                // Save updated results into grid functions
                rE[i4D]  = E;
                unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
                rN[i4D] = N;
            }
        } UTILS_ENDLOOP3(thc_m1_calc_update);
        gsl_root_fsolver_free(gsl_solver_1d);
        gsl_multiroot_fdfsolver_free(gsl_solver_nd);
    }

    // Done with printing
    thc::Printer::stop();

    // Restore GSL error handler
    gsl_set_error_handler(gsl_err);
}
