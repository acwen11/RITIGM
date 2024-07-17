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

#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_printer.hh"
#include "utils.hh"

#include "thc_M1_closure.hh"

using namespace utils;
using namespace thc::m1;

extern "C" void THC_M1_AddToTmunu(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("THC_M1_AddToTmunu");
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
            "[INFO|THC|THC_M1_AddToTmunu]: ",
            "[WARN|THC|THC_M1_AddToTmunu]: ",
            "[ERR|THC|THC_M1_AddToTmunu]: ",
            m1_max_num_msg, m1_max_num_msg);

    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, psi_bssn);
    tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz, fidu_w_lorentz,
            fidu_velx, fidu_vely, fidu_velz);

#pragma omp parallel
    {
        gsl_root_fsolver * gsl_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
        UTILS_LOOP3_DYN(thc_m1_tmunu,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            if (thc_m1_mask[ijk]) {
                continue;
            }

            tensor::metric<4> g_dd;
            tensor::inv_metric<4> g_uu;
            tensor::generic<CCTK_REAL, 4, 1> n_d;
            geom.get_metric(ijk, &g_dd);
            geom.get_inv_metric(ijk, &g_uu);
            geom.get_normal_form(ijk, &n_d);

            CCTK_REAL const W = fidu_w_lorentz[ijk];
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
            tensor::symmetric2<CCTK_REAL, 4, 2> P_dd;
            tensor::symmetric2<CCTK_REAL, 4, 2> rT_dd;

            // To de-densitize the Tmunu
            CCTK_REAL const iV = 1.0/std::pow(psi_bssn[ijk], 6);

            for (int ig = 0; ig < nspecies*ngroups; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                pack_F_d(betax[ijk], betay[ijk], betaz[ijk],
                         rFx[i4D], rFy[i4D], rFz[i4D], &F_d);

                CCTK_REAL mychi = 0.5;
                calc_closure(
                        cctkGH, i, j, k, ig,
                        closure_fun, gsl_solver, g_dd, g_uu, n_d,
                        W, u_u, v_d, proj_ud, rE[i4D], F_d,
                        &mychi, &P_dd);

                assemble_rT(n_d, rE[i4D], F_d, P_dd, &rT_dd);

                eTtt[ijk] += rT_dd(0,0)*iV;
                eTtx[ijk] += rT_dd(0,1)*iV;
                eTty[ijk] += rT_dd(0,2)*iV;
                eTtz[ijk] += rT_dd(0,3)*iV;
                eTxx[ijk] += rT_dd(1,1)*iV;
                eTxy[ijk] += rT_dd(1,2)*iV;
                eTxz[ijk] += rT_dd(1,3)*iV;
                eTyy[ijk] += rT_dd(2,2)*iV;
                eTyz[ijk] += rT_dd(2,3)*iV;
                eTzz[ijk] += rT_dd(3,3)*iV;
            }
        } UTILS_ENDLOOP3(thc_m1_tmunu);
        gsl_root_fsolver_free(gsl_solver);
    }

    // Done with printing
    thc::Printer::stop();

    // Restore GSL error handler
    gsl_set_error_handler(gsl_err);
}
