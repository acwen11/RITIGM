//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, Sebastiano Bernuzzi <sebastiano.bernuzzi@uni-jena.de>
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


#ifndef THC_M1_SOURCES_HH
#define THC_M1_SOURCES_HH

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "utils_tensor.hh"

#define THC_M1_SOURCE_OK        0
#define THC_M1_SOURCE_THIN      1
#define THC_M1_SOURCE_EDDINGTON 2
#define THC_M1_SOURCE_EQUIL     3
#define THC_M1_SOURCE_SCAT      4
#define THC_M1_SOURCE_FAIL      5

using namespace utils;

namespace thc {
namespace m1 {

// Solves the implicit problem
// .  q^new = q^star + dt S[q^new]
// The source term is S^a = (eta - ka J) u^a - (ka + ks) H^a and includes
// also emission.
int source_update(
        cGH const * cctkGH,
        int const i, int const j, int const k, int const ig,
        closure_t closure,
        gsl_root_fsolver * gsl_solver_1d,
        gsl_multiroot_fdfsolver * gsl_solver_nd,
        CCTK_REAL const cdt,
        CCTK_REAL const alp,
        tensor::metric<4> const & g_dd,
        tensor::inv_metric<4> const & g_uu,
        tensor::generic<CCTK_REAL, 4, 1> const & n_d,
        tensor::generic<CCTK_REAL, 4, 1> const & n_u,
        tensor::generic<CCTK_REAL, 4, 2> const & gamma_ud,
        tensor::generic<CCTK_REAL, 4, 1> const & u_d,
        tensor::generic<CCTK_REAL, 4, 1> const & u_u,
        tensor::generic<CCTK_REAL, 4, 1> const & v_d,
        tensor::generic<CCTK_REAL, 4, 1> const & v_u,
        tensor::generic<CCTK_REAL, 4, 2> const & proj_ud,
        CCTK_REAL const W,
        CCTK_REAL const Eold,
        tensor::generic<CCTK_REAL, 4, 1> const & Fold_d,
        CCTK_REAL const Estar,
        tensor::generic<CCTK_REAL, 4, 1> const & Fstar_d,
        CCTK_REAL const eta,
        CCTK_REAL const kabs,
        CCTK_REAL const kscat,
        CCTK_REAL * chi,
        CCTK_REAL * Enew,
        tensor::generic<CCTK_REAL, 4, 1> * Fnew_d);

} // namespace m1
} // namespace thc

#endif
