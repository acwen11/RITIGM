#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <vector>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <util_Table.h>

#include <typeinfo>

#define length(X) (sizeof(X)/sizeof(*X))
#define POW2(X) ((X)*(X))
#define POW3(X) (POW2(X)*(X))
#define SQ(X) POW2(X)

using namespace std;

namespace {
inline CCTK_REAL spatialdet(CCTK_REAL const hxx, CCTK_REAL const hxy,
        CCTK_REAL const hxz, CCTK_REAL const hyy, CCTK_REAL const hyz,
        CCTK_REAL const hzz) {
  return -SQ(hxz)*hyy + 2*hxy*hxz*hyz - hxx*SQ(hyz) - SQ(hxy)*hzz + hxx*hyy*hzz;
}
}

struct Field {
  char const * name;
  CCTK_REAL * data;
};

template <typename T> std::string type_name();

extern "C" void ZelmaniTracers_InterpBase(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( (cctk_iteration-1) % evolve_every != 0) return;

  // figure out how many tracers each process has
  int group = CCTK_GroupIndex("ZelmaniTracers::tracerevol");
  cGroupDynamicData data;
  int retval = CCTK_GroupDynamicData (cctkGH, group, &data); (void)retval;
  int myntracers = data.lsh[0];
  int myproc = CCTK_MyProc(cctkGH);


  int const interp_handle = CCTK_InterpHandle (interpolator);
  assert (interp_handle >= 0);
  int const options_handle =
    Util_TableCreateFromString (interpolator_options);
  assert (options_handle >= 0);
  int const coords_handle = CCTK_CoordSystemHandle (coordinate_system);
  assert (coords_handle >= 0);

  int npoints = myntracers;

  void const *const interp_coords[] = {
      reinterpret_cast<void *>(tx),
      reinterpret_cast<void *>(ty),
      reinterpret_cast<void *>(tz)
  };

  CCTK_INT const input_array_indices[] = {
      CCTK_VarIndex("ADMBase::alp"),
      CCTK_VarIndex("ADMBase::betax"),
      CCTK_VarIndex("ADMBase::betay"),
      CCTK_VarIndex("ADMBase::betaz"),
      CCTK_VarIndex("HydroBase::vel[0]"),
      CCTK_VarIndex("HydroBase::vel[1]"),
      CCTK_VarIndex("HydroBase::vel[2]")
  };
  int const ninputs = length(input_array_indices);

  CCTK_INT const output_array_types[] = {
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL,
      CCTK_VARIABLE_REAL
  };
  assert(ninputs == length(output_array_types));

  void * output_arrays[] = {
      reinterpret_cast<void *>(talp),
      reinterpret_cast<void *>(tbetax),
      reinterpret_cast<void *>(tbetay),
      reinterpret_cast<void *>(tbetaz),
      reinterpret_cast<void *>(tvx),
      reinterpret_cast<void *>(tvy),
      reinterpret_cast<void *>(tvz)
  };
  assert(ninputs == length(output_arrays));

  // check if tracers goes outside of the domain before interpolation
  CCTK_REAL physical_min[3], physical_max[3], interior_min[3], interior_max[3], exterior_min[3], exterior_max[3], thespacing[3];

  GetDomainSpecification(3, physical_min, physical_max, interior_min, interior_max, exterior_min, exterior_max, thespacing);

  for (int i = 0; i < myntracers; i++) {
    // ******* initialize tracers_out_of_domain, keep it to make sure tracers working correctly *******
    if (!(tracers_out_of_domain[i] ==0 || tracers_out_of_domain[i] == 1)) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "tracers_out_of_domain at %d are problematic:\n"
      "(%g, %g, %g), (alp=%g, betax=%g, betay=%g, betaz=%g, vx=%g, vy=%g, vz=%g) %d\n",
      i, tx[i], ty[i], tz[i], talp[i], tbetax[i], tbetay[i], tbetaz[i],
      tvx[i], tvy[i], tvz[i], tracers_out_of_domain[i] );
      tracers_out_of_domain[i] = 0;
    }
    // ************************************************************************************************

    // assume symmetry exist in the z direction
    if ((tx[i] >= physical_max[0] || tx[i] <= physical_min[0] ||
        ty[i] >= physical_max[1] || ty[i] <= physical_min[1] ||
        tz[i] >= physical_max[2] || tz[i] <= -physical_max[2]) && (!tracers_out_of_domain[i]) ) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "tracer %d at x=%g, y=%g, z=%g is outside of the physical domain, freeze it.\n"
      "the corresponding interpolated quantities:\n" 
      "(alp=%g, betax=%g, betay=%g, betaz=%g, vx=%g, vy=%g, vz=%g)\n",
      i, tx[i], ty[i], tz[i], talp[i], tbetax[i], tbetay[i], tbetaz[i],
      tvx[i], tvy[i], tvz[i] );
      tracers_out_of_domain[i] = 1;
    }

    // make sure interpolation is working and will not produce warnings and spam the err and out files
    if (tracers_out_of_domain[i]) {
      tx[i] = max(min(tx[i], physical_max[0]), physical_min[0]);
      ty[i] = max(min(ty[i], physical_max[1]), physical_min[1]);
      tz[i] = max(min(tz[i], physical_max[2]), -physical_max[2]);
    }
  }


  // Actual interpolation
  int const ierr =
    CCTK_InterpGridArrays (cctkGH, 3,
                           interp_handle, options_handle, coords_handle,
                           npoints, CCTK_VARIABLE_REAL, interp_coords,
                           ninputs, input_array_indices,
                           ninputs, output_array_types, output_arrays);
  // debug
  // if (ierr==CCTK_ERROR_INTERP_POINT_OUTSIDE) {
  //   for (int i = 0; i < myntracers; i++) {
  //     if (tx[i] >= physical_max[0] || tx[i] <= physical_min[0] ||
  //         ty[i] >= physical_max[1] || ty[i] <= physical_min[1] ||
  //         tz[i] >= physical_max[2] || tz[i] <= -physical_max[2] ) {
  //       CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING, "interpolation failed for tracer %d because of out-of-domain.\n"
  //       "Locate at x=%g, y=%g, z=%g, the corresponding interpolated quantities:\n" 
  //       "(alp=%g, betax=%g, betay=%g, betaz=%g, vx=%g, vy=%g, vz=%g)\n"
  //       "the physical boundaries: x(%g,%g); y(%g,%g); z(%g,%g); %d\n",
  //         i, tx[i], ty[i], tz[i], talp[i], tbetax[i], tbetay[i], tbetaz[i],
  //         tvx[i], tvy[i], tvz[i], physical_min[0], physical_max[0], physical_min[1], physical_max[1], physical_min[2], physical_max[2], 
  //         tracers_out_of_domain[i]); 
  //     }
  //   }
  // }
  // .....
  assert (ierr==0);

  Util_TableDestroy (options_handle);

  return;
}

extern "C" void ZelmaniTracers_InterpExtra(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(outTracers_every == 0) {
    if(cctk_iteration != outTracers_iteration) {
      return;
    }
  }
  else {
    if(cctk_iteration % outTracers_every != 0 && cctk_iteration != outTracers_iteration) {
      return;
    }
  }

  // figure out how many tracers each process has
  int group = CCTK_GroupIndex("ZelmaniTracers::tracerevol");
  cGroupDynamicData data;
  int retval = CCTK_GroupDynamicData (cctkGH, group, &data); (void)retval;
  int myntracers = data.lsh[0];


  int const interp_handle = CCTK_InterpHandle (interpolator);
  assert (interp_handle >= 0);
  int const options_handle =
    Util_TableCreateFromString (interpolator_options);
  assert (options_handle >= 0);
  int const coords_handle = CCTK_CoordSystemHandle (coordinate_system);
  assert (coords_handle >= 0);

  int npoints = myntracers;
  void const *const interp_coords[] = { reinterpret_cast<void *>(tx),
      reinterpret_cast<void *>(ty), reinterpret_cast<void *>(tz) };

  // Fields to interpolate in any case
  Field const fields_base[] = {
    {"HydroBase::rho",         trho},
    // {"HydroBase::eps",         teps},
    // {"HydroBase::press",       tpress},
    {"ADMBase::gxx",           tgxx},
    {"ADMBase::gxy",           tgxy},
    {"ADMBase::gxz",           tgxz},
    {"ADMBase::gyy",           tgyy},
    {"ADMBase::gyz",           tgyz},
    {"ADMBase::gzz",           tgzz},
  };
  int const nfields_base = length(fields_base);

  // Thermodynamics fields (might not have storage)
  Field const fields_tdm[] = {
    // Optional fields go here
    // {"HydroBase::entropy",     tentropy},
    {"HydroBase::temperature", ttemp},
    {"HydroBase::Y_e",         tye},
    // {"HydroBase::Abar",        tabar},
  };
  int const nfields_tdm = length(fields_tdm);

  // Leakage quantities (only if tracers_with_leakage is yes)
  Field const fields_leak[] = {
    {"THC_LeakageBase::::Q_eff_nue",  tlum0},
    {"THC_LeakageBase::::Q_eff_nua",  tlum1},
    {"THC_LeakageBase::::Q_eff_nux",  tlum2},
    {"THC_LeakageBase::::eave_local[0]", teave0},
    {"THC_LeakageBase::::eave_local[1]", teave1},
    {"THC_LeakageBase::::eave_local[2]", teave2},
  };
  int const nfields_leak = length(fields_leak);

  // Select fields to interpolate
  vector<CCTK_INT> input_array_indices;
  vector<CCTK_INT> output_array_types;
  vector<void *> output_arrays;
  for(int i = 0; i < nfields_base; ++i) {
    input_array_indices.push_back(CCTK_VarIndex(fields_base[i].name));
    output_array_types.push_back(CCTK_VARIABLE_REAL);
    output_arrays.push_back(reinterpret_cast<void *>(fields_base[i].data));
  }
  for(int i = 0; i < nfields_tdm; ++i)
  {
    if(CCTK_ActiveTimeLevelsVN(cctkGH, fields_tdm[i].name) > 0) {
      input_array_indices.push_back(CCTK_VarIndex(fields_tdm[i].name));
      output_array_types.push_back(CCTK_VARIABLE_REAL);
      output_arrays.push_back(reinterpret_cast<void *>(fields_tdm[i].data));
    }
  }
  if(tracers_with_leakage && CCTK_IsThornActive("ZelmaniLeak"))
  {
    for(int i = 0; i < nfields_leak; ++i)
    {
      input_array_indices.push_back(CCTK_VarIndex(fields_leak[i].name));
      output_array_types.push_back(CCTK_VARIABLE_REAL);
      output_arrays.push_back(reinterpret_cast<void *>(fields_leak[i].data));
    }
  }
  int const ninputs = input_array_indices.size();

  // Actual interpolation
  int const ierr =
    CCTK_InterpGridArrays (cctkGH, 3,
                           interp_handle, options_handle, coords_handle,
                           npoints, CCTK_VARIABLE_REAL, interp_coords,
                           ninputs, &input_array_indices[0],
                           ninputs, &output_array_types[0], &output_arrays[0]);
  assert (ierr==0);

  Util_TableDestroy (options_handle);

  // Compute derived quantities
  for(int i = 0; i < myntracers; ++i) {
    // CCTK_REAL const v_x = tgxx[i]*tvx[i] + tgxy[i]*tvy[i] + tgxz[i]*tvz[i];
    // CCTK_REAL const v_y = tgxy[i]*tvx[i] + tgyy[i]*tvy[i] + tgyz[i]*tvz[i];
    // CCTK_REAL const v_z = tgxz[i]*tvx[i] + tgyz[i]*tvy[i] + tgzz[i]*tvz[i];
    // CCTK_REAL const v2  = v_x*tvx[i] + v_y*tvy[i] + v_z*tvz[i];
    // tw_lorentz[i] = 1.0/sqrt(1.0 - v2);
    // teninf[i] = - tw_lorentz[i]*(v_x*tbetax[i] + v_y*tbetay[i] +
    //         v_z*tbetaz[i] - talp[i]) - 1.0;
    CCTK_REAL const hsml3 = hsml_num_neighbours * tmass[i] /
        ( 4.0/3.0 * M_PI * trho[i] );
    thsml[i] = pow(hsml3, 1.0/3.0);

    // CCTK_REAL const det = spatialdet(tgxx[i], tgxy[i], tgxz[i], tgyy[i],
    //         tgyz[i], tgzz[i]);
    // tdens[i] = trho[i]*tw_lorentz[i]*sqrt(det);
  }

  return;
}

extern "C" void ZelmaniTracers_EvolveTracers(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( (cctk_iteration-1) % evolve_every != 0) return;
  CCTK_REAL dtime = evolve_every * CCTK_DELTA_TIME;

  CCTK_REAL lapse_excision_hydro;
  // lapse_excision_hydro = *(const CCTK_REAL *) CCTK_ParameterGet("lapse_excision", "IllinoisGRMHD", NULL);
  // Parameter no longer exists, disabling for now
  lapse_excision_hydro = 0.0;

  // figure out how many tracers each process has
  int group = CCTK_GroupIndex("ZelmaniTracers::tracerevol");
  cGroupDynamicData data;
  int retval = CCTK_GroupDynamicData (cctkGH, group, &data); (void)retval;
  int myntracers = data.lsh[0];

  // rotate timelevels
  memcpy(twx_old, twx, myntracers*sizeof(CCTK_REAL));
  memcpy(twy_old, twy, myntracers*sizeof(CCTK_REAL));
  memcpy(twz_old, twz, myntracers*sizeof(CCTK_REAL));

  if(*first_step == 1 || CCTK_Equals(evolution_method, "Euler")) {
#pragma omp parallel for
    for(int i = 0; i < myntracers; ++i) {
      // initialize tracers_out_of_domain // move initialization to setup routines
      // if (*first_step == 1) tracers_out_of_domain[i]=false;

      // freeze tracers that are out of the domain
      if (tracers_out_of_domain[i]) continue;

      CCTK_REAL norm = sqrt(pow(tvx[i], 2) + pow(tvy[i], 2) + pow(tvz[i], 2));

      if ((talp[i] <= lapse_excision_hydro*(1.0+1.0e-3)) && BH_puncture) {
        if (!CCTK_isfinite(norm) ) {
          CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING, "tracer #%d entering excision region with non-finite velocity (%g,%g,%g)",
                      i, double(tvx[i]), double(tvy[i]), double(tvz[i]));
          tvx[i] = 0.0; tvy[i] = 0.0; tvz[i] = 0.0;
        }
      }
      twx[i] = talp[i]*tvx[i] - tbetax[i];
      twy[i] = talp[i]*tvy[i] - tbetay[i];
      twz[i] = talp[i]*tvz[i] - tbetaz[i];
      tx[i]  = tx[i] + dtime*twx[i];
      ty[i]  = ty[i] + dtime*twy[i];
      tz[i]  = tz[i] + dtime*twz[i];
    }
    *first_step = 0;
  }
  else {
#pragma omp parallel for
    for(int i = 0; i < myntracers; ++i) {
      // freeze tracers that are out of the domain
      if (tracers_out_of_domain[i]) continue;

      CCTK_REAL norm = sqrt(pow(tvx[i], 2) + pow(tvy[i], 2) + pow(tvz[i], 2));

      if ((talp[i] <= lapse_excision_hydro*(1.0+1.0e-3)) && BH_puncture) {
        if (!CCTK_isfinite(norm) ) {
          CCTK_VWarn(2, __LINE__, __FILE__,CCTK_THORNSTRING, "tracer #%d entering excision region with non-finite velocity (%g,%g,%g)",
                      i, double(tvx[i]), double(tvy[i]), double(tvz[i]));
          tvx[i] = 0.0; tvy[i] = 0.0; tvz[i] = 0.0;
        }
      }
      twx[i] = talp[i]*tvx[i] - tbetax[i];
      twy[i] = talp[i]*tvy[i] - tbetay[i];
      twz[i] = talp[i]*tvz[i] - tbetaz[i];
      tx[i]  = tx[i] + 0.5*dtime*(3*twx[i] - twx_old[i]);
      ty[i]  = ty[i] + 0.5*dtime*(3*twy[i] - twy_old[i]);
      tz[i]  = tz[i] + 0.5*dtime*(3*twz[i] - twz_old[i]);
    }
  }

  return;
}
