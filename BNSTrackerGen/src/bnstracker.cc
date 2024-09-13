#include "bnstracker.h"
#include <cassert>
#include <cmath>

using namespace Pizza;
using namespace Pizza::BNSTrackerGen;

bns_tracker::bns_tracker(int use_reflevel_, bool pi_sym_, int dist_num_points_)
: use_reflevel(use_reflevel_), pi_sym(pi_sym_), 
  dist_num_points(dist_num_points_)
{}

pz_real bns_tracker::sym_weight(vec_u& pos) const
{
  assert(pi_sym);
  pz_real s = pos(0) / dw(0);
  if (s<= -0.5) return 0.0;
  else if (s<0.5) return (2.0*s + 1.0);
  else return 2.0;
}

void bns_tracker::get_dens(const gpos& i, pz_real& rmd, pz_real& crmd)
{
  mats_l glo;
  glo_pc(i) >> glo;
  pz_real vc  = sqrt(det(glo));
  rmd         = rmd_pc(i);
  crmd        = vc * lfac_pc(i) * rmd;
}

void bns_tracker::track_individual_pisym(pz_real phase_est)
{
  assert(pi_sym);
  pz_cmplx  rotate = exp(pz_cmplx(0,-phase_est));
  pz_cmplx avg_pos_1(0), avg_pos_2(0), pos_max_1(0), pos_max_2(0);
  pz_real wsum_1(0), maxrmd_1(-1);
  for (region::iterator i(r_int());i;++i) {
    vec_u pos   = coord(i);
    pz_real weight, rmd; 
    get_dens(i, rmd, weight);
    pz_cmplx cp(pos(0), pos(1));
    cp         *= rotate;
    if (cp.real() <= 0) {
      cp = -cp;
    }
    avg_pos_1 += weight * cp;
    wsum_1    += weight;
    if (rmd > maxrmd_1) {
      maxrmd_1  = rmd;
      pos_max_1 = cp;
    }
  }
  avg_pos_1      = global_sum(avg_pos_1);
  wsum_1         = global_sum(wsum_1);

  avg_pos_1      = (avg_pos_1 / wsum_1) / rotate;
  avg_pos_2      = -avg_pos_1;

  pz_real mcnt_1 = (maxrmd_1 < global_max(maxrmd_1)) ? 0.0 : 1.0;
  pos_max_1      = pos_max_1 * mcnt_1;  
  mcnt_1         = global_sum(mcnt_1);
  pos_max_1      = global_sum(pos_max_1);
  pos_max_1      = (pos_max_1 / mcnt_1) / rotate;
  pos_max_2      = -pos_max_1;


  star1.set_phase_sep(remove_phase_jump(phase_est, 
                        std::arg(avg_pos_1)), abs(avg_pos_1));
  star2.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                        std::arg(avg_pos_2)), abs(avg_pos_2));

  star1_maxloc.set_phase_sep(remove_phase_jump(phase_est, 
                               std::arg(pos_max_1)), abs(pos_max_1));
  star2_maxloc.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                               std::arg(pos_max_2)), abs(pos_max_2));
}

void bns_tracker::track_individual_nosym(pz_real phase_est)
{
  assert(!pi_sym);
  pz_cmplx  rotate = exp(pz_cmplx(0,-phase_est));
  pz_cmplx avg_pos_1(0), avg_pos_2(0), pos_max_1(0), pos_max_2(0);
  pz_real wsum_1(0), wsum_2(0), maxrmd_1(-1), maxrmd_2(-1);
  for (region::iterator i(r_int());i;++i) {
    vec_u pos   = coord(i);
    pz_real weight, rmd; 
    get_dens(i, rmd, weight);
    pz_cmplx cp(pos(0)-cms_x, pos(1)-cms_y);
    cp         *= rotate;
    if (cp.real() > 0) {
      avg_pos_1 += weight * cp;
      wsum_1    += weight;
      if (rmd > maxrmd_1) {
        maxrmd_1  = rmd;
        pos_max_1 = cp;
      }
    }
    else {
      avg_pos_2 += weight * cp;
      wsum_2    += weight;
      if (rmd > maxrmd_2) {
        maxrmd_2  = rmd;
        pos_max_2 = cp;
      }
    }
  }
  avg_pos_1      = global_sum(avg_pos_1);
  wsum_1         = global_sum(wsum_1);
  avg_pos_2      = global_sum(avg_pos_2);
  wsum_2         = global_sum(wsum_2);
  avg_pos_1      = (avg_pos_1 / wsum_1) / rotate;
  avg_pos_2      = (avg_pos_2 / wsum_2) / rotate;

  pz_real mcnt_1 = (maxrmd_1 < global_max(maxrmd_1)) ? 0.0 : 1.0;
  pos_max_1      = pos_max_1 * mcnt_1;
  mcnt_1         = global_sum(mcnt_1);
  pos_max_1      = global_sum(pos_max_1);
  pos_max_1      = (pos_max_1 / mcnt_1) / rotate;

  pz_real mcnt_2 = (maxrmd_2 < global_max(maxrmd_2)) ? 0.0 : 1.0;
  pos_max_2      = pos_max_2 * mcnt_2;
  mcnt_2         = global_sum(mcnt_2);
  pos_max_2      = global_sum(pos_max_2);
  pos_max_2      = (pos_max_2 / mcnt_2) / rotate;

  star1.set_phase_sep(remove_phase_jump(phase_est, 
                        std::arg(avg_pos_1)), abs(avg_pos_1));
  star2.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                        std::arg(avg_pos_2)), abs(avg_pos_2));

  star1_maxloc.set_phase_sep(remove_phase_jump(phase_est, 
                               std::arg(pos_max_1)), abs(pos_max_1));
  star2_maxloc.set_phase_sep(remove_phase_jump(phase_est + M_PI, 
                               std::arg(pos_max_2)), abs(pos_max_2));
}

void bns_tracker::track_average_pisym(const double prev_phase)
{
  assert(pi_sym);
  pz_cmplx avg_pos(0);
  pz_real avg_r(0), wsum(0);
  for (region::iterator i(r_int());i;++i) {
    vec_u pos   = coord(i);
    pz_real rmd, crmd;
    get_dens(i, rmd, crmd);
    pz_real weight = crmd * sym_weight(pos);
    pz_cmplx cp(pos(0), pos(1));
    pz_real d      = abs(cp);
    pz_cmplx cp2   = d > 0 ? (cp*(cp/d)) : 0.0;
    avg_r      += weight * d;
    avg_pos    += weight * cp2;
    wsum       += weight;
  }
  avg_r           = global_sum(avg_r);
  avg_pos         = global_sum(avg_pos);
  wsum            = global_sum(wsum);
  avg_r          /= wsum;
  avg_pos        /= wsum;
  pz_real new_phase  = std::arg(avg_pos) / 2.0;
  pz_real new_sep    = avg_r;
 
  avg.set_phase_sep(remove_phase_jump(prev_phase*2.0, 
                      new_phase*2.0) / 2.0 , new_sep);  
  cms_x = 0.0;
  cms_y = 0.0;
  cms_z = 0.0; //TODO: track z-movement also for pi-symmetric case
}

void bns_tracker::track_average_nosym(const double prev_phase)
{
  assert(!pi_sym);
  pz_cmplx avg_pos(0), cms_pos(0);
  pz_real cms_zpos(0);
  pz_real avg_r(0), wsum(0);
  pz_real rmd(0), crmd(0);
  
  for (region::iterator i(r_int());i;++i) {
    vec_u pos         = coord(i);
    get_dens(i, rmd, crmd);
    pz_cmplx cp(pos(0), pos(1));
    cms_pos           += crmd * cp;
    cms_zpos          += crmd * pos(2);
    wsum              += crmd;
  }
  cms_pos           = global_sum(cms_pos);
  cms_zpos          = global_sum(cms_zpos);
  wsum              = global_sum(wsum);
  cms_pos           /= wsum;
  cms_zpos          /= wsum;

  for (region::iterator i(r_int());i;++i) {
    vec_u pos         = coord(i);
    get_dens(i, rmd, crmd);
    pz_cmplx cp       = pz_cmplx(pos(0), pos(1)) - cms_pos;
    pz_real d         = abs(cp);
    pz_cmplx cp2      = d > 0 ? (cp*(cp/d)) : 0.0;
    avg_r             += crmd * d;
    avg_pos           += crmd * cp2;
  }
  avg_r             = global_sum(avg_r);
  avg_pos           = global_sum(avg_pos);
  avg_r             /= wsum;
  avg_pos           /= wsum;
  pz_real new_phase = std::arg(avg_pos) / 2.0;
  pz_real new_sep   = avg_r;
 
  avg.set_phase_sep(remove_phase_jump(prev_phase*2.0, 
                        new_phase*2.0) / 2.0 , new_sep);  

  cms_x             = cms_pos.real();
  cms_y             = cms_pos.imag();
  cms_z             = cms_zpos;
}

const bns_locations& bns_tracker::track_minalp()
{
  pz_real mx(0),my(0),mz(0);
  pz_real minalp(0);
  bool first(true);
  for (region::iterator i(r_int());i;++i) {
    vec_u pos   = coord(i);
    pz_real alp = alp_pc(i);
    if (first || (alp < minalp)) {
      first = false;
      minalp  = alp;
      mx = pos(0);
      my = pos(1);
      mz = pos(2);
    }
  }
  
  pz_real mcnt   = (minalp > global_min(minalp)) ? 0.0 : 1.0;
  mx             = mx * mcnt;
  my             = my * mcnt;
  mz             = mz * mcnt;
  mcnt           = global_sum(mcnt);
  assert(mcnt>0);
  mx             = global_sum(mx);
  my             = global_sum(my);
  mz             = global_sum(mz);
  
  minalp_x       = mx / mcnt;
  minalp_y       = my / mcnt;
  minalp_z       = mz / mcnt;
  
  return *this;
}

pz_real bns_tracker::remove_phase_jump(pz_real old_phase, pz_real new_phase)
{
  while (new_phase > old_phase + M_PI) new_phase -= 2.0*M_PI;
  while (new_phase < old_phase - M_PI) new_phase += 2.0*M_PI;
  return new_phase;
}

void orb_param::set_phase_sep(pz_real phi, pz_real r)
{
  sep   = r;
  phase = phi;
  x     = r * cos(phi);
  y     = r * sin(phi);
}

const bns_locations& bns_tracker::track_stars(const double prev_phase) 
{
  if (pi_sym) {
    track_average_pisym(prev_phase);
    track_individual_pisym(avg.phase);
  }
  else {
    track_average_nosym(prev_phase);
    track_individual_nosym(avg.phase);
  }

  const double dx     = star1.x - star2.x;
  const double dy     = star1.y - star2.y;
  coord_sep_bc        = sqrt(dx * dx + dy * dy);

  const double dx_md  = star1_maxloc.x - star2_maxloc.x;
  const double dy_md  = star1_maxloc.y - star2_maxloc.y;
  coord_sep_md        = sqrt(dx_md * dx_md + dy_md * dy_md);

  return *this;
}


