#ifndef BNSTRACKER_H
#define BNSTRACKER_H
#include "pizzacactus.h"

namespace Pizza {
namespace BNSTrackerGen {

enum merger_state_t {UNKNOWN=-1, INSPIRAL=0, MERGING=1, COLLAPSING=2};

struct orb_param {
  pz_real phase;
  pz_real sep;
  pz_real x;
  pz_real y;
  void set_phase_sep(pz_real phi, pz_real r);
};

struct bns_locations {
  orb_param star1;
  orb_param star1_maxloc;
  orb_param star2;
  orb_param star2_maxloc;
  orb_param avg;
  pz_real coord_sep_bc;
  pz_real coord_sep_md;
  pz_real cms_x;
  pz_real cms_y;
  pz_real cms_z;
  pz_real minalp_x;
  pz_real minalp_y;
  pz_real minalp_z;
};

class bns_tracker : cactus_grid, bns_locations {
  scalar_pc rmd_pc;
  scalar_pc alp_pc;
  scalar_pc lfac_pc;
  mats_l_pc glo_pc;

  int use_reflevel;
  bool pi_sym;
  int dist_num_points;

  static pz_real remove_phase_jump(pz_real old_phase, pz_real new_phase);
  pz_real sym_weight(vec_u& pos) const;
  void get_dens(const gpos& i, pz_real& rmd, pz_real& crmd);
  void track_individual_nosym(pz_real phase_est);
  void track_individual_pisym(pz_real phase_est);
  void track_average_nosym(const double prev_phase);
  void track_average_pisym(const double prev_phase);

  public:
  bns_tracker(int use_reflevel_, bool pi_sym_, int dist_num_points_);
  void prep(const cGH *cgh);
  const bns_locations& track_stars(const double prev_phase);
  const bns_locations& track_minalp();

};

}
}


#endif

