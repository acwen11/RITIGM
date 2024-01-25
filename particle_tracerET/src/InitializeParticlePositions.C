#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

extern void Interpolate_density_many_pts(cGH *cctkGH,int interp_num_points,double *particle_x_temp,double *particle_y_temp,double *particle_z_temp, double *particle_density_temp);

void get_random_position(double &x,double &y,double &z, const int fam) {
  DECLARE_CCTK_PARAMETERS;

  double r=1e100;
  double theta=1e100;
  while((r>seed_particles_inside_sphere__radius[fam]) || (r < seed_particles_outside_sphere__radius[fam]) 
					|| (theta < seed_particles_theta_min[fam]) || (theta > seed_particles_theta_max[fam])) {
    x = 2.0*(drand48()-0.5)*seed_particles_inside_sphere__radius[fam];
    y = 2.0*(drand48()-0.5)*seed_particles_inside_sphere__radius[fam];
    z = 2.0*(drand48()-0.5)*seed_particles_inside_sphere__radius[fam];

    r = sqrt(x*x+y*y+z*z);
    theta = atan2(z, sqrt(x*x+y*y));
  }

  x += seed_particles_inside_sphere__x_coord[fam];
  y += seed_particles_inside_sphere__y_coord[fam];
  z += seed_particles_inside_sphere__z_coord[fam];
}

/*
  Algorithm for seeding particles initially:
  1) Choose random point (x,y,z) within sphere.
  2) Probability of accepting random point = (density(x,y,z)/density_max)^central_condensation_parameter
  3) Go to (1) until all particles are seeded.
*/

void InitializeParticlePositions(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration!=start_tracing_particles_iteration[*particle_family] || !num_particles[*particle_family]) return;

	int pf = *particle_family; 
	int np = num_particles[pf];

  double *particle_x_temp  = (double *)malloc(sizeof(double)*np);
  double *particle_y_temp  = (double *)malloc(sizeof(double)*np);
  double *particle_z_temp  = (double *)malloc(sizeof(double)*np);
  double *particle_density_temp = (double *)malloc(sizeof(double)*np);

  srand48(42);

  int which_particle = *num_active;
  int total_trials = 0;
  //Technically, this algorithm is nondeterministic. However it should complete within a few iterations.
  for(int iter=0;iter<100000;iter++) {
    for(int i=0;i<np;i++) {
      // Find all particles whose positions still need to be set:
      get_random_position(particle_x_temp[i],particle_y_temp[i],particle_z_temp[i],pf);
    }
    Interpolate_density_many_pts(cctkGH,np,particle_x_temp,particle_y_temp,particle_z_temp, particle_density_temp);

    for(int i=0;i<np;i++) {
      double random_number_zero_to_one = drand48();
      if(random_number_zero_to_one < pow(particle_density_temp[i]/density_max,central_condensation_parameter)) {
        // Accept particle!
        particle_position_x[which_particle] = particle_x_temp[i];
        particle_position_y[which_particle] = particle_y_temp[i];
        particle_position_z[which_particle] = particle_z_temp[i];
        which_particle++;
      }
      total_trials++;
      if(which_particle == *num_active + np) {
        // If we've already seeded all the particles, break out of the loop!
        iter=1000000;
        i=np+100;
        CCTK_INFO("SHOULD BE ALL DONE!");
      }
    }

    if(verbose>=1 && iter!=1000000)
      CCTK_VINFO("Iteration #%d: Need to specify %d more particle location(s). Central condensation parameter = %e. Success rate = %.2e. Need ~ %d more iterations.",
                 iter,num_particles-which_particle,central_condensation_parameter,(double)which_particle/(double)total_trials, (int)((double)total_trials/(double)which_particle) - iter - 1);
    if(iter==99999)
      CCTK_WARN(CCTK_WARN_ABORT, "Hit iteration limit.");
  }
  free(particle_x_temp);
  free(particle_y_temp);
  free(particle_z_temp);
  free(particle_density_temp);

	// Update global counting vars.
	*particle_family += 1;
	*num_active += np;	
}
