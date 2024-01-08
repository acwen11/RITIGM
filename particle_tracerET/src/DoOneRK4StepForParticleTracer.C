#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/*
  Do one RK4 step for particle tracer.
*/

void Interpolate_velocities_at_particle_positions(cGH *cctkGH,int interp_num_points,double *particle_x,double *particle_y,double *particle_z,  double *particle_velx,double *particle_vely,double *particle_velz);

const double DEBUG_A = -1.31;

void debug_substitute_analytic_expressionET(int num_particles,double *x_posn,double *y_posn,double *z_posn,  double *particle_velx,double *particle_vely,double *particle_velz) {

  for(int i=0;i<num_particles;i++) {
    double xc = x_posn[i];
    double yc = y_posn[i];
    double zc = z_posn[i];
    particle_velx[i] = DEBUG_A*xc;
    particle_vely[i] = DEBUG_A*yc;
    particle_velz[i] = DEBUG_A*zc;
  }
}

void DoOneRK4StepForParticleTracerET(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(update_RK4_freq<=0 || !num_particles) return;
  if( *initial_number_of_active_refinement_levels > GetRefinementLevels(cctkGH) ) {
    // If this is ever triggered then it means
    // that the user has *removed* a refinement
    // level, which is bad. Error out.
    CCTK_VERROR("The current number of active refinement levels is %d, indicating\n"
                "that one or more refinement levels have been deactivated since\n"
                "initially %d active refinement levels were detected.\n"
                "particle_tracerET does not support this behavior.",
                GetRefinementLevels(cctkGH),
                *initial_number_of_active_refinement_levels);
  }

  if(cctk_iteration%update_RK4_freq==0 && cctk_iteration>=start_tracing_particles_iteration) {

    const double dt = *RK4_delta_time;

    /* We must solve the following 3 ODEs to get the next particle position
       dx/dt = vx
       dy/dt = vy
       dz/dt = vz

       For one variable:
       y' = f(y_n,t_n)
    */

    if(*RK4IterationCounter==1) {
      /*
	RK4 step 1:
	k_1 = dt f(y_n,t_n) = dt f(y_n) <- no explicit time dependence in f
      */
      Interpolate_velocities_at_particle_positions(cctkGH,num_particles,particle_position_x,particle_position_y,particle_position_z, particle_velx,particle_vely,particle_velz);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_substitute_analytic_expressionET(num_particles,particle_position_x,particle_position_y,particle_position_z, particle_velx,particle_vely,particle_velz);
      /***************************************************************/

      for(int i=0;i<num_particles;i++) {
	particle_position_x_k1[i] = dt*particle_velx[i];
	particle_position_y_k1[i] = dt*particle_vely[i];
	particle_position_z_k1[i] = dt*particle_velz[i];
      }
      *RK4IterationCounter=2;
      CCTK_VINFO("**** It: %d, t: %e  - RK4 Substep 1 (dt_RK4 = %e, dt_sim = %e) ****", cctk_iteration, cctk_time, dt, CCTK_DELTA_TIME);
      return;

    } else if(*RK4IterationCounter==2) {
      double *shifted_x_posn = (double *)malloc(sizeof(double)*num_particles);
      double *shifted_y_posn = (double *)malloc(sizeof(double)*num_particles);
      double *shifted_z_posn = (double *)malloc(sizeof(double)*num_particles);

      /*
	RK4 step 2:
	k_2 = dt f(y_n + 0.5*k_1,t_n + 0.5 dt) = dt f(y_n + 0.5*k_1) <- no explicit time dependence in f
      */

      for(int i=0;i<num_particles;i++) {
	shifted_x_posn[i] = particle_position_x[i] + 0.5*particle_position_x_k1[i];
	shifted_y_posn[i] = particle_position_y[i] + 0.5*particle_position_y_k1[i];
	shifted_z_posn[i] = particle_position_z[i] + 0.5*particle_position_z_k1[i];
      }

      Interpolate_velocities_at_particle_positions(cctkGH,num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_substitute_analytic_expressionET(num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);
      /***************************************************************/

      for(int i=0;i<num_particles;i++) {
	particle_position_x_k2[i] = dt*particle_velx[i];
	particle_position_y_k2[i] = dt*particle_vely[i];
	particle_position_z_k2[i] = dt*particle_velz[i];
      }
      CCTK_VINFO("**** It: %d, t: %e  - RK4 Substeps 2 & 3 (dt_RK4 = %e, dt_sim = %e) ****", cctk_iteration, cctk_time, dt, CCTK_DELTA_TIME);

      /*
	RK4 step 3:
	k_3 = dt f(y_n + 0.5*k_2,t_n + 0.5 dt) = dt f(y_n + 0.5 k_2) <- no explicit time dependence in f
      */

      for(int i=0;i<num_particles;i++) {
	shifted_x_posn[i] = particle_position_x[i] + 0.5*particle_position_x_k2[i];
	shifted_y_posn[i] = particle_position_y[i] + 0.5*particle_position_y_k2[i];
	shifted_z_posn[i] = particle_position_z[i] + 0.5*particle_position_z_k2[i];
      }

      Interpolate_velocities_at_particle_positions(cctkGH,num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_substitute_analytic_expressionET(num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);
      /***************************************************************/

      for(int i=0;i<num_particles;i++) {
	particle_position_x_k3[i] = dt*particle_velx[i];
	particle_position_y_k3[i] = dt*particle_vely[i];
	particle_position_z_k3[i] = dt*particle_velz[i];
      }

      free(shifted_x_posn);
      free(shifted_y_posn);
      free(shifted_z_posn);

      *RK4IterationCounter=4;
      return;
    } else if(*RK4IterationCounter==4) {
      double *shifted_x_posn = (double *)malloc(sizeof(double)*num_particles);
      double *shifted_y_posn = (double *)malloc(sizeof(double)*num_particles);
      double *shifted_z_posn = (double *)malloc(sizeof(double)*num_particles);

      CCTK_VINFO("**** It: %d, t: %e  - RK4 Substeps 4 & 1 (dt_RK4 = %e, dt_sim = %e) ****", cctk_iteration, cctk_time, dt, CCTK_DELTA_TIME);

      /*
	RK4 step 4:
	k_4 = dt f(y_n + k_3,t_n + dt) = dt f(y_n + k_3) <- no explicit time dependence in f
      */

      for(int i=0;i<num_particles;i++) {
	shifted_x_posn[i] = particle_position_x[i] + particle_position_x_k3[i];
	shifted_y_posn[i] = particle_position_y[i] + particle_position_y_k3[i];
	shifted_z_posn[i] = particle_position_z[i] + particle_position_z_k3[i];
      }

      Interpolate_velocities_at_particle_positions(cctkGH,num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_substitute_analytic_expressionET(num_particles,shifted_x_posn,shifted_y_posn,shifted_z_posn, particle_velx,particle_vely,particle_velz);
      /***************************************************************/

      for(int i=0;i<num_particles;i++) {
	particle_position_x_k4[i] = dt*particle_velx[i];
	particle_position_y_k4[i] = dt*particle_vely[i];
	particle_position_z_k4[i] = dt*particle_velz[i];
      }

      /******************************/
      /* DEBUG MODE: compute errors */
      if(debug) {
	double sumx,sumy,sumz,sumxe,sumye,sumze;
	sumx=sumy=sumz=sumxe=sumye=sumze=0;

	for(int i=0;i<num_particles;i++) {
	  double orig_x = particle_position_x[i];
	  double orig_y = particle_position_y[i];
	  double orig_z = particle_position_z[i];

	  //dx/dt = DEBUG_A*x -> x = x_0 e^{DEBUG_A t} = 1 + (DEBUG_A) t + (DEBUG_A t)^2/2 + (DEBUG_A t)^3/6 + (DEBUG_A t)^4/24 + O(t^5). RK4 should get all but the O(t^5) term exact.
	  double a_deltat=DEBUG_A*dt;
	  double x_exact_to_4th_order = orig_x * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);
	  double y_exact_to_4th_order = orig_y * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);
	  double z_exact_to_4th_order = orig_z * (1. + a_deltat + pow(a_deltat,2)/2. + pow(a_deltat,3)/6. + pow(a_deltat,4)/24.);

	  double x_num = (orig_x + (1.0/6.0)*(particle_position_x_k1[i] + 2.0*particle_position_x_k2[i] + 2.0*particle_position_x_k3[i] + particle_position_x_k4[i]));
	  double y_num = (orig_y + (1.0/6.0)*(particle_position_y_k1[i] + 2.0*particle_position_y_k2[i] + 2.0*particle_position_y_k3[i] + particle_position_y_k4[i]));
	  double z_num = (orig_z + (1.0/6.0)*(particle_position_z_k1[i] + 2.0*particle_position_z_k2[i] + 2.0*particle_position_z_k3[i] + particle_position_z_k4[i]));
	  sumx += x_exact_to_4th_order - x_num;
	  sumy += y_exact_to_4th_order - y_num;
	  sumz += z_exact_to_4th_order - z_num;

	  sumxe+= fabs(x_exact_to_4th_order);
	  sumye+= fabs(y_exact_to_4th_order);
	  sumze+= fabs(z_exact_to_4th_order);
	}
	CCTK_VINFO("DEBUG OUTPUT: Error in x=%e. Error in y=%e. Error in z=%e.",sumx/sumxe,sumy/sumye,sumz/sumze);
	CCTK_VINFO(" (Errors within ~8x machine precision -> particle tracer RK4 routine is coded correctly.)");
      }
      /******************************/


      for(int i=0;i<num_particles;i++) {
	/*
	  Update particle positions now. For RK4:
	  x_{n+1} = x_n + 1/6 ( k_1 + 2 k_2 + 2 k_3 + k_4 )
	*/
	particle_position_x[i] += (1.0/6.0)*(particle_position_x_k1[i] + 2.0*particle_position_x_k2[i] + 2.0*particle_position_x_k3[i] + particle_position_x_k4[i]);
	particle_position_y[i] += (1.0/6.0)*(particle_position_y_k1[i] + 2.0*particle_position_y_k2[i] + 2.0*particle_position_y_k3[i] + particle_position_y_k4[i]);
	particle_position_z[i] += (1.0/6.0)*(particle_position_z_k1[i] + 2.0*particle_position_z_k2[i] + 2.0*particle_position_z_k3[i] + particle_position_z_k4[i]);

      }

      free(shifted_x_posn);
      free(shifted_y_posn);
      free(shifted_z_posn);

      /*
	Ready for RK4 step 1 again!
	k_1 = dt f(y_n,t_n)
      */
      Interpolate_velocities_at_particle_positions(cctkGH,num_particles,particle_position_x,particle_position_y,particle_position_z, particle_velx,particle_vely,particle_velz);

      /***************************************************************/
      /* DEBUG MODE: overwrite numerical data with analytic solution */
      if(debug) debug_substitute_analytic_expressionET(num_particles,particle_position_x,particle_position_y,particle_position_z, particle_velx,particle_vely,particle_velz);
      /***************************************************************/

      for(int i=0;i<num_particles;i++) {
	particle_position_x_k1[i] = dt*particle_velx[i];
	particle_position_y_k1[i] = dt*particle_vely[i];
	particle_position_z_k1[i] = dt*particle_velz[i];
      }

      *RK4IterationCounter=2;
      return;
    }
  }
}
