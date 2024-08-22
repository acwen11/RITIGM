#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_IOMethods.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "CactusBase/IOUtil/src/ioutil_Utils.h"

void Interpolate_four_velocities_at_particle_positions(
      cGH *cctkGH,
      int interp_num_points,
      double *particle_x,
      double *particle_y,
      double *particle_z,
      double *particle_u4U0,
      double *particle_u4U1,
      double *particle_u4U2,
      double *particle_u4U3,
      double *particle_u4D0,
      double *particle_u4D1,
      double *particle_u4D2,
      double *particle_u4D3 );

void Interpolate_hydro_vars_at_particle_positions(
      cGH *cctkGH,
      int interp_num_points,
      double *particle_x,
      double *particle_y,
      double *particle_z,
      double *particle_rho,
      double *particle_T,
      double *particle_Ye,
      double *particle_W);

void Interpolate_metric_vars_at_particle_positions(
      cGH *cctkGH,
      int interp_num_points,
      double *particle_x,
      double *particle_y,
      double *particle_z,
      double *particle_psi_bssn,
      double *particle_alp);

/********************************************************************
 ********************    Macro Definitions   ************************
 ********************************************************************/
inline void CREATE_OUTDIR(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
  {
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }

    /* create the directory */
    int izzle = IOUtil_CreateDirectory (cctkGH, actual_dir, 0, 0);
    if(izzle < 0)
      CCTK_VWARN(CCTK_WARN_ALERT, "file_output_routines.C Problem creating directory '%s' for output", actual_dir);
    else if(izzle >= 0 && verbose==1)
      CCTK_VINFO("Outputting to directory '%s'", actual_dir);
  }

inline void OUTDIR_NAME(CCTK_ARGUMENTS,const char *input_dir,char *&actual_dir)
{
    DECLARE_CCTK_PARAMETERS;
    const char *_dir = input_dir;


    /* check whether "dir" was set; if not default to "IO::out_dir" */
    if (*_dir == 0)
      {
        _dir = out_dir;
      }

    /* omit the directory name if it's the current working dir */
    if (strcmp (_dir, ".") == 0)
      {
        actual_dir = strdup ("");
      }
    else
      {
        actual_dir = (char *)malloc (strlen (_dir) + 2);
        sprintf (actual_dir, "%s/", _dir);
      }
}

void particle_tracerET_file_output_routine_Startup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  char *actual_dir;
  CREATE_OUTDIR(CCTK_PASS_CTOC,outdir,actual_dir);
}

void add_header_if_file_is_empty(FILE *fp, const char *header_message) {
  fseek(fp, 0, SEEK_END);
  /* Print header if file size is zero. */
  if(!ftell(fp)) {
    char *header_buffer = (char *)malloc(sizeof(char)*500); // <- Tiny amount of memory.
    sprintf(header_buffer, "%s", header_message);
    fprintf(fp, "%s", header_buffer);
    free(header_buffer);
  }
}

void add_number_of_particles_and_quantities_if_file_is_empty(
      FILE *fp,
      const int num_particles,
      const int num_quantities ) {
  fseek(fp, 0, SEEK_END);
  /* Print header if file size is zero. */
  if(!ftell(fp)) {
    fwrite(&num_quantities, sizeof(int), 1, fp);
    fwrite(&num_particles , sizeof(int), 1, fp);
  }
}

void particle_tracerET_file_output_ascii(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char *actual_dir;
  OUTDIR_NAME(CCTK_PASS_CTOC,outdir,actual_dir);

	int np_tot = 0;
	for (int pf=0; pf<*particle_family; pf++) {
		const int np = num_particles[pf];

		char *filename = (char *)malloc(strlen(actual_dir)+strlen("particles_batch10_u4D.asc")+10);
		sprintf (filename, "%sparticles_batch%d.asc", actual_dir,pf);
		FILE *file = fopen(filename,"a+");
		if(!file) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
		add_header_if_file_is_empty(file, "# Col. 1: Time (cctk_time). Subsequent columns: x,y,z coordinate of particle i\n");

		/* Next print one line of data, corresponding to a single point in time. */
		char *buffer = (char *) malloc(sizeof(char)*5000000); // <- Allocate 5MB; sizeof(char) is 1.
		sprintf(buffer,"%e",cctk_time);
		for(int which_particle=np_tot;which_particle<np_tot+np;which_particle++) {
			sprintf(buffer, "%s %e %e %e", buffer, particle_position_x[which_particle], particle_position_y[which_particle], particle_position_z[which_particle]);
		}
		fprintf(file,"%s\n",buffer);
		fclose(file);

		if( output_four_velocity_u4U ) {
			sprintf(filename, "%sparticles_batch%d_u4U.asc", actual_dir,pf);
			if(!(file = fopen(filename, "a+"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_header_if_file_is_empty(file, "# Col. 1: Time (cctk_time). Subsequent columns: u^{t},u^{x},u^{y},u^{z} of particle i\n");
			sprintf(buffer,"%e",cctk_time);
			for(int which_particle=np_tot;which_particle<np_tot+np;which_particle++) {
				sprintf(buffer, "%s %e %e %e %e", buffer, particle_u4U0[which_particle], particle_u4U1[which_particle], particle_u4U2[which_particle], particle_u4U3[which_particle]);
			}
			fprintf(file, "%s\n", buffer);
			fclose(file);
		}
		if( output_four_velocity_u4D ) {
			sprintf(filename, "%sparticles_batch%d_u4D.asc", actual_dir,pf);
			if(!(file = fopen(filename, "a+"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_header_if_file_is_empty(file, "# Col. 1: Time (cctk_time). Subsequent columns: u_{t},u_{x},u_{y},u_{z} of particle i\n");
			sprintf(buffer,"%e",cctk_time);
			for(int which_particle=np_tot;which_particle<np_tot+np;which_particle++) {
				sprintf(buffer, "%s %e %e %e %e", buffer, particle_u4D0[which_particle], particle_u4D1[which_particle], particle_u4D2[which_particle], particle_u4D3[which_particle]);
			}
			fprintf(file, "%s\n", buffer);
			fclose(file);
		}
		if( output_hydro ) {
			sprintf(filename, "%sparticles_batch%d_hydro.asc", actual_dir,pf);
			if(!(file = fopen(filename, "a+"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_header_if_file_is_empty(file, "# Col. 1: Time (cctk_time). Subsequent columns: rho, T, Ye, W of particle i\n");
			sprintf(buffer,"%e",cctk_time);
			for(int which_particle=np_tot;which_particle<np_tot+np;which_particle++) {
				sprintf(buffer, "%s %e %e %e %e", buffer, particle_rho[which_particle], particle_T[which_particle], particle_Ye[which_particle], particle_W[which_particle]);
			}
			fprintf(file, "%s\n", buffer);
			fclose(file);
		}

		if( output_metric ) {
			sprintf(filename, "%sparticles_batch%d_metric.asc", actual_dir,pf);
			if(!(file = fopen(filename, "a+"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_header_if_file_is_empty(file, "# Col. 1: Time (cctk_time). Subsequent columns: psi_bssn, lapse of particle i\n");
			sprintf(buffer,"%e",cctk_time);
			for(int which_particle=np_tot;which_particle<np_tot+np;which_particle++) {
				sprintf(buffer, "%s %e %e", buffer, particle_psi_bssn[which_particle], particle_alp[which_particle]);
			}
			fprintf(file, "%s\n", buffer);
			fclose(file);
		}
		free(buffer);
		free(filename);
		np_tot += np;
	}
}

void particle_tracerET_file_output_binary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  char *actual_dir;
  OUTDIR_NAME(CCTK_PASS_CTOC, outdir, actual_dir);

	int np_tot = 0;
	for (int pf=0; pf<*particle_family; pf++) {
		const int np = num_particles[pf];

		char *filename = (char *)malloc(strlen(actual_dir)+strlen("particles_batch10_u4D.bin")+10);
		sprintf (filename, "%sparticles_batch%d.bin", actual_dir,pf);
		FILE *file = fopen(filename,"a+b");
		if(!file) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
		int num_quantities = 3;
		add_number_of_particles_and_quantities_if_file_is_empty(file, np, num_quantities);
		fwrite(&cctk_time         				 , sizeof(CCTK_REAL), 1            , file);
		fwrite(particle_position_x + np_tot, sizeof(CCTK_REAL), np, file);
		fwrite(particle_position_y + np_tot, sizeof(CCTK_REAL), np, file);
		fwrite(particle_position_z + np_tot, sizeof(CCTK_REAL), np, file);
		fclose(file);

		if( output_four_velocity_u4U ) {
			num_quantities = 4;
			sprintf(filename, "%sparticles_batch%d_u4U.bin",actual_dir,pf);
			if(!(file = fopen(filename, "a+b"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_number_of_particles_and_quantities_if_file_is_empty(file, np, num_quantities);
			fwrite(&cctk_time   , sizeof(CCTK_REAL), 1 , file);
			fwrite(particle_u4U0 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4U1 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4U2 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4U3 + np_tot, sizeof(CCTK_REAL), np, file);
			fclose(file);
		}
		if( output_four_velocity_u4D ) {
			sprintf(filename, "%sparticles_batch%d_u4D.bin", actual_dir,pf);
			if(!(file = fopen(filename, "a+b"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_number_of_particles_and_quantities_if_file_is_empty(file, np, num_quantities);
			fwrite(&cctk_time   , sizeof(CCTK_REAL), 1 , file);
			fwrite(particle_u4D0 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4D1 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4D2 + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_u4D3 + np_tot, sizeof(CCTK_REAL), np, file);
			fclose(file);
		}
		if( output_hydro ) {
			num_quantities = 4;
			sprintf(filename, "%sparticles_batch%d_hydro.bin", actual_dir,pf);
			if(!(file = fopen(filename, "a+b"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_number_of_particles_and_quantities_if_file_is_empty(file, np, num_quantities);
			fwrite(&cctk_time   , sizeof(CCTK_REAL), 1 , file);
			fwrite(particle_rho + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_T   + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_Ye  + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_W   + np_tot, sizeof(CCTK_REAL), np, file);
			fclose(file);
		}
		if( output_metric ) {
			num_quantities = 2;
			sprintf(filename, "%sparticles_batch%d_metric.bin", actual_dir,pf);
			if(!(file = fopen(filename, "a+b"))) CCTK_VWARN(CCTK_WARN_ABORT, "Cannot open output file '%s'", filename);
			add_number_of_particles_and_quantities_if_file_is_empty(file, np, num_quantities);
			fwrite(&cctk_time   , sizeof(CCTK_REAL), 1 , file);
			fwrite(particle_psi_bssn + np_tot, sizeof(CCTK_REAL), np, file);
			fwrite(particle_alp + np_tot, sizeof(CCTK_REAL), np, file);
			fclose(file);
		}
		free(filename);
		np_tot += np;
	}
}

void particle_tracerET_file_output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if(update_RK4_freq<=0 || !max_num_particles) return;

  if( (cctk_iteration==start_tracing_particles_iteration[0] || cctk_iteration%(*file_output_freq)==0) ) {

    if( output_four_velocity_u4D || output_four_velocity_u4U ) {
      if(verbose>1) CCTK_VINFO("**** It: %d - Interpolating four-velocities ****", cctk_iteration);
      Interpolate_four_velocities_at_particle_positions(cctkGH,
                                                        *num_active,
                                                        particle_position_x,
                                                        particle_position_y,
                                                        particle_position_z,
                                                        particle_u4U0,
                                                        particle_u4U1,
                                                        particle_u4U2,
                                                        particle_u4U3,
                                                        particle_u4D0,
                                                        particle_u4D1,
                                                        particle_u4D2,
                                                        particle_u4D3);
    }

		if ( output_hydro ) {
      if(verbose>1) CCTK_VINFO("**** It: %d - Interpolating hydro vars ****", cctk_iteration);
			Interpolate_hydro_vars_at_particle_positions(
						cctkGH,
						*num_active,
						particle_position_x,
						particle_position_y,
						particle_position_z,
						particle_rho,
						particle_T,
						particle_Ye,
						particle_W);
		}

		if ( output_metric ) {
      if(verbose>1) CCTK_VINFO("**** It: %d - Interpolating metric vars ****", cctk_iteration);
			Interpolate_metric_vars_at_particle_positions(
						cctkGH,
						*num_active,
						particle_position_x,
						particle_position_y,
						particle_position_z,
						particle_psi_bssn,
						particle_alp);
		}

    if(verbose>0) CCTK_VINFO("**** It: %d - Outputting to file ****", cctk_iteration);

    if( CCTK_MyProc(cctkGH)==0 ) {
      if( CCTK_Equals(output_file_format, "ascii") )
        particle_tracerET_file_output_ascii(CCTK_PASS_CTOC);
      else if( CCTK_Equals(output_file_format, "binary") )
        particle_tracerET_file_output_binary(CCTK_PASS_CTOC);
    }
  }
}
