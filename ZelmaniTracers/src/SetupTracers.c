#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <util_Table.h>
#include <Interpolate_density_many_pts.h>

#define SQ(X) ((X)*(X))

static inline void spher2cart(CCTK_REAL const rad, CCTK_REAL const theta,
        CCTK_REAL const phi, CCTK_REAL * x, CCTK_REAL * y, CCTK_REAL * z,
        CCTK_REAL const x0, CCTK_REAL const y0, CCTK_REAL const z0) {
    *x = x0 + rad*sin(theta)*cos(phi);
    *y = y0 + rad*sin(theta)*sin(phi);
    *z = z0 + rad*cos(theta);
}

static inline CCTK_REAL spatialdet(CCTK_REAL const hxx, CCTK_REAL const hxy,
        CCTK_REAL const hxz, CCTK_REAL const hyy, CCTK_REAL const hyz,
        CCTK_REAL const hzz) {
    return
        -SQ(hxz)*hyy + 2*hxy*hxz*hyz - hxx*SQ(hyz) - SQ(hxy)*hzz + hxx*hyy*hzz;
}

static inline CCTK_REAL mkrand(CCTK_REAL const a, CCTK_REAL const b) {
    CCTK_REAL const r = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX);
    return (b - a)*r + a;
}

void get_random_position(double *x,double *y,double *z, const int fam) {
  DECLARE_CCTK_PARAMETERS;

  const double rmin = rad_min[fam];
  const double rmax = rad_max[fam];
  // const double thmin = seed_particles_theta_min[fam] * M_PI / 180; // convert to radians
  // const double thmax = seed_particles_theta_max[fam] * M_PI / 180;
	const double thmin = 0.0;
	const double thmax = M_PI;
	// We want to uniformly sample in space:
	double r3 = mkrand(pow(rmin,3), pow(rmax,3));
	//double r3 = drand48() * (pow(rmax,3) - pow(rmin,3)) + pow(rmin,3);
	double costh = mkrand(cos(thmin), cos(thmax));
	// double costh = drand48() * (cos(thmax) - cos(thmin)) + cos(thmin);
	double phi = mkrand(0, 2 * M_PI);
	// double phi = drand48() * 2 * M_PI;

	double r = pow(r3, 1.0/3.0);
	double theta = acos(costh);

	*x = r * sin(theta) * cos(phi);
	*y = r * sin(theta) * sin(phi);
	*z = r * costh;

	if (fam == 0){
		*x += center_ns1[0];
		*y += center_ns1[1];
		*z += center_ns1[2];
	}
	else {
		*x += center_ns2[0];
		*y += center_ns2[1];
		*z += center_ns2[2];
	}
}

// This stores the conserved density on a spherical grid shared across all the
// processes
typedef struct {
    int const nrad;
    int const ntheta;
    int const nphi;
    int const npoints;

    CCTK_REAL const x0;
    CCTK_REAL const y0;
    CCTK_REAL const z0;
    CCTK_REAL const rmin;
    CCTK_REAL const rmax;

    CCTK_REAL const drad;
    CCTK_REAL const dtheta;
    CCTK_REAL const dphi;

    CCTK_REAL * x;
    CCTK_REAL * y;
    CCTK_REAL * z;

    CCTK_REAL * dens;
    CCTK_REAL * vol;
    CCTK_REAL total_mass;
    CCTK_REAL norm_const;
} DensityDistribution;

static DensityDistribution * ddist = NULL;

static inline CCTK_REAL ddist_rad(int const ir) {
    return ddist->rmin + ir*ddist->drad;
}
static inline CCTK_REAL ddist_theta(int const it) {
    return (it + 0.5)*ddist->dtheta;
}
static inline CCTK_REAL ddist_phi(int const ip) {
    return ip*ddist->dphi;
}

static inline int ddist_idx(int const ir, int const it, int const ip) {
    return ip + ddist->nphi*(it + ddist->ntheta*ir);
}
static inline void ddist_inv_idx(int const idx, int * ir, int * it, int * ip) {
    *ip = idx % ddist->nphi;
    *it = ((idx - (*ip))/ddist->nphi) % ddist->ntheta;
    *ir = (idx - (*ip) - (*it)*ddist->nphi) /(ddist->nphi*ddist->ntheta);
}

// Initialize the density distribution on the given domain
static void ddist_init(int ns_idx) {
    DECLARE_CCTK_PARAMETERS

    assert(NULL == ddist);
    ddist = malloc(sizeof(DensityDistribution)); assert(ddist);

    *(int *)&ddist->nrad    = nradialshells;
    *(int *)&ddist->ntheta  = ntheta;
    *(int *)&ddist->nphi    = nphi;
    *(int *)&ddist->npoints = nradialshells * ntheta * nphi;

    if(CCTK_Equals(TracerProblem, "BinaryNS")) {
        if(0 == ns_idx) {
            *(CCTK_REAL *)&ddist->x0 = center_ns1[0];
            *(CCTK_REAL *)&ddist->y0 = center_ns1[1];
            *(CCTK_REAL *)&ddist->z0 = center_ns1[2];

            *(CCTK_REAL *)&ddist->rmin = rad_min[0];
            *(CCTK_REAL *)&ddist->rmax = rad_max[0];
        }
        else {
            *(CCTK_REAL *)&ddist->x0 = center_ns2[0];
            *(CCTK_REAL *)&ddist->y0 = center_ns2[1];
            *(CCTK_REAL *)&ddist->z0 = center_ns2[2];

            *(CCTK_REAL *)&ddist->rmin = rad_min[1];
            *(CCTK_REAL *)&ddist->rmax = rad_max[1];
        }

        assert(!CCTK_Equals(symm, "octant"));
        *(CCTK_REAL *)&ddist->dphi = 2*M_PI/ddist->nphi;
    }
    else if(CCTK_Equals(TracerProblem, "StellarCore")) {
        *(CCTK_REAL *)&ddist->x0 = 0;
        *(CCTK_REAL *)&ddist->y0 = 0;
        *(CCTK_REAL *)&ddist->z0 = 0;

        *(CCTK_REAL *)&ddist->rmin = rad_min[0];
        *(CCTK_REAL *)&ddist->rmax = rad_max[0];

        if(CCTK_Equals(symm, "octant")) {
            *(CCTK_REAL *)&ddist->dphi = M_PI_2/ddist->nphi;
        }
        else if(CCTK_Equals(symm, "quadrant")) {
            *(CCTK_REAL *)&ddist->dphi = M_PI/ddist->nphi;
        }
        else {
            *(CCTK_REAL *)&ddist->dphi = 2*M_PI/ddist->nphi;
        }
    }
    else {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Unkown tracer problem: \"%s\"", TracerProblem);
        CCTK_ERROR(msg);
    }

    if(CCTK_Equals(symm, "full")) {
        *(CCTK_REAL *)&ddist->dtheta = M_PI/ddist->ntheta;
    }
    else {
        *(CCTK_REAL *)&ddist->dtheta = M_PI_2/ddist->ntheta;
    }

    *(CCTK_REAL *)&ddist->drad = (ddist->rmax - ddist->rmin)/(ddist->nrad - 1);

    ddist->x = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(ddist->x);
    ddist->y = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(ddist->y);
    ddist->z = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(ddist->z);

    for(int ir = 0; ir < ddist->nrad; ++ir) {
        CCTK_REAL const rad = ddist_rad(ir);
        for(int it = 0; it < ddist->ntheta; ++it) {
            CCTK_REAL const th = ddist_theta(it);
            for(int ip = 0; ip < ddist->nphi; ++ip) {
                CCTK_REAL const phi = ddist_phi(ip);
                int const idx = ddist_idx(ir, it, ip);
                spher2cart(rad, th, phi, &ddist->x[idx], &ddist->y[idx],
                        &ddist->z[idx], ddist->x0, ddist->y0, ddist->z0);
            }
        }
    }

    ddist->vol  = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(ddist->vol);
    ddist->dens = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(ddist->dens);
    ddist->total_mass = 0;
}

// Computes the density distribution on the domain
static void ddist_calc(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * rhot = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(rhot);
    CCTK_REAL * Wlrt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(Wlrt);
    CCTK_REAL * hxxt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hxxt);
    CCTK_REAL * hxyt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hxyt);
    CCTK_REAL * hxzt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hxzt);
    CCTK_REAL * hyyt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hyyt);
    CCTK_REAL * hyzt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hyzt);
    CCTK_REAL * hzzt = malloc(ddist->npoints*sizeof(CCTK_REAL)); assert(hzzt);

    int const interp_handle = CCTK_InterpHandle(interpolator);
    assert(interp_handle >= 0);
    int const options_handle = Util_TableCreateFromString(interpolator_options);
    assert(options_handle >= 0);
    int const coords_handle = CCTK_CoordSystemHandle(coordinate_system);
    assert(coords_handle >= 0);

    int npoints = ddist->npoints;

    void const * const interp_coords[3] = {ddist->x, ddist->y, ddist->z};
    CCTK_INT const input_array_indices[8] = {
        CCTK_VarIndex("HydroBase::rho"),
        CCTK_VarIndex("HydroBase::w_lorentz"),
        CCTK_VarIndex("ADMBase::gxx"),
        CCTK_VarIndex("ADMBase::gxy"),
        CCTK_VarIndex("ADMBase::gxz"),
        CCTK_VarIndex("ADMBase::gyy"),
        CCTK_VarIndex("ADMBase::gyz"),
        CCTK_VarIndex("ADMBase::gzz")
    };
    int ninputs = 8;

    CCTK_INT const output_array_types[8] = {
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL
    };

    void * const output_arrays[8] = {
      rhot, Wlrt, hxxt, hxyt, hxzt, hyyt, hyzt, hzzt
    };

    int const ierr = CCTK_InterpGridArrays(cctkGH, 3, interp_handle,
            options_handle, coords_handle, npoints, CCTK_VARIABLE_REAL,
            interp_coords, ninputs, input_array_indices, ninputs,
            output_array_types, output_arrays);
    assert(ierr==0);

    Util_TableDestroy(options_handle);

    for(int i = 0; i < ddist->npoints; ++i) {
        CCTK_REAL const det = spatialdet(hxxt[i], hxyt[i], hxzt[i],
                hyyt[i], hyzt[i], hzzt[i]);
        ddist->dens[i] = rhot[i]*Wlrt[i]*sqrt(det);
    }

    ddist->total_mass = 0;
    ddist->norm_const = 0;
    for(int ir = 0; ir < ddist->nrad; ++ir) {
        CCTK_REAL const rad = ddist_rad(ir);
        for(int it = 0; it < ddist->ntheta; ++it) {
            CCTK_REAL const th = ddist_theta(it);
            for(int ip = 0; ip < ddist->nphi; ++ip) {
                int const idx = ddist_idx(ir, it, ip);

                CCTK_REAL const det = spatialdet(hxxt[idx], hxyt[idx],
                        hxzt[idx], hyyt[idx], hyzt[idx], hzzt[idx]);
                ddist->dens[idx] = rhot[idx]*Wlrt[idx]*sqrt(det);

                ddist->vol[idx] = SQ(rad) * sin(th);
                ddist->norm_const  = fmax(ddist->norm_const,
                        ddist->dens[idx] * ddist->vol[idx]);

                ddist->vol[idx] = ddist->vol[idx] *
                    ddist->drad * ddist->dtheta * ddist->dphi;
                ddist->total_mass += ddist->vol[idx] * ddist->dens[idx];
            }
        }
    }

    free(hzzt);
    free(hyzt);
    free(hyyt);
    free(hxzt);
    free(hxyt);
    free(hxxt);
    free(Wlrt);
    free(rhot);
}

// Evaluates the density distribution at a given point
static CCTK_REAL ddist_eval(CCTK_REAL const rad, CCTK_REAL const theta,
        CCTK_REAL const phi) {
    int ir[2], it[2], ip[2];
    CCTK_REAL wr[2], wt[2], wp[2];

    /* Interpolation coefficients in the radial direction */
    assert(rad >= ddist->rmin);
    assert(rad <= ddist->rmax);
    ir[0] = (rad - ddist->rmin)/ddist->drad;
    ir[1] = ir[0] + 1;
    assert(ir[1] < ddist->nrad);
    wr[0] = (ddist_rad(ir[1]) - rad)/ddist->drad;
    wr[1] = (rad - ddist_rad(ir[0]))/ddist->drad;
    assert(wr[0] >= 0);
    assert(wr[0] <= 1);
    assert(wr[1] >= 0);
    assert(wr[1] <= 1);

    /* Interpolation coefficients in the theta direction */
    assert(theta >= ddist_theta(0) - 0.5*ddist->dtheta);
    assert(theta <= ddist_theta(ddist->ntheta-1) + 0.5*ddist->dtheta);
    if(theta <= ddist_theta(0)) {
        it[0] = 0;
        it[1] = 1;
        wt[0] = 1;
        wt[1] = 0;
    }
    else if(theta >= ddist_theta(ddist->ntheta-1)) {
        it[0] = ddist->ntheta - 2;
        it[1] = ddist->ntheta - 1;
        wt[0] = 0;
        wt[1] = 1;
    }
    else {
        it[0] = (theta - ddist_theta(0))/ddist->dtheta;
        it[1] = it[0] + 1;
        wt[0] = (ddist_theta(it[1]) - theta)/ddist->dtheta;
        wt[1] = (theta - ddist_theta(it[0]))/ddist->dtheta;
    }
    assert(wt[0] >= 0);
    assert(wt[0] <= 1);
    assert(wt[1] >= 0);
    assert(wt[1] <= 1);

    /* Interpolation coefficients in the phi direction */
    assert(phi >= ddist_phi(0));
    assert(phi <= ddist_phi(ddist->nphi));
    ip[0] = (phi - ddist_phi(0))/ddist->dphi;
    if(ip[0] > ddist->nphi - 2) {
        ip[0] = ddist->nphi - 1;
        ip[1] = 0;
    }
    else {
        ip[1] = ip[0] + 1;
    }
    wp[0] = (phi - ddist_phi(ip[0]))/ddist->nphi;
    wp[1] = 1 - wp[0];
    assert(wp[0] >= 0);
    assert(wp[0] <= 1);
    assert(wp[1] >= 0);
    assert(wp[1] <= 1);

    /* Does the interpolation */
    CCTK_REAL dens = 0;
    for(int sr = 0; sr < 2; ++sr)
    for(int st = 0; st < 2; ++st)
    for(int sp = 0; sp < 2; ++sp) {
        dens += wr[sr] * wt[st] * wp[sp] *
            ddist->dens[ddist_idx(ir[sr], it[st], ip[sp])];
    }

    return dens;
}

static void ddist_free() {
    free(ddist->dens);
    free(ddist->z);
    free(ddist->y);
    free(ddist->x);
    free(ddist);
    ddist = NULL;
}


void ZelmaniTracers_SetupTracers(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    *first_step = 1;

    int group = CCTK_GroupIndex("ZelmaniTracers::tracerevol");
    cGroupDynamicData data;
    (void)CCTK_GroupDynamicData (cctkGH, group, &data);
    int myntracers = data.lsh[0];

    bool const two_domains = CCTK_Equals(TracerProblem, "BinaryNS") &&
        !CCTK_Equals(symm, "quadrant");
    int const ndomains = two_domains ? 2 : 1;
    int const siz = two_domains ? myntracers/2 : myntracers;
    int const myoffset = two_domains ? data.lbnd[0]/2 : data.lbnd[0];
    int const toffset = data.lbnd[0];
    int myproc = CCTK_MyProc(cctkGH);

    // initialize tracers_out_of_domain
    for (int i=0; i<myntracers; i++) {
        tracers_out_of_domain[i] = 0 ;
    }

    if(CCTK_Equals(seed_method, "randomize_dens")) {
        srand(CCTK_MyProc(cctkGH));
        for(int d = 0; d < ndomains; ++d) {
            ddist_init(d);
            ddist_calc(CCTK_PASS_CTOC);

            /* We add ~10% tolerance to the density maximum to avoid issues
             * when evaluating the density at the poles or on the equatorial
             * plane (which lie out of the grid in theta) */
            CCTK_REAL const inorm = 0.9/ddist->norm_const;
            CCTK_REAL const pmass = ddist->total_mass / ntracers;

            int i = d*siz;
            while(i < (d + 1)*siz) {
                CCTK_REAL const rad = mkrand(ddist->rmin, ddist->rmax);
                CCTK_REAL const theta = mkrand(
                        ddist_theta(0) - 0.5*ddist->dtheta,
                        ddist_theta(ddist->ntheta - 1) + 0.5*ddist->dtheta);
                CCTK_REAL const phi = mkrand(ddist_phi(0),
                        ddist_phi(ddist->nphi));
                CCTK_REAL const rr = mkrand(0, 1);

                CCTK_REAL const prob = ddist_eval(rad, theta,
                        phi)*SQ(rad)*sin(theta)*inorm;
                assert(prob >= 0);
                assert(prob <= 1);
                if(rr < prob) {
                    spher2cart(rad, theta, phi, &tx[i], &ty[i], &tz[i],
                            ddist->x0, ddist->y0, ddist->z0);
										CCTK_VINFO("rad = %e; theta = %e; phi = %e; x = %e; y = %e, z = %e; x0 = %e, y0 = %e, z0 = %e",
                    		rad, theta, phi, tx[i], ty[i], tz[i],
                        ddist->x0, ddist->y0, ddist->z0);
                    tmass[i] = pmass;
                    ++i;
                }
            }

            ddist_free();
        }
    }
		else if (CCTK_Equals(seed_method, "randomize_vol")){

				srand(CCTK_MyProc(cctkGH));
        for(int d = 0; d < ndomains; ++d) {
						double *particle_x_temp  = (double *)malloc(sizeof(double)*siz);
						double *particle_y_temp  = (double *)malloc(sizeof(double)*siz);
						double *particle_z_temp  = (double *)malloc(sizeof(double)*siz);
						double *particle_density_temp = (double *)malloc(sizeof(double)*siz);
						double *particle_volform_temp = (double *)malloc(sizeof(double)*siz);

						int which_particle = d * siz;
						int total_trials = 0;
						// Technically, this algorithm is nondeterministic. However it should complete within a few iterations.
						for(int iter=0;iter<100000;iter++) {
								for(int i=0;i<siz;i++) {
										// Find all particles whose positions still need to be set:
										get_random_position(&particle_x_temp[i],&particle_y_temp[i],&particle_z_temp[i],d);
								}
								Interpolate_density_many_pts(cctkGH,siz,particle_x_temp,particle_y_temp,particle_z_temp, particle_density_temp, particle_volform_temp);	
								for(int i=0;i<siz;i++) {
										if(particle_density_temp[i] > init_dens_min){
												// Accept particle!
												tx[which_particle] = particle_x_temp[i];
												ty[which_particle] = particle_y_temp[i];
												tz[which_particle] = particle_z_temp[i];
												tmass[which_particle] = particle_density_temp[i] * particle_volform_temp[i];
												which_particle++;
										}
										total_trials++;
										if(which_particle == (d + 1) * siz) {
												// If we've already seeded all the particles, break out of the loop!
												iter=1000000;
												i=siz+100;
												CCTK_INFO("SHOULD BE ALL DONE!");
										}
								}
								if(iter!=1000000)
									CCTK_VINFO("Iteration #%d: Need to specify %d more particle location(s). Success rate = %.2e. Need ~ %d more iterations.",
														 iter,siz-which_particle,(double)which_particle/(double)total_trials, (int)((double)total_trials/(double)which_particle) - iter - 1);
								if(iter==99999)
									CCTK_WARN(CCTK_WARN_ABORT, "Hit iteration limit.");
						}
						free(particle_x_temp);
						free(particle_y_temp);
						free(particle_z_temp);
						free(particle_density_temp);
				}
		}
    else if (CCTK_Equals(seed_method, "on_grid")){
        for(int d = 0; d < ndomains; ++d) {
            ddist_init(d);
            ddist_calc(CCTK_PASS_CTOC);

            for(int i = 0; i < siz; ++i) {
                tx[i + d*siz]    = ddist->x[i + myoffset];
                ty[i + d*siz]    = ddist->y[i + myoffset];
                tz[i + d*siz]    = ddist->z[i + myoffset];
                tmass[i + d*siz] = ddist->dens[i + myoffset]*ddist->vol[i + myoffset];
            }

            ddist_free();
        }
    }
    else {
        char msg[BUFSIZ];
        snprintf(msg, BUFSIZ, "Unkown seeding method: \"%s\"", seed_method);
        CCTK_ERROR(msg);
    }
}
