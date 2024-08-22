#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define length(X) (sizeof(X)/sizeof(*X))
#define MYH5CHECK(val) assert(val >= 0)

// Stores all the tracers info
typedef struct {
    CCTK_REAL * x_coords;
    CCTK_REAL * y_coords;
    CCTK_REAL * z_coords;
    CCTK_REAL * mb;
    const int npoints;
} myTracers;

static myTracers * trcrs = NULL;

static void trcrs_free() {
    free(trcrs->x_coords);
    free(trcrs->y_coords);
    free(trcrs->z_coords);
    free(trcrs->mb);
    free(trcrs);
    trcrs = NULL;
}


int Read_Positions(const char *filename, bool manual_tracers_mass) {

    double * coords = NULL;
    double * mass_b = NULL;
    int rank, rank_m;

    int const npoints;

    /* Decl. HDF vars */
    hid_t fid, dataid, fapl, dataid_m;
    hsize_t *dim = NULL, *dim_m=NULL, size, size_m ;
    hid_t datatype, dataspace, dataspace_m;
    herr_t hErrVal;

    printf("Reading from file: %s\n",filename);

    /* Load the library -- not required for most platforms. */
    hErrVal = H5open(); MYH5CHECK(hErrVal);

		/*
    * Open file, thereby enforcing proper file close
    * semantics
    */
    fapl = H5Pcreate(H5P_FILE_ACCESS); MYH5CHECK(fapl);
    H5Pset_fclose_degree(fapl,H5F_CLOSE_STRONG);
    fid = H5Fopen(filename, H5F_ACC_RDONLY,fapl); MYH5CHECK(fid);
    H5Pclose(fapl);

    /*
    * open and read the coords dataset with the initial positions
    */
    dataid = H5Dopen(fid,"/initial_pos/coords",H5P_DEFAULT); MYH5CHECK(dataid);
    dataspace = H5Dget_space(dataid); MYH5CHECK(dataspace);
    rank = H5Sget_simple_extent_ndims(dataspace); MYH5CHECK(rank);
    dim = malloc(rank*sizeof(hsize_t));
		H5Sget_simple_extent_dims(dataspace, dim, NULL); 
    size = dim[0]*dim[1];
    coords = malloc(size*sizeof(double));


    /*
    * Coordinates are stored as double
    */
    H5Dread(dataid, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords);
    H5Dclose(dataid);
    H5Sclose(dataspace);

    /*
    * open and read the mass_b dataset for tracers
    */
    if (manual_tracers_mass) {
        dataid_m = H5Dopen(fid,"/initial_pos/mass_b",H5P_DEFAULT); MYH5CHECK(dataid_m);
        dataspace_m = H5Dget_space(dataid_m); MYH5CHECK(dataspace_m);
        rank_m = H5Sget_simple_extent_ndims(dataspace_m); MYH5CHECK(rank_m);
        dim_m = malloc(rank_m*sizeof(hsize_t));
        H5Sget_simple_extent_dims(dataspace_m, dim_m, NULL);
        size_m = dim_m[0];
        mass_b = malloc(size_m*sizeof(double));

        H5Dread(dataid_m, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mass_b);
        H5Dclose(dataid_m);
        H5Sclose(dataspace_m);
    }
    H5Fclose(fid);

		printf("Successfully read data from file: %s\n",filename);

    /*
    * Allocate memory for the coordinate vectors
    */
    trcrs->x_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->x_coords);
    trcrs->y_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->y_coords);
    trcrs->z_coords = malloc(dim[0]*sizeof(CCTK_REAL)); assert(trcrs->z_coords);
    *(int *)&trcrs->npoints  = dim[0];


    int index = 0;

		printf("Begin setting up tracers\n");
    for(int i = 0; i < size; i = i+3)
    {
        //printf("Index %i \n", index);
        trcrs->x_coords[index] = coords[i];
        trcrs->y_coords[index] = coords[i+1];
        trcrs->z_coords[index] = coords[i+2];
        /*
         * DEBUG
         *
        if (index < 50){
            printf("X %8.2f \n",trcrs->x_coords[index]);
            printf("Y %8.2f \n",trcrs->y_coords[index]);
            printf("Z %8.2f \n",trcrs->z_coords[index]);
        }
        */
        index++;
    }

    if (manual_tracers_mass) {
        trcrs->mb = malloc(dim_m[0]*sizeof(CCTK_REAL)); assert(trcrs->mb);
        for (int j = 0; j < size_m; j++) {
            trcrs->mb[j] = mass_b[j];
        }
        free(mass_b) ;
        free(dim_m) ;
        mass_b = NULL ;
        dim_m = NULL ;
    }
    
    // Free coords and dim
    free(coords);
    free(dim);
    coords = NULL;
    dim = NULL;

    //printf("X %8.2f Y %8.2f Z %8.2f\n",trcrs->x_coords[0],trcrs->y_coords[0],trcrs->z_coords[0]);
    printf("ZelmaniTracers:: Number of tracers %i \n", index);
    printf("ZelmaniTracers:: Done reading the file.\n");

    return 0;

}

void Manual_SetupTracers(CCTK_ARGUMENTS) {

    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int  Err = 0;

    printf("ZelmaniTracers:: Enter Manual_SetupTracers\n");
    assert(NULL == trcrs);
    trcrs = malloc(sizeof(myTracers)); assert(trcrs);

    Err =  Read_Positions(tracers_pos_file, manual_tracers_mass); assert(Err==0);
    
    int group = CCTK_GroupIndex("ZelmaniTracers::tracerevol");
    cGroupDynamicData data;
    (void)CCTK_GroupDynamicData (cctkGH, group, &data);
    int myntracers = data.lsh[0];
		/* 
  		assumes that tracers are distributed over many procs.
  		This is not true here. so myoffset=0. 
  	*/

    int myoffset = data.lbnd[0];
    const int siz = myntracers;
 
    int npoints = trcrs->npoints;

    /* For now we assume that the tracers are "massless".
     * We assign to them only an initial position.
     */

    for(int i = 0 ; i < siz; ++i) {
                tx[i]    = trcrs->x_coords[i+myoffset];
                ty[i]    = trcrs->y_coords[i+myoffset];
                tz[i]    = trcrs->z_coords[i+myoffset];
                if (manual_tracers_mass) {
                    tmass[i] = trcrs->mb[i+myoffset];
                }
                else {
                    tmass[i] = 0.; //Tracers do not have a mass :)
                }
		// initialize tracers_out_of_domain
		tracers_out_of_domain[i] = 0 ;
    };

    trcrs_free();

    printf("ZelmaniTracers:: Exit Manual_SetupTracers\n");

}
