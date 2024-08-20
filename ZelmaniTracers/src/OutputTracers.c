#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define H5Acreate_vers 2
#define H5Gcreate_vers 2
#define H5Gopen_vers   2

#include <hdf5.h>
#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <util_String.h>

#define STRLEN 1024
#define MYH5CHECK(val) assert(val >= 0)

#define length(X) (sizeof(X)/sizeof(*X))

typedef struct {
  double box_size;
  int flag_cooling;
  int flag_double_precision;
  int flag_feedback;
  int flag_ic_info;
  int flag_metals;
  int flag_sfr;
  int flag_stellar_age;
  double hubble_param;
  double mass_table[6];
  int num_files_per_snapshot;
  int num_part_this_file[6];
  unsigned int num_part_total[6];
  unsigned int num_part_total_high_word[6];
  double Omega0;
  double OmegaLambda;
  double redshift;
  double time;
} Header;

typedef struct {
  char const * name;
  CCTK_REAL * data;
} Field;

static void WriteHeader(hid_t file_id, Header * header)
{
  hid_t group_id, dspace_id, attr_id;

  hsize_t const size_one = 1;
  hsize_t const size_six = 6;

  group_id = H5Gcreate(file_id, "Header", H5P_DEFAULT, H5P_DEFAULT,
          H5P_DEFAULT); MYH5CHECK(group_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "BoxSize", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->box_size);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Cooling", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_cooling);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_DoublePrecision", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_double_precision);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Feedback", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_feedback);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_IC_Info", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_ic_info);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Metals", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_metals);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_Sfr", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_sfr);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Flag_StellarAge", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->flag_stellar_age);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "HubbleParam", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->hubble_param);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "MassTable", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->mass_table);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumFilesPerSnapshot", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->num_files_per_snapshot);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_ThisFile", H5T_NATIVE_INT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_INT, &header->num_part_this_file);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_Total", H5T_NATIVE_UINT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_UINT, &header->num_part_total);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_six, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "NumPart_Total_HighWord", H5T_NATIVE_UINT,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_UINT, &header->num_part_total_high_word);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Omega0", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->Omega0);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "OmegaLambda", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->OmegaLambda);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Redshift", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->redshift);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  dspace_id = H5Screate_simple(1, &size_one, NULL); MYH5CHECK(dspace_id);
  attr_id = H5Acreate(group_id, "Time", H5T_NATIVE_DOUBLE,
          dspace_id, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(attr_id);
  H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &header->time);
  H5Aclose(attr_id);
  H5Sclose(dspace_id);

  H5Gclose(group_id);
}

static void WriteIDs(hid_t file_id, unsigned int const * IDs, int const siz)
{
  hid_t group_id, dspace_id, dset_id;
  herr_t h5err;

  hsize_t hsiz = siz;
  group_id = H5Gcreate(file_id, "PartType0", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT); MYH5CHECK(group_id);
  dspace_id = H5Screate_simple(1, &hsiz, NULL); MYH5CHECK(dspace_id);
  dset_id = H5Dcreate2(group_id, "ParticleIDs", H5T_NATIVE_UINT, dspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(dset_id);
  h5err = H5Dwrite(dset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, IDs); MYH5CHECK(h5err);

  H5Dclose(dset_id);
  H5Sclose(dspace_id);
  H5Gclose(group_id);
}

static void WriteField(hid_t file_id, float const * var, char const * name,
    int ncomp, int const siz)
{
  hid_t group_id, dspace_id, dset_id;
  herr_t h5err;

  hsize_t hsiz = siz;
  group_id = H5Gopen(file_id, "/PartType0", H5P_DEFAULT); MYH5CHECK(group_id);
  if(ncomp == 1)
  {
    dspace_id = H5Screate_simple(1, &hsiz, NULL); MYH5CHECK(dspace_id);
  }
  else if(ncomp == 3)
  {
    hsize_t dspace_siz[2] = {hsiz, 3};
    dspace_id = H5Screate_simple(2, dspace_siz, NULL); MYH5CHECK(dspace_id);
  }
  else
  {
    char msg[STRLEN];
    snprintf(msg, STRLEN, "Trying to write a field of unsupported ncomp: %d",
        ncomp);
    CCTK_ERROR(msg);
  }
  dset_id = H5Dcreate2(group_id, name, H5T_NATIVE_FLOAT, dspace_id,
      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); MYH5CHECK(dset_id);
  h5err = H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, var); MYH5CHECK(h5err);

  H5Dclose(dset_id);
  H5Sclose(dspace_id);
  H5Gclose(group_id);
}

static void WriteVisItFile(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  char *fn;
  FILE *fh;
  static int first_time = 1;
  int truncate, nprocs;

  nprocs = CCTK_nProcs(cctkGH);

  fn = NULL;
  Util_asprintf(&fn, "%s/%s.visit", out_dir, tracersfile);
  assert(fn != NULL);

  // open file and write header if needed
  truncate = first_time && IO_TruncateOutputFiles(cctkGH);
  fh = fopen(fn, truncate ? "w" : "a");
  if(fh == NULL)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Could not open file '%s' for writing: %s\n", fn,
                strerror(errno));
  }
  if(truncate)
    fprintf(fh, "!NBLOCKS %d\n", one_file_per_timestep ? 1 : nprocs);

  // one line per output file
  if(one_file_per_timestep)
  {
    fprintf(fh, "%s.%d.hdf5\n", tracersfile, cctk_iteration);
  }
  else
  {
    for(int i = 0 ; i < nprocs ; i++)
    {
      fprintf(fh, "%s.%d.%d.hdf5\n", tracersfile, cctk_iteration, i);
    }
  }

  fclose(fh);
  free(fn);

  first_time = 0;
}

typedef void (*WriteVar)(cGH const *, int, int, CCTK_REAL * const restrict *
        restrict, char const *, hid_t);

static void WriteSingleVariableChunked(cGH const * cctkGH, int ntracers_on_this_proc,
    int components, CCTK_REAL * const restrict * restrict data, char const * label,
    hid_t file_id)
{
  float * local_data = malloc(ntracers_on_this_proc * components * sizeof(float));
  for(int i = 0; i < ntracers_on_this_proc; ++i)
    for(int c = 0; c < components; ++c)
      local_data[i*components + c] = (float)data[c][i];
  WriteField(file_id, local_data, label, components, ntracers_on_this_proc);
  free(local_data);
}

static void CollectAndWriteSingleVariable(const cGH *cctkGH, int ntracers_on_this_proc,
                                          int components, CCTK_REAL * const restrict * restrict data,
                                          const char *label,
                                          hid_t file_id)
{
  DECLARE_CCTK_PARAMETERS;

  const int nProcs = CCTK_nProcs(cctkGH);
  const int myProc = CCTK_MyProc(cctkGH);

  // offsets and counts stay zero on the non-Root procs
  int *offsets = calloc(nProcs, sizeof(*offsets));
  int *counts = calloc(nProcs, sizeof(*counts));
  float *local_data, *global_data;

  int count = ntracers_on_this_proc * components;
  int count_total;

  const MPI_Comm comm = CCTK_IsFunctionAliased("GetMPICommWorld") ?
                         *(const MPI_Comm*)GetMPICommWorld(cctkGH) :
                         MPI_COMM_WORLD;

  assert(counts);
  assert(offsets);

  // local size on each process
  MPI_Gather(&count, 1, MPI_INT, counts, 1, MPI_INT, 0, comm);

  // turns size into offsets
  count_total = 0;
  for(int i = 0 ; i < nProcs ; i++)
  {
    offsets[i] = count_total;
    count_total += counts[i];
  }
  assert(count_total == ntracers*components || myProc != 0);

  // data buffers to collect data in
  local_data = malloc(count * sizeof(*local_data));
  global_data = malloc(count_total * sizeof(*global_data));
  assert(local_data);
  assert(global_data);
  for(int i = 0 ; i < ntracers_on_this_proc ; i++)
    for(int c = 0 ; c < components ; c++)
      local_data[i*components+c] = (float)data[c][i];

  // get all data
  MPI_Gatherv(local_data, count, MPI_FLOAT, global_data,
              counts, offsets, MPI_FLOAT, 0, comm);

  // write to disk
  if(myProc == 0)
    WriteField(file_id, global_data, label, components, ntracers);

  // clean up
  free(global_data);
  free(local_data);
  free(counts);
  free(offsets);
}

void ZelmaniTracers_OutputTracers(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_ARGUMENTS

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

  const int nProcs = CCTK_nProcs(cctkGH);
  const int myProc = CCTK_MyProc(cctkGH);
  int ierr;

  int ntracers_on_this_proc;
  int group;
  cGroupDynamicData groupdata;

  // collect information about tracers get data
  group = CCTK_GroupIndex("ZelmaniTracers::tracervars");
  ierr = CCTK_GroupDynamicData(cctkGH, group, &groupdata);
  if(ierr != 0)
  {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Failed to obtain information about tracervars group: %d",
                ierr);
  }
  ntracers_on_this_proc = groupdata.lsh[0];
  int const ntracers_to_write = one_file_per_timestep ? ntracers :
    ntracers_on_this_proc;

  // open file for the output
  char * fn = NULL;
  hid_t file_id = -1;
  if(one_file_per_timestep)
  {
    if(myProc == 0)
    {
      Util_asprintf(&fn, "%s/%s.%d.hdf5", out_dir, tracersfile,
          cctk_iteration);
      assert(fn);
      file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        MYH5CHECK(file_id);
    }
  }
  else
  {
    Util_asprintf(&fn, "%s/%s.%d.%d.hdf5", out_dir, tracersfile,
        cctk_iteration, myProc);
    assert(fn);
      file_id = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        MYH5CHECK(file_id);
  }

  // Construct and write header
  Header header;
  memset(&header, 0, sizeof(Header));
  header.num_files_per_snapshot = one_file_per_timestep ? 1 : nProcs;
  header.num_part_this_file[0] = ntracers_to_write;
  header.num_part_total[0] = ntracers;
  header.time = cctk_time;
  if(file_id >= 0)
    WriteHeader(file_id, &header);

  // Construct the index
  unsigned int * ID = malloc(ntracers_to_write*sizeof(*ID));
  if(one_file_per_timestep)
  {
    for(int i = 0; i < ntracers; ++i)
     ID[i] = i;
  }
  else
  {
    for(int i = 0; i < ntracers_to_write; ++i)
      ID[i] = i + groupdata.lbnd[0];
  }
  if(file_id >= 0)
    WriteIDs(file_id, ID, ntracers_to_write);

  // All of the scalar fields to output
  Field fields_base[] = {
    {"LorentzFactor",           tw_lorentz},
    {"KineticEnergyAtInfinity", teninf},
    {"Masses",                  tmass},
    {"SmoothingLength",         thsml},
    {"Density",                 trho},
    {"Temperature",             ttemp},
    {"Ye",                      tye},
    {"Lapse",                   talp},
  };
  int const n_fields_base = length(fields_base);

  Field fields_leak[] = {
    {"NeuLuminosity0",          tlum0},
    {"NeuLuminosity1",          tlum1},
    {"NeuLuminosity2",          tlum2},
    {"NeuAveEnergy0",           teave0},
    {"NeuAveEnergy1",           teave1},
    {"NeuAveEnergy2",           teave2},
  };
  int const n_fields_leak = length(fields_leak);

  Field fields_extra[] = {
    {"InternalEnergy",          teps},
    {"Pressure",                tpress},
    {"Entropy",                 tentropy},
    {"Abar",                    tabar},
    {"gxx",                     tgxx},
    {"gxy",                     tgxy},
    {"gxz",                     tgxz},
    {"gyy",                     tgyy},
    {"gyz",                     tgyz},
    {"gzz",                     tgzz},
    {"ConservedDensity",        tdens},
  };
  int const n_fields_extra = length(fields_extra);

  // Choose which function to use for writing the data
  WriteVar writer;
  if(one_file_per_timestep)
    writer = CollectAndWriteSingleVariable;
  else
    writer = WriteSingleVariableChunked;

  // Output vector vars
  CCTK_REAL * Pos[3] = {tx, ty, tz};
  writer(cctkGH, ntracers_on_this_proc, 3, Pos, "Coordinates", file_id);
  CCTK_REAL * Vel[3] = {tvx, tvy, tvz};
  writer(cctkGH, ntracers_on_this_proc, 3, Vel, "Velocities", file_id);
  if(output_extra_fields)
  {
    CCTK_REAL * Shift[3] = {tbetax, tbetay, tbetaz};
    writer(cctkGH, ntracers_on_this_proc, 3, Shift, "Shift", file_id);
  }

  // Output scalar vars
  for(int f = 0; f < n_fields_base; ++f)
  {
    writer(cctkGH, ntracers_on_this_proc, 1, &fields_base[f].data,
        fields_base[f].name, file_id);
  }
  if(tracers_with_leakage)
  {
    for(int f = 0; f < n_fields_leak; ++f)
    {
      writer(cctkGH, ntracers_on_this_proc, 1, &fields_leak[f].data,
          fields_leak[f].name, file_id);
    }
  }
  if(output_extra_fields)
  {
    for(int f = 0; f < n_fields_extra; ++f)
    {
      writer(cctkGH, ntracers_on_this_proc, 1, &fields_extra[f].data,
          fields_extra[f].name, file_id);
    }
  }

  // Write VisIt file
  if(myProc == 0 && write_VisIt_file)
    WriteVisItFile(CCTK_PASS_CTOC);

  // Cleanup
  if(file_id >= 0)
  {
    H5Fclose(file_id);
  }

  free(ID);
  free(fn);
}
