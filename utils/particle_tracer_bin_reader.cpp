#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>

using namespace std;

void read_coords(int num_batches){

	cout << "Reading coords..." << endl;
	for(int nb=0; nb<num_batches; nb++){
		cout << "Batch " << nb << "." << endl;
    char *in_filename = (char *)malloc(strlen("particles_batch10_u4D.bin")+10);
    sprintf (in_filename, "particles_batch%d.bin", nb);
		FILE *in_file = fopen(in_filename, "rb");
		assert(in_file);

    char *out_filename = (char *)malloc(strlen("particles_batch10_u4D.asc")+10);
    sprintf (out_filename, "particles_batch%d.asc", nb);
    FILE *out_file = fopen(out_filename,"a+");
		assert(out_file);

		
		// Read in binary
		int np;
		int nq;
		fread(&nq, sizeof(int), 1, in_file);
		fread(&np, sizeof(int), 1, in_file);
		cout << "np = " << np << endl;
		cout << "nq = " << nq << endl;

		char* buffer = (char*) malloc (sizeof(char)*5000000); // TODO: Actually set this size intelligently. But this number works for the thorn...
		double time;
		while(fread(&time, sizeof(double), 1, in_file) != 0) {
			//cout << "Reading file..."  << endl;
			double* xcoords = new double[np];
			double* ycoords = new double[np];
			double* zcoords = new double[np];

			fread(xcoords, sizeof(double), np, in_file);
			fread(ycoords, sizeof(double), np, in_file);
			fread(zcoords, sizeof(double), np, in_file);

			// Write to ascii
			//cout << "Writing file..."  << endl;
			sprintf(buffer, "%e", time);	
			for(int ip=0; ip<np; ip++){
				sprintf(buffer, "%s %e %e %e", buffer, xcoords[ip], ycoords[ip], zcoords[ip]);
			}
    	fprintf(out_file, "%s\n", buffer);

			free(xcoords);
			free(ycoords);
			free(zcoords);
		}

    fclose(out_file);
    fclose(in_file);
		free(buffer);
		free(in_filename);
		free(out_filename);
	}
}

void read_u4u(int num_batches){

	cout << "Reading u4U..." << endl;
	for(int nb=0; nb<num_batches; nb++){
		cout << "Batch " << nb << "." << endl;
    char *in_filename = (char *)malloc(strlen("particles_batch10_u4D.bin")+10);
    sprintf (in_filename, "particles_batch%d_u4U.bin", nb);
		FILE *in_file = fopen(in_filename, "rb");
		assert(in_file);

    char *out_filename = (char *)malloc(strlen("particles_batch10_u4D.asc")+10);
    sprintf (out_filename, "particles_batch%d_u4U.asc", nb);
    FILE *out_file = fopen(out_filename,"a+");
		assert(out_file);

		
		// Read in binary
		int np;
		int nq;
		fread(&nq, sizeof(int), 1, in_file);
		fread(&np, sizeof(int), 1, in_file);
		cout << "np = " << np << endl;
		cout << "nq = " << nq << endl;

		char* buffer = (char*) malloc (sizeof(char)*5000000); // TODO: Actually set this size intelligently. But this number works for the thorn...
		double time;
		while(fread(&time, sizeof(double), 1, in_file) != 0) {
			//cout << "Reading file..."  << endl;
			double* u4U0 = new double[np];
			double* u4U1 = new double[np];
			double* u4U2 = new double[np];
			double* u4U3 = new double[np];

			fread(u4U0, sizeof(double), np, in_file);
			fread(u4U1, sizeof(double), np, in_file);
			fread(u4U2, sizeof(double), np, in_file);
			fread(u4U3, sizeof(double), np, in_file);

			// Write to ascii
			//cout << "Writing file..."  << endl;
			sprintf(buffer, "%e", time);	
			for(int ip=0; ip<np; ip++){
				sprintf(buffer, "%s %e %e %e %e", buffer, u4U0[ip], u4U1[ip], u4U2[ip], u4U3[ip]);
			}
    	fprintf(out_file, "%s\n", buffer);

			free(u4U0);
			free(u4U1);
			free(u4U2);
			free(u4U3);
		}

    fclose(out_file);
    fclose(in_file);
		free(buffer);
		free(in_filename);
		free(out_filename);
	}
}

void read_u4d(int num_batches){

	cout << "Reading u4D..." << endl;
	for(int nb=0; nb<num_batches; nb++){
		cout << "Batch " << nb << "." << endl;
    char *in_filename = (char *)malloc(strlen("particles_batch10_u4D.bin")+10);
    sprintf (in_filename, "particles_batch%d_u4D.bin", nb);
		FILE *in_file = fopen(in_filename, "rb");
		assert(in_file);

    char *out_filename = (char *)malloc(strlen("particles_batch10_u4D.asc")+10);
    sprintf (out_filename, "particles_batch%d_u4D.asc", nb);
    FILE *out_file = fopen(out_filename,"a+");
		assert(out_file);

		
		// Read in binary
		int np;
		int nq;
		fread(&nq, sizeof(int), 1, in_file);
		fread(&np, sizeof(int), 1, in_file);
		cout << "np = " << np << endl;
		cout << "nq = " << nq << endl;

		char* buffer = (char*) malloc (sizeof(char)*5000000); // TODO: Actually set this size intelligently. But this number works for the thorn...
		double time;
		while(fread(&time, sizeof(double), 1, in_file) != 0) {
			//cout << "Reading file..."  << endl;
			double* u4D0 = new double[np];
			double* u4D1 = new double[np];
			double* u4D2 = new double[np];
			double* u4D3 = new double[np];

			fread(u4D0, sizeof(double), np, in_file);
			fread(u4D1, sizeof(double), np, in_file);
			fread(u4D2, sizeof(double), np, in_file);
			fread(u4D3, sizeof(double), np, in_file);

			// Write to ascii
			//cout << "Writing file..."  << endl;
			sprintf(buffer, "%e", time);	
			for(int ip=0; ip<np; ip++){
				sprintf(buffer, "%s %e %e %e %e", buffer, u4D0[ip], u4D1[ip], u4D2[ip], u4D3[ip]);
			}
    	fprintf(out_file, "%s\n", buffer);

			free(u4D0);
			free(u4D1);
			free(u4D2);
			free(u4D3);
		}

    fclose(out_file);
    fclose(in_file);
		free(buffer);
		free(in_filename);
		free(out_filename);
	}
}

void read_hydro(int num_batches){

	cout << "Reading hydro..." << endl;
	for(int nb=0; nb<num_batches; nb++){
		cout << "Batch " << nb << "." << endl;
    char *in_filename = (char *)malloc(strlen("particles_batch10_hydro.bin")+10);
    sprintf (in_filename, "particles_batch%d_hydro.bin", nb);
		FILE *in_file = fopen(in_filename, "rb");
		assert(in_file);

    char *out_filename = (char *)malloc(strlen("particles_batch10_hydro.asc")+10);
    sprintf (out_filename, "particles_batch%d_hydro.asc", nb);
    FILE *out_file = fopen(out_filename,"a+");
		assert(out_file);

		
		// Read in binary
		int np;
		int nq;
		fread(&nq, sizeof(int), 1, in_file);
		fread(&np, sizeof(int), 1, in_file);
		cout << "np = " << np << endl;
		cout << "nq = " << nq << endl;

		char* buffer = (char*) malloc (sizeof(char)*5000000); // TODO: Actually set this size intelligently. But this number works for the thorn...
		double time;
		while(fread(&time, sizeof(double), 1, in_file) != 0) {
			//cout << "Reading file..."  << endl;
			double* pt_rho = new double[np];
			double* pt_T = new double[np];
			double* pt_Ye = new double[np];
			double* pt_W = new double[np];

			fread(pt_rho, sizeof(double), np, in_file);
			fread(pt_T, sizeof(double), np, in_file);
			fread(pt_Ye, sizeof(double), np, in_file);
			fread(pt_W, sizeof(double), np, in_file);

			// Write to ascii
			//cout << "Writing file..."  << endl;
			sprintf(buffer, "%e", time);	
			for(int ip=0; ip<np; ip++){
				sprintf(buffer, "%s %e %e %e %e", buffer, pt_rho[ip], pt_T[ip], pt_Ye[ip], pt_W[ip]);
			}
    	fprintf(out_file, "%s\n", buffer);

			free(pt_rho);
			free(pt_T);
			free(pt_Ye);
			free(pt_W);
		}

    fclose(out_file);
    fclose(in_file);
		free(buffer);
		free(in_filename);
		free(out_filename);
	}
}

void read_metric(int num_batches){

	cout << "Reading metric..." << endl;
	for(int nb=0; nb<num_batches; nb++){
		cout << "Batch " << nb << "." << endl;
    char *in_filename = (char *)malloc(strlen("particles_batch10_metric.bin")+10);
    sprintf (in_filename, "particles_batch%d_metric.bin", nb);
		FILE *in_file = fopen(in_filename, "rb");
		assert(in_file);

    char *out_filename = (char *)malloc(strlen("particles_batch10_metric.asc")+10);
    sprintf (out_filename, "particles_batch%d_metric.asc", nb);
    FILE *out_file = fopen(out_filename,"a+");
		assert(out_file);

		
		// Read in binary
		int np;
		int nq;
		fread(&nq, sizeof(int), 1, in_file);
		fread(&np, sizeof(int), 1, in_file);
		cout << "np = " << np << endl;
		cout << "nq = " << nq << endl;

		char* buffer = (char*) malloc (sizeof(char)*5000000); // TODO: Actually set this size intelligently. But this number works for the thorn...
		double time;
		while(fread(&time, sizeof(double), 1, in_file) != 0) {
			//cout << "Reading file..."  << endl;
			double* psi = new double[np];
			double* alp = new double[np];

			fread(psi, sizeof(double), np, in_file);
			fread(alp, sizeof(double), np, in_file);

			// Write to ascii
			//cout << "Writing file..."  << endl;
			sprintf(buffer, "%e", time);	
			for(int ip=0; ip<np; ip++){
				sprintf(buffer, "%s %e %e", buffer, psi[ip], alp[ip]);
			}
    	fprintf(out_file, "%s\n", buffer);

			free(psi);
			free(alp);
		}

    fclose(out_file);
    fclose(in_file);
		free(buffer);
		free(in_filename);
		free(out_filename);
	}
}

int main(){
	using namespace std;
	int num_batches;
	bool out_u4u;
	bool out_u4d;
	bool out_hydro;
	bool out_metric;

	// Read in a file formatted such that:
	// # batches (starting from 0)
	// 1/0 (output u4u, yes or no)
	// 1/0 (output u4d, yes or no)
	// 1/0 (output hydro, yes or no)
	// 1/0 (output metric vars, yes or no)
	
	cout << "Reading config file." << endl;
	
	ifstream pt_config;
	pt_config.open("particles_config.txt");
	assert(pt_config.is_open());

	pt_config >> num_batches;
	pt_config >> out_u4u;
	pt_config >> out_u4d;
	pt_config >> out_hydro;
	pt_config >> out_metric;

	num_batches = num_batches + 1;

	// Read in coordinates	
	read_coords(num_batches);

	// Read u4u
	if (out_u4u)
		read_u4u(num_batches);	

	// Read u4d
	if (out_u4d)
		read_u4d(num_batches);	

	// Read hydro 
	if (out_hydro)
		read_hydro(num_batches);	

	// Read metric
	if (out_metric)
		read_metric(num_batches);	

	return 0;
}
