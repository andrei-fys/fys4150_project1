#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

void file_writer(char* , double* , double* , int  ); 

int main (int argc, char* argv[])
{
	char *output_filename;
	//iofstream ofile;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	double *analytical_solution = new double[N];

	output_filename=argv[2];

	double h=1.0/(N+1);
	double precalc_exp=1.0-exp(-10.0);
	for (int i=1; i<N; i++) {
		grid_points[i] =  i*h;
		analytical_solution[i] = 1.0-precalc_exp*grid_points[i]-exp(-10.0*grid_points[i]);
	}
	
	file_writer(output_filename, grid_points, analytical_solution, N);

}

void file_writer(char* filename, double* g_points, double* a_solution, int n ) {
	ofstream ofile;
	ofile.open(filename);
	for (int i=1; i<n; i++) {
		ofile << g_points[i] << "," << a_solution[i] <<  endl;
	}
	ofile.close();
return;

}
