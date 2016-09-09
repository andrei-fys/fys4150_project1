#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

int main (int argc, char* argv[])
{
	int N, i; //N - number of grid points, i - dummy index
	double h; // step lengh
	char *output_filename;
	ofstream ofile;
	double *grid_points = new double[N]; 
	double *analytical_solution = new double[N];
	double precalc_exp; //precalculated value for 1-exp(-10)

	N = atof(argv[1]);
	output_filename=argv[2];

	h=1.0/(N+1);
	precalc_exp=1.0-exp(-10.0);
	for (i=1; i<N; i++) {
		grid_points[i] =  i*h;
		analytical_solution[i] = 1.0-precalc_exp*grid_points[i]-exp(-10.0*grid_points[i]);
	}
	
	ofile.open(output_filename);
	for (i=1; i<N; i++) {
		ofile << grid_points[i] << "," << analytical_solution[i] <<  endl;
	}
	ofile.close();

}

