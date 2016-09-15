#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "lib.h"
#include "time.h"

using namespace std;

void file_writer(char* , double* , double* , int  ); 

int main (int argc, char* argv[])
{
	char *output_filename_LU;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 

	output_filename_LU=argv[2];

	double h=1.0/(N+1);
	double h_squared_100=h*h*100.0;
	
	for (int i=0; i<N; i++) {
		grid_points[i] = (i+1)*h;
	}
	double *b_tilda = new double[N];
	
	for (int i=0;i<N;i++) {                                                     
		b_tilda[i]=h_squared_100*exp(-10.0*grid_points[i]);                     
	} // RHS of the equation
	
	if (N <= 1000) {
		//LU decomposition
		//RHS is precalculated before and stored in b_tilda array
		double ** AA = new double*[N];
		clock_t LU_start, LU_finish;
		for (int i = 0; i < N; i++){
			AA[i] = new double[N];
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++){
			//AA[i][j] = 0.0;
			if (i == j) AA[i][j] = 2.0;
			if (abs(i - j) == 1) AA[i][j] = -1.0;
			}
		}
		int *indx = new int[N];
		double d;
		//ludcmp(double **a, int n, int *indx, double *d)
		
		//LU decomposition start
		LU_start = clock();
		ludcmp(AA, N, indx, &d);
		lubksb(AA, N, indx, b_tilda);
		LU_finish = clock();
		//LU decomposition end
		
		file_writer(output_filename_LU, grid_points, b_tilda, N);              
		cout << "Time of LU " << ((double) (LU_finish-LU_start)/CLOCKS_PER_SEC) << endl;
	}
	


	delete [] b_tilda;
//	for (int i = 0; i < N; i++){
//		delete[] AA[i];
//		delete[] AA;
//	}
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
