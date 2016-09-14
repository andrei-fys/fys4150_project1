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
	char *output_filename_computed;
	char *output_filename_calculated;
	char *output_filename_error;
	char *output_filename_thomas;
	char *output_filename_LU;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	double *analytical_solution = new double[N];

	output_filename_computed=argv[2];
	output_filename_calculated=argv[3];
	output_filename_error=argv[4];
	output_filename_thomas=argv[5];
	output_filename_LU=argv[6];

	double h=1.0/(N+1);
	double h_squared_100=h*h*100.0;
	double precalc_exp=1.0-exp(-10.0);
	for (int i=1; i<N; i++) {
		grid_points[i] =  i*h;
		analytical_solution[i] = 1.0-precalc_exp*grid_points[i]-exp(-10.0*grid_points[i]);
	}
	
	file_writer(output_filename_calculated, grid_points, analytical_solution, N);

	double *b = new double[N];			// main diagonal
	double *a = new double[N-1];		//lower diagonal
	double *c = new double[N-1];		//upper diagonal
	double *b_prime = new double[N];
	double *b_tilda = new double[N];
	double *b_prime_tilda = new double[N];
	double *V = new double[N];
	clock_t g_start, g_finish;
	for (int i=0;i<=N-2;i++) {
		a[i]=-1.0;
		c[i]=-1.0;
	}
	for (int i=0;i<=N-1;i++) {
		b[i]=2.0;
		b_tilda[i]=h_squared_100*exp(-10.0*grid_points[i+1]);
	}	// RHS of the equation
	b[N-1]=2.0;
	b_prime[0]=b[0];
	b_prime_tilda[0]=b_tilda[0];
	
	g_start = clock();
	// GAUSS BF START (Fiat)
	// forward subst.
	for (int j=1;j<=N-1;j++){
		b_prime[j]=b[j]-a[j-1]*c[j-1]/b_prime[j-1];
		b_prime_tilda[j]=b_tilda[j]-a[j-1]*b_prime_tilda[j-1]/b_prime[j-1];
	}
	// end of forward subst.
	// backward subst.
	V[N-1]=b_prime_tilda[N-1]/b_prime[N-1];
	for (int k=N-2;k>=0;k--) {
		V[k]=(b_prime_tilda[k+1]-c[k]*V[k+1])/b_prime[k];
	}
	// end of backward subst.
	// GAUSS BF END
	g_finish=clock();
	file_writer(output_filename_computed, grid_points, V, N);

	double *relative_error = new double[N];
	double *relative_error_log10 = new double[N];
	for (int i=0;i<=N-1;i++) {
		relative_error[i]=(abs((V[i]-analytical_solution[i])/analytical_solution[i]));
		relative_error_log10[i]=log10(relative_error[i]);
	}
	
	file_writer(output_filename_error, grid_points, relative_error, N);
	
	double thomas_a = -1.0;
	double thomas_c = -1.0;
	double thomas_b = 2.0;
	double ac = 1.0;
	double *thomas_b_prime = new double[N];
	double *thomas_b_prime_tilda = new double[N];
	double *thomas_V = new double[N];
	clock_t t_start, t_finish;
	
	thomas_b_prime[0]=thomas_b;
	thomas_b_prime_tilda[0]=b_tilda[0];
	
	t_start = clock();
	//THOMAS START (Ferrari)
	//forward subst.
	for (int j=1;j<=N-1;j++){
		thomas_b_prime[j]=thomas_b-ac/thomas_b_prime[j-1];
		thomas_b_prime_tilda[j]=b_tilda[j]-thomas_a*thomas_b_prime_tilda[j-1]/thomas_b_prime[j-1];
	}
	//end of forward subst.
	//backward subst.
	thomas_V[N-1]=thomas_b_prime_tilda[N-1]/b_prime[N-1];
	for (int k=N-2;k>=0;k--) {
		thomas_V[k]=(thomas_b_prime_tilda[k+1]-thomas_c*thomas_V[k+1])/thomas_b_prime[k];
	}
	//end of backward subst.
	//THOMAS END
	t_finish = clock();
	
	file_writer(output_filename_thomas, grid_points, thomas_V, N);

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
	
	cout << "Time of Gauss " << ((double) (g_finish-g_start)/CLOCKS_PER_SEC) << endl;
	cout << "Time of Thomas " << ((double) (t_finish-t_start)/CLOCKS_PER_SEC) << endl;
	cout << "Time of LU " << ((double) (LU_finish-LU_start)/CLOCKS_PER_SEC) << endl;
	
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
