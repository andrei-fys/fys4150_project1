#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;

void file_writer(char* , double* , double* , int  ); 

int main (int argc, char* argv[])
{
	char *output_filename_computed;
	char *output_filename_calculated;
	//iofstream ofile;
	int N = atof(argv[1]);
	double *grid_points = new double[N]; 
	double *analytical_solution = new double[N];

	output_filename_computed=argv[2];
	output_filename_calculated=argv[3];

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

	for (int i=0;i<=N-2;i++) {
		a[i]=-1.0;
		c[i]=-1.0;
	}

	for (int i=0;i<=N-1;i++) {
		b[i]=2.0;
		b_tilda[i]=h_squared_100*exp(-10.0*grid_points[i+1]);
	}
	b[N-1]=2.0;
	b_prime[0]=b[0];
	b_prime_tilda[0]=b_tilda[0];
	for (int j=1;j<=N-1;j++){
		b_prime[j]=b[j]-a[j-1]*c[j-1]/b_prime[j-1];
		b_prime_tilda[j]=b_tilda[j]-a[j-1]*b_prime_tilda[j-1]/b_prime[j-1];
		//cout << "b_prime_tilda " << b_prime_tilda[j] ;
		//cout << "b_prime " << b_prime[j] << endl;

	}
	V[N-1]=b_prime_tilda[N-1]/b_prime[N-1];
	//cout << "#########" << V[N-1] << endl;
	for (int k=N-2;k>=0;k--) {
		V[k]=(b_prime_tilda[k+1]-c[k]*V[k+1])/b_prime[k];
		//cout << "k  " << k << endl;
		//cout << "V  " << V[k] << endl;
	}
	file_writer(output_filename_computed, grid_points, V, N);
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
