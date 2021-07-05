/*  This is an implementation of the Levenberg-Marquardt-Algorithm as given in

"W. H. Press, S. A. Teukolksy et al.: Numerical Recipes in C: The Art of
"Scientific Computing, 2nd ed., CUP, 1992

This program takes as input a file of the format (x_i, y_i, sigma_i) (sigma's
are optional) and uses Levenberg-Marquardt-Algorithm to fit power law with
corrections (optional) of the form

y(x) = A * x **(-B) * (1 + C * x ** (-D))

and returns parameters A,B,C,D with their errors and goodness-of-fit estimates.
For description of algorithm see Numerical Recipes as above.

It is highly non-optimised yet and primarily written for easy readibility and debugging.

For questions: benjamin (b.walter16@imperial.ac.uk) */

#include <stdio.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <math.h>
#include <string.h>

#include "nrutil_double.c"
#include "nrutil_double.h"

long double swap;
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define LAMDA_STRETCH 10 // Levenberg-Marquardt stretching factor
#define COLS 11
#define ITMAX 100 // For gammaq etc.
#define EPS 3.0e-7
#define FPMIN 1.0e-30


void fpowerlaw(long double, long double [], long double *, long double [], int);
void fpowerlaw_correction_fixed_exponent(long double x, long double a[], long double *y, long double dyda[], int na);
void fpowerlaw_correction_variable_exponent(long double x, long double a[], long double *y, long double dyda[], int na);
void mrqmin(long double x[], long double y[], long double sig[], int ndata, long double a[], int ia[], int ma, long double **covar, long double **alpha, long double *chisq, void (*funcs)(long double, long double[], long double *, long double [], int), long double *alamda);
int count_lines(char *filename, int xcol, int ycol, int scol, long double min, long double max);
void read_data(char *filename, int xcol, int ycol, int scol, long double min, long double max, long double xdata[], long double ydata[], long double sigma[]);
void printhelp(void);
void gcf(long double *gammcf, long double a, long double x, long double *gln);
void gser(long double *gamser, long double a, long double x, long double *gln);
long double gammq(long double a, long double x);
long double gammln(long double xx);

int verbose = 1;
int columns=1;

int main(int argc, char *argv[]) {  
	printf("Test 1\n");
	if(verbose) printf("# PowerFit \n# Levenberg-Marquardt-Method \n# (From Numerical Recipes in C: The Art of Scientific Computing)\n# Fits Power Law y(x) = A * x ** -B to data of format (x_i, y_i, sigma_i). \n# Returns parameters, deviations, and goodness-of-fit.\n# Syntax: ./powerfit [Options] input_data.dat \n");     
	
	// Variable Declaration
	int i, j;	
	char *filename;
	long double data_min = -1;
	long double data_max = -1;
	int x_col = -1, y_col = -1, sigma_col = -1;
	long double a_init_guess = 1.0, b_init_guess = 2.0;
	long double *xdata, *ydata, *sigma;
	int data_count;
	// Getopt routine
	opterr = 0;
	int c = 0;
	while( (c = getopt (argc, argv, "m:M:x:y:s:A:B:h:") ) != -1)
		switch(c)
			{
				case 'm':
					data_min = atof(optarg);
					break;
				case 'M':
					data_max = atof(optarg);
					break;
				case 'x':
					x_col = atoi(optarg);
					break;
				case 'y':
					y_col = atoi(optarg);
					break;
				case 's':
					sigma_col = atoi(optarg);
					break;
				case 'A':
					a_init_guess = atof(optarg);
					break;
				case 'B':
					b_init_guess = atof(optarg);
					break;
				case 'h':
					printhelp();
					break;
			default:
				exit(EXIT_FAILURE);
			}

	if ( optind == argc - 1)
	{
		filename = argv[optind];
	}
{
	// Setting default values
	if(x_col == -1) x_col = 1;
	if(y_col == -1) y_col = 2;
	// Read file into program
	data_count = count_lines(filename, x_col, y_col, sigma_col, data_min, data_max);
	xdata = vector(1, data_count);
	ydata = vector(1, data_count);
	sigma = vector(1, data_count);
	if(sigma_col == -1){
		for(i = 1; i <= data_count; i++) {sigma[i] = 1.0;}
		printf("No errorbars have been given. Calculation is proceeded assuming sigma = 1.0 for all i.\n");
	}
	read_data(filename, x_col, y_col, sigma_col, data_min, data_max, xdata, ydata, sigma);

	for(i = 1; i <= data_count; i++){
		printf("x[i] = %Lf, %Lf, %Lf\n",xdata[i],ydata[i],sigma[i]);
	}
	/* !!! Check whether m < x[i] <= M!!! */
	
	// Generate noisy data
	/*long double *xdata, *ydata, *sigma;
	xdata = vector(1,ndata);
	ydata = vector(1, ndata);
	sigma = vector(1, ndata);
		for(i = 1; i <= ndata; i++){
			xdata[i] = 0.1*(long double) i;
			ydata[i] = 100.0 * pow( xdata[i], -2.0) + 5.0*(1 - 2 * ((long double) rand()) / (long double) RAND_MAX);
			//ydata[i] = 1.0 + 2.0*xdata[i] + 3.0*xdata[i]*xdata[i];
			sigma[i] = 1.0;
			//printf("#DATA %f\t%f\t%f\n", xdata[i], ydata[i], sigma[i]);
		}
	*/
	// Initialise parameters
	// Simple Power Law Fit
	int ndata = data_count;
	int ma = 2;
	long double *a;
	a = vector(1, ma);
	// Initial Guess
	a[1] = a_init_guess;	
	a[2] = b_init_guess;
	int *ia;
	ia = ivector(1,2);
	ia[1] = 1;
	ia[2] = 1;
	long double *alamda;
	alamda = (long double *) malloc(sizeof(long double));
	*alamda = -1.0;
	/*
	long double **alpha = malloc(ma * sizeof(long double *));
	for(i = 0; i < ma; i++)
		alpha[i] = malloc(ma*sizeof(long double));

	long double **covar = malloc(ma*sizeof(long double *));
	for(i = 0; i < ma; i++)
		covar[i] = malloc(ma*sizeof(long double));
	*/

	long double **alpha;
	alpha = matrix(1,ma, 1, ma);
	for(i = 1; i <= ma; i++) 
		alpha[i][i] = 1.0;
	long double **covar;
	covar = matrix(1,ma,1,ma);
	long double *chisq;
	chisq = (long double *) malloc(sizeof(long double));
	void (*funcPtr)(long double, long double [], long double *, long double [], int) = &fpowerlaw;

	/* Preparing Iteration */
	int iteration = 0;
	*chisq = 1e7;
	long double chisq_previous;
	long double delta_chisq = 1.0;
	
	printf("\n");
	printf("Pure Powerlaw Fit\n");
	printf("-------------------\n");
	printf("Y = A*x**B\n");
	printf("-------------------\n");

	do
	{
		iteration++;
		chisq_previous = *chisq;
		mrqmin(xdata, ydata, sigma, ndata, a, ia, ma, covar, alpha, chisq, funcPtr, alamda);
		if(*chisq < chisq_previous) delta_chisq = chisq_previous - *chisq;printf("iteration %i, alamda %Lf, chi_old %Lf, chi_new %Lf, delta_chi %Lf, a[1] %Lf, a[2] %Lf \n", iteration, *alamda, chisq_previous, *chisq, delta_chisq, a[1], a[2]);
				
	}while(delta_chisq > 0.01 && *alamda < 1e15);

	printf("\n");
	if(delta_chisq < 0.01) printf("Fit converged after %d iterations (Delta Chi^2 %f)\n", iteration, delta_chisq);
	if(delta_chisq > 0.01)  printf("Fit did not converg. Interrupted after %d iterations (Delta Chi^2 %f, lamda = %f)\n", iteration, delta_chisq, *alamda);
	// After convergence set alamda = 0
	*alamda = 0;
	mrqmin(xdata, ydata, sigma, ndata, a, ia, ma, covar, alpha, chisq, funcPtr, alamda);

	// Output
	printf("\n");
	printf("Chi^2: %Lf, DOF %i, Iterations %i\n", *chisq, ndata, iteration);
	printf("****************\n");
	printf("A = %Lf \t (+/- %Lf)\n", a[1], sqrtl(covar[1][1]));
	printf("B = %Lf \t (+/- %Lf)\n", a[2], sqrtl(covar[2][2]));
	printf("Covariance Matrix\n");
	for (i = 1; i <= ma; i++){
		for(j = 1; j <= ma; j++){
			printf("%Lf\t", covar[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	long double goodnessOfFit = gammq( 0.5*((long double) ndata - ma), (0.5*(*chisq)) );
	printf("Goodness of Fit %Lf\n", goodnessOfFit );
	printf("\n");

/* #################################################################*/

	printf("++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("Powerlaw Fit with corrections (flexible exponent)\n");
	printf("-------------------\n");
	printf("Y = A*X**B  + C*X**(B-1/2))\n");
	printf("-------------------\n");
	ma = 3;
	long double **alphaCorr1;
	long double **covarCorr1;
	alphaCorr1 = matrix(1,ma, 1, ma);
	for(i = 1; i <= ma; i++) 
		alphaCorr1[i][i] = 1.0;
	covarCorr1 = matrix(1,ma,1,ma);
	long double *b1;
	b1 = vector(1, ma);
	// Initial Guess
	b1[1] = a[1];	
	b1[2] = a[2];
	b1[3] = 0.01*a[1];
	int *iia1;
	iia1 = ivector(1,ma);
	iia1[1] = 1;
	iia1[2] = 1;
	iia1[3] = 1;
	iia1[4] = 1;
	*alamda = -1.0;
	*chisq = 1e7;
	
	void (*funcPtrCorr1)(long double, long double [], long double *, long double [], int) = &fpowerlaw_correction_fixed_exponent;
	iteration = 0;

	do
	{
		iteration++;
		chisq_previous = *chisq;
		mrqmin(xdata, ydata, sigma, ndata, b1, iia1, ma, covarCorr1, alphaCorr1, chisq, funcPtrCorr1, alamda);
		if(*chisq < chisq_previous) delta_chisq = chisq_previous - *chisq;
	}while(delta_chisq > 0.01 && *alamda < 1e18);
	
	printf("\n");
	if(delta_chisq < 0.01) printf("Fit converged after %d iterations (Delta Chi^2 %Lf)\n", iteration, delta_chisq);
	if(delta_chisq > 0.01)  printf("Fit did not converg. Interrupted after %d iterations (Delta Chi^2 %Lf, lamda = %Lf)\n", iteration, delta_chisq, *alamda);
	// After convergence set alamda = 0
	*alamda = 0;
	mrqmin(xdata, ydata, sigma, ndata, b1, iia1, ma, covarCorr1, alphaCorr1, chisq, funcPtrCorr1, alamda);

	// Output
	printf("\n");
	printf("Chi^2: %Lf, DOF %d, Iterations %i\n", *chisq, ndata, iteration);
	printf("****************\n");
	printf("A = %Lf \t (+/- %Lf)\n", b1[1], sqrtl(covarCorr1[1][1]));
	printf("B = %Lf \t (+/- %Lf)\n", b1[2], sqrtl(covarCorr1[2][2]));
	printf("C = %Lf \t (+/- %Lf)\n", b1[3], sqrtl(covarCorr1[3][3]));
	printf("Covariance Matrix\n");
	for (i = 1; i <= ma; i++){
		for(j = 1; j <= ma; j++){
			printf("%Lf\t", covarCorr1[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	long double goodnessOfFit1 = gammq( 0.5*((long double) ndata - ma), (0.5*(*chisq)) );
	printf("Goodness of Fit %Lf\n", goodnessOfFit1 );

/* ########################################################## */

	printf("++++++++++++++++++++++++++++++++++++++++++++\n");
	printf("Powerlaw Fit with corrections (flexible exponent)\n");
	printf("-------------------\n");
	printf("Y = A*X**B  + C*X**D)\n");
	printf("-------------------\n");
	ma = 4;
	long double **alphaCorr2;
	long double **covarCorr2;
	alphaCorr2 = matrix(1,ma, 1, ma);
	for(i = 1; i <= ma; i++) 
		alphaCorr2[i][i] = 1.0;
	covarCorr2 = matrix(1,ma,1,ma);
	long double *b2;
	b2 = vector(1, ma);
	// Initial Guess
	b2[1] = b1[1];	
	b2[2] = b1[2];
	b2[3] = b1[3];
	b2[4] = b1[2]-0.5;
	int *iia2;
	iia2 = ivector(1,4);
	iia2[1] = 1;
	iia2[2] = 1;
	iia2[3] = 1;
	iia2[4] = 1;
	*alamda = -1.0;
	*chisq = 1e7;
	
	void (*funcPtrCorr2)(long double, long double [], long double *, long double [], int) = &fpowerlaw_correction_variable_exponent;
	iteration = 0;

	do
	{
		iteration++;
		chisq_previous = *chisq;

		mrqmin(xdata, ydata, sigma, ndata, b2, iia2, ma, covarCorr2, alphaCorr2, chisq, funcPtrCorr2, alamda);
		
		if(*chisq < chisq_previous) delta_chisq = chisq_previous - *chisq;
		

		//printf("iteration %i, alamda %f, chi_old %f, chi_new %f, delta_chi %f, a[1] %f, a[2] %f \n", iteration, *alamda, chisq_previous, *chisq, delta_chisq, a[1], a[2]);
				
	}while(delta_chisq > 0.01 && *alamda < 1e18);
	printf("\n");
	if(delta_chisq < 0.01) printf("Fit converged after %d iterations (Delta Chi^2 %Lf)\n", iteration, delta_chisq);
	if(delta_chisq > 0.01)  printf("Fit did not converg. Interrupted after %d iterations (Delta Chi^2 %Lf, lamda = %Lf)\n", iteration, delta_chisq, *alamda);
	// After convergence set alamda = 0
	*alamda = 0;
	mrqmin(xdata, ydata, sigma, ndata, b2, iia2, ma, covarCorr2, alphaCorr2, chisq, funcPtrCorr2, alamda);

	// Output
	printf("\n");
	printf("Chi^2: %Lf, DOF %d, Iterations %i\n", *chisq, ndata, iteration);
	printf("****************\n");
	printf("A = %Lf \t (+/- %Lf)\n", b2[1], sqrtl(covarCorr2[1][1]));
	printf("B = %Lf \t (+/- %Lf)\n", b2[2], sqrtl(covarCorr2[2][2]));
	printf("C = %Lf \t (+/- %Lf)\n", b2[3], sqrtl(covarCorr2[3][3]));
	printf("D = %Lf \t (+/- %Lf)\n", b2[4], sqrtl(covarCorr2[4][4]));
	printf("Covariance Matrix\n");
	for (i = 1; i <= ma; i++){
		for(j = 1; j <= ma; j++){
			printf("%Lf\t", covarCorr2[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	long double goodnessOfFit2 = gammq( 0.5*((long double) ndata - ma), (0.5*(*chisq)) );
	//long double goodnessOfFit2 = 0.0;
	printf("Goodness of Fit %Lf\n", goodnessOfFit2 );

	// LaTeX- Output
	printf("\nLatex-output (table)\n");
	printf("Range = [%Lf, %Lf] \\\\\n \\hline ", data_min, data_max);
	printf("A = %Lf $( \\pm %Lf)$ & B = %Lf $( \\pm %Lf)$ & & & Q = %Lf\\\\ \n", a[1], sqrtl(covar[1][1]), a[2], sqrtl(covar[2][2]), goodnessOfFit);
	printf("A = %Lf $( \\pm %Lf)$ & B = %Lf $( \\pm %Lf)$ & C = %Lf $( \\pm %Lf)$ & & Q = %Lf\\\\ \n", b1[1], sqrtl(covarCorr1[1][1]), b1[2], sqrtl(covarCorr1[2][2]), b1[3], sqrtl(covarCorr1[3][3]), goodnessOfFit1);
	printf("A = %Lf $( \\pm %Lf)$ & B = %Lf $( \\pm %Lf)$ & C = %Lf $( \\pm %Lf)$ & D =  %Lf $( \\pm %Lf)$ & Q = %Lf\\\\ \n", b2[1], sqrtl(covarCorr2[1][1]), b2[2], sqrtl(covarCorr2[2][2]), b2[3], sqrtl(covarCorr2[3][3]), b2[4], sqrtl(covarCorr2[4][4]), goodnessOfFit2);

	printf("\n Gnuplot-output \n");
	printf("plot [%Lf: %Lf] '%s' u %i:%i:%i w e, %Lf*x**(%Lf) title 'A*x**B', %Lf*x**(%Lf) + (%Lf)*x**(%Lf - 0.5) title 'A*x**B + C*x**(B-1/2)', %Lf*x**(%Lf) + (%Lf)*x**(%Lf) title 'A*X**B+C*X**D'\n", data_min, data_max, filename, x_col, y_col, sigma_col, a[1], a[2], b1[1], b1[2], b1[3], b1[2], b2[1], b2[2], b2[3], b2[4]);

}
	return 0; 
}

void mrqmin(long double x[], long double y[], long double sig[], int ndata, long double a[], int ia[], int ma, long double **covar, long double **alpha, long double *chisq, void (*funcPtr)(long double, long double[], long double *, long double [], int), long double *alamda){
	/* 
	Function parameters
	=================================
	x[1..ndata]	-	x-data array
	y[1..ndata]	-	y-data array
	sig[1..ndata]	-	standard deviation of y-data
	a[1..ma]	- Fit parameters (A, B, C, D..)
	ia[1..ma]	-	set '0' or '1' depending on whether parameter should be used for fit
	covar[1..ma][1..ma]	-	Covariance matrix
	alpha[1..ma][1..ma]	-	Marquardt-matrix
	*funcPtr(x, y, yfit, dyda, ma)	-	Fitting function and derivatives
	yfit	- fit function
	dyda[1..ma]	-	derivatives of fit function
	alamda	-	Marquardt-scalar 
	*/

	void covsrt(long double **covar, int ma, int ia[], int mfit);
	void gaussj(long double **a, int n, long double **b, int m);
	void mrqcof(long double x[], long double y[], long double sig[], int ndata, long double a[], int ia[], int ma, long double **alpha, long double beta[], long double *chisq, void (*funcs)(long double, long double [], long double *, long double [], int));
	
	int j, k, l;
	static int mfit;
	static long double ochisq, *atry, *beta, *da, **oneda;
	// Initialise (alamda initially set < 0.)
	if(*alamda < 0.0){		
		atry = vector(1, ma);
		beta = vector(1, ma);
		da = vector(1, ma);
		for (mfit = 0, j = 1; j <=ma; j++)
			if (ia[j]) mfit++;
		oneda = matrix(1, mfit, 1, 1);
		*alamda = 0.001;
		mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta, chisq, funcPtr); // Insert functions here
		
		ochisq = (*chisq);
		for(j = 1; j <= ma; j++) atry[j] = a[j];
	}
	// Fitting matrix is 'pumped up' on the diagonals by alamda
	for(j = 1; j <= mfit; j++) {
		for (k = 1; k <= mfit; k++) covar[j][k] = alpha[j][k];
			covar[j][j] = alpha[j][j] * (1.0 + (*alamda));
			oneda[j][1] = beta[j];
	}
	// Matrix equation is solved
	gaussj(covar, mfit, oneda, 1);
	for(j = 1; j <= mfit; j++) da[j] = oneda[j][1];
	// After convergence, alamda is set to '0'. Evaluate then covariance matrix
	if(*alamda == 0.0){
		covsrt(covar, ma, ia, mfit);
		covsrt(alpha,ma,ia,mfit);
		free_matrix(oneda, 1, mfit, 1, 1);
		free_vector(da, 1, ma);
		free_vector(beta, 1, ma);
		free_vector(atry, 1, ma);
		return;
	}
	// Try new parameter set
	for(j = 0, l = 1; l<=ma; l++)
		if (ia[l]) atry[l] = a[l] + da[++j];
	mrqcof(x, y, sig, ndata, atry, ia, ma, covar, da, chisq, funcPtr); // Insert functions here
	// If success accept new solution
	if (*chisq < ochisq){
		*alamda *= 0.1;
		ochisq = (*chisq);
		for (j =1; j <= mfit; j++) {
			for (k=1; k <= mfit; k++) alpha[j][k] = covar[j][k];
				beta[j] = da[j];
		}
		for(l=1; l<=ma; l++) a[l] = atry[l];
	} // If unsuccessful, increase alamda
	else {
		*alamda *= 10;
		*chisq = ochisq;
	}
}
void fpowerlaw(long double x, long double a[], long double *y, long double dyda[], int na) {
	/* Power Law y(x) = a1 * x ** -a2 */
	*y = a[1] * powl(x, a[2]);
	dyda[1] = powl(x, a[2]);
	dyda[2] = a[1] * logl(x) * powl(x, a[2]);
	
	/*
	*y = 1 + a[1]*x + a[2]*x*x;
	dyda[1] = x;
	dyda[2] = 2*x;
	*/
}

void fpowerlaw_correction_fixed_exponent(long double x, long double a[], long double *y, long double dyda[], int na) {
	// A*t^beta + B^(beta - 1/2)
	*y = a[1] * powl(x, a[2]) + a[3]*powl(x,a[2]-0.5);
	dyda[1] = powl(x, a[2]);
	dyda[2] = a[1] * logl(x) * powl(x, a[2]);
	dyda[3] = powl(x, a[2] - 0.5);
}

void fpowerlaw_correction_variable_exponent(long double x, long double a[], long double *y, long double dyda[], int na) {
	/* Power Law y(x) = a1 * x ** a2 + a3* x ** a4 */
	*y = a[1] * powl(x, a[2])+ a[3]*powl(x,a[4]);
	dyda[1] = powl(x, a[2]);
	dyda[2] = a[1] * logl(x) * powl(x, a[2]);
	dyda[3] = powl(x, a[4]);
	dyda[4] = a[3]*logl(x) * powl(x, a[4]);
	
	/*
	*y = 1 + a[1]*x + a[2]*x*x;
	dyda[1] = x;
	dyda[2] = 2*x;
	*/
}
void mrqcof(long double x[], long double y[], long double sig[], int ndata, long double a[], int ia[], int ma, long double **alpha, long double beta[], long double *chisq, void (*funcPtr)(long double, long double [], long double *, long double [], int)){
	/* Used by mrqmin to evaluate linearised fitting matrix alpha, vector beta, and calculate chi^2 */
	int i, j, k, l, m, mfit = 0;
	long double ymod, wt, sig2i, dy, *dyda;
	dyda = vector(1, ma);
	for(j = 1; j <= ma; j++)
		if (ia[j]) mfit++;
	// Initialise alpha, beta (symmetric)
	for(j=1; j <=mfit; j++) {
		for (k=1; k <= j; k++) alpha[j][k] = 0.0;
			beta[j] = 0.0;
	}
	*chisq = 0.0;
	for(i = 1; i <= ndata; i++) {
		(*funcPtr)(x[i], a, &ymod, dyda, ma);
		sig2i = 1.0/(sig[i] * sig[i]);
		dy = y[i] - ymod;
		for (j = 0, l = 1; l <= ma; l++) {
			if (ia[l]) {
				wt = dyda[l] * sig2i;
				for(j++, k = 0, m=1; m <=l; m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m]; 
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i; // Compute chi^2
	}
	// Fill up symmetric side
	for (j = 2; j <=mfit; j++)
		for(k = 1; k < j; k++) alpha[k][j] = alpha[j][k];
	free_vector(dyda, 1, ma);
}
void covsrt(long double **covar, int ma, int ia[], int mfit){
	int i,j,k;
	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}
void gaussj(long double **a, int n, long double **b, int m){
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	long double big,dum,pivinv,temp;
	indxc=ivector(1,n); 
	indxr=ivector(1,n); 
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big = 0.0;
		for (j=1;j<=n;j++) 
			if (ipiv[j] != 1) 
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1, pivot > 1");
				}
		++(ipiv[icol]);
		
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow; 
		indxc[i]=icol; 
		if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix, alpha zero on diagonal");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
			for (ll=1;ll<=n;ll++)
				if (ll != icol) { 
					dum=a[ll][icol];
					a[ll][icol]=0.0;
					for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
					for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
		}
	}
	
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	} 
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

int count_lines(char *filename, int xcol, int ycol, int scol, long double min, long double max)
{
	/* Counts lines of valid data to then generate xdata, ydata, sigma array in main-loop */
	int i = 0, j = 0, k=0;
	FILE *data;
	
	data = fopen(filename, "r");
	
	if(data == NULL){
		printf("File could not be opened. Check filename\n");
		return;
	}
	


	int lines = 0;
	char line[10000];
	char c, buffer_c;
	
	int valid_data;
	char *buffered_value;
	long double x_value;
	int columns_checked = 0;
	buffered_value = (char *)malloc(20*sizeof(char));
	while(fgets(line, sizeof(line), data))
	{
		valid_data = 0;
		if(!isdigit(*line)) continue;
		else valid_data = 1;

		// At first line, count columns
		if(lines == 0 && valid_data == 1 && columns_checked == 0) 
		{
			columns = 1;
			c = line[0];
			while(c != '\n')
			{
				if(c == '\t') columns++; // For each tab, one column more
				c = line[++i];
			}
			columns_checked = 1;
			// if(xcol > columns || ycol > columns) return -1;
		}

		// Check whether value is within allowed range // How access xcol ??
		if(valid_data)
		{
			// Go to next tab whilst filling char into array. If you encounter tab, save char into long double.
			// 2nd idea go to (n-1)th tab then sscanf to next \t into long double
			c = line[0];
			i = 0;
			j=0;
			while(c != '\n') // Till end of line
			{	
				if(j == xcol-1) // Triggered when just before value
				{

					buffer_c = line[i++]; // Go to first value
					k = 0;
					while(buffer_c != '\t') // Scan till next delimiter
					{
						buffered_value[k] = buffer_c; // Write next char into buffered_value
						buffer_c = line[i++];
						k++;
					}
					x_value = (long double) atof(buffered_value);
					buffered_value = (char *)malloc(30*sizeof(char));

					if(min != -1)
					{
						if(max != -1)
						{
							if(x_value >= min && x_value < max) lines++;
						}
						else
						{
							if(x_value >= min) lines ++;
						}
					}
					else
					{
						if(max != -1)
						{
							if(x_value < max) lines++;
						}
						else
						{
							lines++;
						}
					}
					break; // Go to next line
				}
				c = line[i++];
				if(c == '\t') j++;

			}
		}
	}
	printf("File %s contains %i valid lines with %i columns\n", filename, lines, columns);
// Now read everything into x,y,sigma
	//Initialise Arrays

	return lines;
}

void read_data(char *filename, int xcol, int ycol, int scol, long double min, long double max, long double xdata[], long double ydata[], long double sigma[])
{
	int i = 0, j = 0, k=0;
	FILE *data;
	
	data = fopen(filename, "r");
	
	if(data == NULL){
		printf("File could not be opened. Check filename\n");
		return;
	}
	


	int lines = 0;
	char line[10000];
	char c, buffer_c;
	int valid_row = 1;
	int col = 0;

	int valid_data;
	char *buffered_value;
	long double x_value;

	char **buffered_row;
	buffered_row = (char **)malloc(columns*sizeof(char *));
	for(i =0; i <columns; i++)
	{
		buffered_row[i] = (char *)malloc(30*sizeof(char));
	}

	if(!buffered_row) return;
	rewind(data);

	while(fgets(line, sizeof(line), data))
	{
		valid_data = 0;
		if(!isdigit(*line)) continue;
		else valid_data = 1;

		// At first line, count columns
		if(lines == 0 && valid_data == 1) 
		{
			columns = 1;
			c = line[0];
			while(c != '\n')
			{
				if(c == '\t') columns++; // For each tab, one column more
				c = line[++i];
			}
			// if(xcol > columns || ycol > columns) return -1;
		}

		// Check whether value is within allowed range // How access xcol ??
		if(valid_data)
		{
			// Load entire line into buffer
			c = line[0];
			i = 0;
			j=0;
			col = 0;

			while(c != '\n') // Till end of line
			{	
				k = 0;
				
				while( c != '\t' && c != '\n')
				{
					buffered_row[col][k++] = c;
					c = line[++i];
				}
				//printf("#ping read_data '%s', c='%c', col = %i / %i\n", buffered_row[col], c, col, columns);
				col++;
				if(c != '\n') c = line[++i];
			}

			x_value = (long double) atof(buffered_row[xcol-1]);

			if(min != -1)
					{
						if(max != -1)
						{
							if(x_value >= min && x_value < max)
							{
								xdata[valid_row] = (long double) atof(buffered_row[xcol-1]);
								ydata[valid_row] = (long double) atof(buffered_row[ycol-1]);
								if (scol != -1) sigma[valid_row] = (long double) atof(buffered_row[scol-1]);
								valid_row++;
							}
						}
						else
						{
							if(x_value >= min) 
							{
								xdata[valid_row] = (long double) atof(buffered_row[xcol-1]);
								ydata[valid_row] = (long double) atof(buffered_row[ycol-1]);
								if (scol != -1) sigma[valid_row] = (long double) atof(buffered_row[scol-1]);
								valid_row++;
							}
						}
					}
					else
					{
						if(max != -1)
						{
							if(x_value < max) 
							{
								xdata[valid_row] = (long double) atof(buffered_row[xcol-1]);
								ydata[valid_row] = (long double) atof(buffered_row[ycol-1]);
								if (scol != -1) sigma[valid_row] = (long double) atof(buffered_row[scol-1]);
								valid_row++;
							}
						}
						else
						{
							xdata[valid_row] = (long double) atof(buffered_row[xcol-1]);
							ydata[valid_row] = (long double) atof(buffered_row[ycol-1]);
							if (scol != -1) sigma[valid_row] = (long double) atof(buffered_row[scol-1]);
							valid_row++;
						}
					}
			
			// Reset buffered row
			for(i =0; i <columns; i++)
			{
				buffered_row[i] = (char *)malloc(30*sizeof(char));
			}
		}
	}
}



void printhelp(){
	printf("# PowerFit \n# Levenberg-Marquardt-Method \n# (From Numerical Recipes in C: The Art of Scientific Computing)\n# Fits Power Law y(x) = A * x ** -B to data of format (x_i, y_i, sigma_i). \n# Returns parameters, deviations, and goodness-of-fit.\n# Syntax: ./powerfit [Options] input_data.dat \n");     
	
	printf("Correct Syntax: ./powerfit -m {minimum x-value | 0} -M {maximal x-value | infty} -x {xcol | 1} -y {ycol | 2} -s {sigma column | 3} -A {Initial Guess A | 1.0} -B {Initial Guess B | 2.0}\n");
}

long double gammq(long double a, long double x)
{
	void gcf(long double *gammcf, long double a, long double x, long double *gln);
	void gser(long double *gamser, long double a, long double x, long double *gln);
	void nrerror(char error_text[]);
	long double gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) {printf("Invalid arguments in routine gammq x = %Lf, a = %Lf\n", x, a); nrerror("Invalid arguments in routine gammq\n");}
	if (x < (a+1.0)){
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} 
	else {
	gcf(&gammcf,a,x,&gln);
	return gammcf;
	}
}

void gser(long double *gamser, long double a, long double x, long double *gln){
	long double gammln(long double xx);
	void nrerror(char error_text[]);
	int n;
	long double sum,del,ap;
	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else{
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
		++ap;
		del *= x/ap;
		sum += del;
		if (fabsl(del) < fabsl(sum)*EPS) {
		*gamser=sum*expl(-x+a*logl(x)-(*gln));
		return;
		}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}



void gcf(long double *gammcf, long double a, long double x, long double *gln)
{
	long double gammln(long double xx);
	void nrerror(char error_text[]);
	int i;
	long double an,b,c,d,del,h;
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabsl(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabsl(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabsl(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=expl(-x+a*logl(x)-(*gln))*h;
}

long double gammln(long double xx)
{
	long double x,y,tmp,ser;
	static long double cof[6]={76.18009172947146,-86.50532032941677,
	24.01409824083091,-1.231739572450155,
	0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*logl(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+logl(2.5066282746310005*ser/x);
}