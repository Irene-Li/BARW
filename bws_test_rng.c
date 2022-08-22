/*Branching Wiener Sausage - see readme for details and /scripts for compilation and execution notes*/

//todo document possibly subtle things e.g. we use flags on forloops in writing to decide when to write out things like histograms, lattices, avalanches etc. see for loops
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined linux || defined __APPLE__
#include <getopt.h>
#include <unistd.h>
#endif
#include "rng/mt19937ar_clean_bkp.h"
#include "bws.h"
#include "lattice_core.h"

//write the trace
#define WRITE_TRACE_SIZE(N, L, t) if (__verbose==1 && __trace > 1){printf("#MAX_TRACE(%.3Lf,%d,%d):%d\n", t, L, N,__trace);}
#define PRINT_VERSION {printf(VERSION); printf("\n");}
#define PRINT_PARAMS printf("#Parameters={ Seed: %i, Branching rate: %g, p : %g, q : %g, Realisations: %i, Chunk size: %i, Dimension: %i, (Max) Lattice size: %i, Graph Type: %d }\n", seed, h, p, q, N, CHUNK_SIZE, D, MAX_L,__graph_type__);
#define PRINT_OPTIONS printf("#TODO\n");
//we have predetermined write times at 1,2,3,...9,10,20,....90,100,200,.....900,1000,2000,...
#define WRITE_TIME_FUNCTION(i) (MIN_T * pow(MAX_T/MIN_T, (double)i/(BINS-1)))
#define INIT_WRITE_TIMES {int i=0; for(i =0;i < BINS; i++){write_times[i] = WRITE_TIME_FUNCTION(i);}}
//prepare an array to store a histogram of for trace counts i each sample path - it is assumed to be between 1 and L^2 for practicality - check?
#define INIT_TRACE_HISTOGRAM(L) {int l;  for (l=0;l< BINS*write_hist; l++ ){__trace_histogram[l]=0;}}
//prepare and set to 0 the avalanche matrix
#define INIT_AVALACNHES {int l, m; 	for (m = 0; m <= MAX_MOMENTS; m++) {for (l = 0; l < 10; l++) {avalanche_moments[m] = 0;}}}

//write the header for the tabular data - make sure it matches what is in COMMIT_MOMENTS_TRACE
#define WRITE_HEADER {int i;for(i = 0; i <= MAX_MOMENTS; i++){printf("\tM%d",i);} printf("\n");}
//SA:The "Ben Index" is used to map a trace value to a bin - MAXSIZE is currently set to a large number which is something like the expected max value (should maybe change)
//#define BIN_TRACE(t) {__trace_histogram[(int)floor(BINS*log(t)/log(MAXSIZE))]++;}
//SA temp checking if this happens an we need to do something about it
#define BIN_TRACE(t) {int i; i = (int)floor(BINS*log(t)/log(MAXSIZE)); if (i < BINS) {__trace_histogram[i]++;}}
//flush out the histogram in a comment
#define WRITE_TRACE_HIST(L) {printf("#TRACE_HIST_L=%i: ",L); int l; for(l=0;l<BINS*write_hist; l++){printf("%d ", __trace_histogram[l]); } printf("\n");}
/* This is such that numbers don't explode. Else, they go with N^n. Rather cut them into pieces here, and then blow them up later */
// At the end of each run in system size L write the final avalanche size into avalanche_moments.
#define COMPUTE_AVALANCHE_STATS(av,N){int i; long double normalised_avalanche = av / ((double)N); for (i = 0; i <= MAX_MOMENTS*write_avalanches; i++) {avalanche_moments[i] += (long double)pow(normalised_avalanche, i);}}
#define WRITE_AVALANCHES(L) {printf("#AVALANCHES_L=%i: ", L); int i; 	for (i = 1; i <= MAX_MOMENTS*write_avalanches; i++) {printf("\t%.1Lf", avalanche_moments[i-1]); } printf("\n");}
//todo: above we should be passng size

#define _RANDOM_INT(r, n) unsigned x= genrand_int32(); while (x >= RNG_MT_MAX - (RNG_MT_MAX % n)) {x = genrand_int32();} x %= n; r = (int)x;
#define RANDOM_DOUBLE  genrand_real2()
#define RANDOM_ORIENTATION (genrand_int32() % (2*D))
#define EXP_WAIT(n) (1.0/(double)n) * (-log(1-RANDOM_DOUBLE))

// #define POP(p)  p=stack.stk[--stack.top];
// #define PUSH(p) stack.stk[stack.top++] = p;
#define PARTICLE_COUNT stack.top
//todo add a stack with a certain capacity and record trace on it unless we exceed threshold in which case?
#define UPDATE_TRACE(p) __trace++; log_trace(p); if (__compute_rog == 1) {radius_gyr += distance_from_center(p); if (__trace > 1){radius_gyr /= (__trace-1);}}
#define ADD(p) if(TRACE_FLAG!=(lattice[p] & TRACE_FLAG)) {UPDATE_TRACE(p)}; lattice[p] |= ADD_FLAG; PUSH(p, &stack); 
#define REMOVE(p) lattice[p] &= ~CURRENT_FLAG
#define MOVE(from,to) REMOVE(from); ADD(to);
#define STAY(p) ADD(p)
/* pruess 18 July 2021: I think there needs to be a log_trace(p) here, too. */
#define TRY_ADD_IMMOBILE(p)  if(IMMOBILE_FLAG!=(lattice[p] & IMMOBILE_FLAG))__immobileTrace++ ; lattice[p] |= IMMOBILE_FLAG; log_trace(p)

/*Function declaration */
void print_write_times(void);
/*globals*/
int __sample_n__ = -1;//541; //-1;
int write_hist = 0, write_lattice = 0, write_image = 1, write_msd = 1, write_edge = 0, write_hull = 0, write_total = 0, write_final = 1, write_coarse_grain_moments = 1; 
int write_avalanches = 0, write_moments = 1, write_edge_reach = 1; 
int branch_method = 2; 

char *lattice;
SSTACK stack;
SSTACK past_pos_stack; 
//temp - idea here is we can configure the resolution of things like write times and historgrams but then we need to malloc
int BIN_N = BINS;

long double radius_gyr = 0;
// Avalanche (Avalanche statistics is highly not optimised, don't take it serious yet)
long double avalanche_moments[MAX_MOMENTS + 1]; // Here, I'm assuming that there will be 10 different system sizes to be checked!
long double tracer_moments[BINS][MAX_MOMENTS + 1];//NB: adding 1 to number of moments to capture 0th moment!! ALL FOREACH shoud terminate at <=MAX_MOMENTS
long double write_times[BINS]; // bw: changed to double
long double active_moments[BINS][MAX_MOMENTS + 1];
int __trace = 0, __immobileTrace = 0, __maxParticles = 0, __verbose = 0, __compute_rog = 0;
int __trace_histogram[BINS];
int CHUNK_SIZE = 0;
int __graph_type__ = 0;
//this can be removed - is function of something else
int sizeCount = 0;
int MAX_L = 0;

int main()
{
	double rd = 0.0;
	int a = 0;

	for (a = 0; a < 100000; a++) {
		rd = RANDOM_DOUBLE;
		printf("%.3Lf \n", rd);
	}
}