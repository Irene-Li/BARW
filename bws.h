#pragma once
#define MAXSIZE 5000000

//default resolution for histogram intervals
#define BINS 2000

#define LATTICE_BASE 3
//dummy change 2
#define TRUE 1
#define FALSE 0
//#a utils  wrapper for random number stuff would be good
#define RNG_MT_BITS (32)
#define RNG_MAX (~0UL)
#define RNG_TYPE unsigned long
#define RNG_MT_MAX (4294967295)
#define RNG_TYPE_OUT "%lu"

#define TRACE_FLAG (1<<2) // = 4 
#define CURRENT_FLAG (1<<3) // = 8 
#define IMMOBILE_FLAG (1<<4) // = 32 
#define ADD_FLAG ( TRACE_FLAG | CURRENT_FLAG ) 
//#define HOLE (1<<6)//porous in lattice  - reserved!
#define BOUNDARY (1<<7)//check


#define TIME_WRITE_COUNTS 10 // e.g from 0..10 write 10 times, from 100..1000 write 10 times etc.
#define TIME_WRITE_SCALES  9
//#define MAX_T TIME_WRITE_COUNTS * TIME_WRITE_SCALES 
#define MAX_T 20000000000.
#define MIN_T  0.1
#define MAX_MOMENTS 6
//#define CHUNKS_SIZE 100000

struct simple_stack
{
	//int wait[MAXSIZE];//experiment
	int stk[MAXSIZE];
	int top;
};
typedef struct simple_stack SSTACK;




extern char* lattice;
void run_for_realisations(int N, int L, int D, double h, double p, double q, double a, int bcs, int seed);

//Required for RNG correctness
/*void init_genrand(RNG_TYPE);
RNG_TYPE genrand_int32(void);
double genrand_real1(void);//[0,1]
double genrand_real2(void);//[0,1)
*/
void printhelp(void);

int parse_args(int *bcs, int *D, int *L, int *C, int *Ln, int *N, int *seed, int *min_l, int *max_l, double *h, double *p, double *q, double *a, int argc, char *argv[]);
