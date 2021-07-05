#pragma once
#define MAXSIZE 500000

//#a utils  wrapper for random number stuff would be good
#define RNG_MT_BITS (32)
#define RNG_MAX (~0UL)
#define RNG_TYPE unsigned long
#define RNG_MT_MAX (4294967295)
#define RNG_TYPE_OUT "%lu"

#define TRACE_FLAG (1<<2)
#define CURRENT_FLAG (1<<3)
#define IMMOBILE_FLAG (1<<4)
#define ADD_FLAG ( TRACE_FLAG | CURRENT_FLAG )
#define BOUNDARY (1<<7)//check


#define TIME_WRITE_COUNTS 10 // e.g from 0..10 write 10 times, from 100..1000 write 10 times etc.
#define TIME_WRITE_SCALES  8
#define MAX_T TIME_WRITE_COUNTS * TIME_WRITE_SCALES 
#define MAX_MOMENTS 8
#define CHUNKS_SIZE 10000

struct simple_stack
{
	//int wait[MAXSIZE];//experiment
	int stk[MAXSIZE];
	int top;
};
typedef struct simple_stack SSTACK;

extern char* lattice;
void run_for_realisations(int N, int L, int D, double h, double sigma, int bcs);

//Required for RNG correctness
void init_genrand(RNG_TYPE);
RNG_TYPE genrand_int32(void);
double genrand_real1(void);//[0,1]
double genrand_real2(void);//[0,1)

void printhelp(void);


