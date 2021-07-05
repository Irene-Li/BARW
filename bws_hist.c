/* GP comment, meeting 16/03
*Branching Wiener Sausage - see readme for details and /scripts for compilation and execution notes*/

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
#include "lattice_core.h"
#include "bws.h"

//write the trace
#define WRITE_TRACE_SIZE(N, L, t) if (__verbose==1 && __trace > 1){printf("#MAX_TRACE(%.3Lf,%d,%d):%d\n", t, L, N,__trace);}
#define PRINT_VERSION {printf(VERSION); printf("\n");}
#define PRINT_PARAMS printf("#Parameters={ Seed: %i, Hopping rate: %g, Realisations: %i, Chunk size: %i, Dimension: %i, (Max) Lattice size: %i, Graph Type: %d }\n", seed, h, N, CHUNK_SIZE, D, MAX_L,__graph_type__);
#define PRINT_OPTIONS printf("#TODO\n");
//we have predetermined write times at 1,2,3,...9,10,20,....90,100,200,.....900,1000,2000,...
#define WRITE_TIME_FUNCTION(i) (MIN_T * pow(MAX_T/MIN_T, (double)i/(BINS-1)))
//#define INIT_WRITE_TIMES {int i=0; int a =0; int b=0; for(a =-2; a < TIME_WRITE_SCALES-2; a++){for (b = 1; b <= TIME_WRITE_COUNTS; b++){ write_times[i++] = (b)*powl(10,a); }}}
#define INIT_WRITE_TIMES {int i=0; for(i =0;i < BINS; i++){write_times[i] = WRITE_TIME_FUNCTION(i);}}
//moments for just the trace (todo add others) we note that MAX_MOMENTS does not incude the zero'th moment so we add MAX_MOMENTS+1 - this is just a convention
#define INIT_MOMENTS_TRACE {int t; int m; for(m = 0; m <= MAX_MOMENTS; m++){for(t = 0; t < BINS; t++){trace_moments[t][m]=0;}}; INIT_WRITE_TIMES; }
//prepare an array to store a histogram of for trace counts i each sample path - it is assumed to be between 1 and L^2 for practicality - check?
#define INIT_TRACE_HISTOGRAM(L) {int l;  for (l=0;l< BINS*write_hist; l++ ){__trace_histogram[l]=0;}}
//prepare and set to 0 the avalanche matrix
#define INIT_AVALACNHES {int l, m; 	for (m = 0; m <= MAX_MOMENTS; m++) {for (l = 0; l < 10; l++) {avalanche_moments[m] = 0;}}}
//write each recorded obs(t) - in this case just one - extend with an offset for other moments e.g. 1-8 for trace, 1-8 for population count
#define WRITE_MOMENTS_TRACE(x,t) {int i; for(i = 0; i <= MAX_MOMENTS; i++){trace_moments[t][i] += (long double)pow(x,i);}}
//write the header for the tabular data - make sure it matches what is in COMMIT_MOMENTS_TRACE
#define WRITE_HEADER {int i;printf("chunk\tL\tD\tBCs\tt"); for(i = 0; i <= MAX_MOMENTS; i++){printf("\tM%d",i);} printf("\n");}
//flush out moments to standard out - add some header info that will be used in tabular data handling
#define COMMIT_MOMENTS_TRACE(c,L,D,BCs){ int t; int m; for(t = 0; t < BINS; t++){printf("%d\t%d\t%d\t%d\t%.3Lf",c,L,D,BCs,write_times[t]); for(m = 0; m <= MAX_MOMENTS; m++){printf("\t%.1Lf", trace_moments[t][m]);} printf("\n");}; } INIT_MOMENTS_TRACE
//SA:The "Ben Index" is used to map a trace value to a bin - MAXSIZE is currently set to a large number which is something like the expected max value (should maybe change)
#define BIN_TRACE(t) {__trace_histogram[(int)floor(BINS*log(t)/log(MAXSIZE))]++;}
//flush out the histogram in a comment
#define WRITE_TRACE_HIST(L) {printf("#TRACE_HIST_L=%i: ",L); int l; for(l=0;l<BINS*write_hist; l++){printf("%d ", __trace_histogram[l]); } printf("\n");}
/* This is such that numbers don't explode. Else, they go with N^n. Rather cut them into pieces here, and then blow them up later */
// At the end of each run in system size L write the final avalanche size into avalanche_moments.
#define COMPUTE_AVALANCHE_STATS(av,N){int i; long double normalised_avalanche = av / ((double)N); for (i = 0; i <= MAX_MOMENTS*write_avalanches; i++) {avalanche_moments[i] += (long double)pow(normalised_avalanche, i);}}
#define WRITE_AVALANCHES(L) {printf("#AVALANCHES_L=%i: ", L); int i; 	for (i = 1; i <= MAX_MOMENTS*write_avalanches; i++) {printf("\t%.1Lf", avalanche_moments[i-1]); } printf("\n");}
//todo: above we should be passng size

#define _RANDOM_INT(r, n)unsigned x= genrand_int32(); while (x >= RNG_MT_MAX - (RNG_MT_MAX % n)) {x = genrand_int32();} x %= n; r = (int)x;
#define RANDOM_DOUBLE  genrand_real1()
#define RANDOM_ORIENTATION (genrand_int32() % (2*D))
#define POP_RANDOM_PARTICLE(p)int r,top; _RANDOM_INT(r, stack.top); p = stack.stk[r]; stack.stk[r] = POP(top); p;
#define EXP_WAIT(n) (1.0/(double)n) * (-log(1-RANDOM_DOUBLE))

#define POP(p)  p=stack.stk[--stack.top]; p
#define PUSH(p) stack.stk[stack.top++] = p;
#define PARTICLE_COUNT stack.top
//todo add a stack with a certain capacity and record trace on it unless we exceed threshold in which case?
#define UPDATE_TRACE(p) __trace++; log_trace(p); if (__compute_rog == 1) {radius_gyr += distance_from_center(p); if (__trace > 1){radius_gyr /= (__trace-1);}}
#define ADD(p) if(TRACE_FLAG!=(lattice[p] & TRACE_FLAG))UPDATE_TRACE(p); lattice[p] |= ADD_FLAG; PUSH(p); //push(&stack,u);
#define REMOVE(p) lattice[p] &= ~CURRENT_FLAG
#define MOVE(from,to) REMOVE(from); ADD(to);
#define STAY(p) ADD(p)
#define TRY_ADD_IMMOBILE(p)  if(IMMOBILE_FLAG!=(lattice[p] & IMMOBILE_FLAG))__immobileTrace++ ;lattice[p] |= IMMOBILE_FLAG

/*Function declaration */
void print_write_times(void);
/*globals*/
int __sample_n__ = -1;//541; //-1;
int write_hist = 0, write_lattice = 0, write_image = 0, write_hull = 0, write_avalanches = 0;

char *lattice;
SSTACK stack;
//temp - idea here is we can configure the resolution of things like write times and historgrams but then we need to malloc
int BIN_N = BINS;

long double radius_gyr = 0;
// Avalanche (Avalanche statistics is highly not optimised, don't take it serious yet)
long double avalanche_moments[MAX_MOMENTS + 1]; // Here, I'm assuming that there will be 10 different system sizes to be checked!
long double trace_moments[BINS][MAX_MOMENTS + 1];//NB: adding 1 to number of moments to capture 0th moment!! ALL FOREACH shoud terminate at <=MAX_MOMENTS
long double write_times[BINS]; // bw: changed to double
int __trace = 0, __immobileTrace = 0, __maxParticles = 0, __verbose = 0, __compute_rog = 0;
int __trace_histogram[BINS];
int CHUNK_SIZE = 0;
int __graph_type__ = 0;
//this can be removed - is function of something else
int sizeCount = 0;
int MAX_L = 0;

void spit_out_image(int L, long double t, int n, int new_lines) {
	int a = 0;
	printf("#LATTICE(%.3Lf,%d, %d): ", t, L, n);
	for (a = 0; a < L*L; a++) {
		if (a % L == 0 && new_lines == TRUE)
			printf("\n");
		if (TRACE_FLAG == (lattice[a] & TRACE_FLAG)) {
			printf("%03d ", a);
		}
		else {
			//this is for testing and viz
			//printf("    ");
		}
	}
	printf("\n");
}

int main(int argc, char *argv[])
{
	// Set buffer to zero so that process has to write into stdout immediately
	#if defined linux || defined __APPLE__
		setlinebuf(stdout);
	#endif

	double sigma, h;

	int BCs, D, l, Ln, L, C, N, seed, min_l, max_l;
	if (parse_args(&BCs, &D, &L, &C, &Ln, &N, &seed, &min_l, &max_l, &h, argc, argv) == TRUE) {
		init_genrand(seed);
		CHUNK_SIZE = (int)(N / C);
		sigma = (1 - h) / 2;
		PRINT_PARAMS
		printf("#Version: ");
		printf("#");
		#if defined linux || defined __APPLE__
				PRINT_VERSION
		#endif

		allocate_lattice(&lattice, MAX_L, D, __graph_type__);
		INIT_WRITE_TIMES
		WRITE_HEADER;
		//we either run a single L if specified, otherwise go through the motions
		if (L > 0) { run_for_realisations(N, L, D, h, sigma, BCs); }
		else {
			for (l = min_l; l <= max_l; l++) {
				L = pow(2, l) - 1;
				run_for_realisations(N, L, D, h, sigma, BCs);
			}
		}
	}
}

inline void run_for_realisations(int N, int L, int D, double h, double sigma, int bcs) {
	printf("# Running for L = %i\n", L);
	INIT_MOMENTS_TRACE;
	INIT_TRACE_HISTOGRAM(L);
	INIT_AVALACNHES

	double __avalanche = 0; // where should this be? changing SA for per run
	int _n = 0, chunk = 0;
	for (_n = 0; _n < N; _n++) {
		double rd = 0.0;
		long double time = 0.0, last_time = 0.0;
		int wait_for_for_a_big_one = 1000;
		int write_time_index = 0, pos = 0, ro = 0, next = 0;
		int lattice_flush_index = 0;
		init_lattice(bcs, L, D);
		stack.top = 0; __trace = 0; __immobileTrace = 0, __maxParticles = 0; __avalanche = 0.0;

		int cen = get_center(L, D);
		ADD(cen);
		// kick out this do while loop

		do {
			time += (EXP_WAIT(PARTICLE_COUNT));

			__avalanche += (time - last_time)/*dt*/ * ((double)PARTICLE_COUNT);

			while ((write_times[write_time_index] < time) && (write_time_index <= BINS - 1)) 			{
				WRITE_MOMENTS_TRACE(__trace, write_time_index);write_time_index++;
/*				if (__trace > 20) {
					//write the hull IF we have set a particular global option
					if ((write_time_index % 10 == 0) &&(__trace>200) )//(_n == __sample_n__))
					{
						if (write_image == 1) spit_out_image(L, time, _n, FALSE);
						//if (write_hull == 1) right_turning_walk(L, _n, time);
						wait_for_for_a_big_one += 500;
					}
*/				
			}
			
			if (write_times[BINS - 1] < time) 
			{ 
				printf("##Recorded excessive time"); 
				break; 
			}//todo see why this sometimes happens

			POP_RANDOM_PARTICLE(pos);
			rd = RANDOM_DOUBLE;//to choose sub-process...

			if (rd <= sigma) {//branch locally (and stay)
				STAY(pos);
				ADD(pos);
			}

			else if (rd <= sigma + h) { //hop
				ro = RANDOM_ORIENTATION;
				next = diffuse(pos, ro);

				while (next == -2)
					next = diffuse(pos, RANDOM_ORIENTATION);//reflected from holes
				if (next == -1) {//-1 illegal - dead for open boundary - not put back on stack, reset flag on lattice
					REMOVE(pos);
					continue;
				}
				MOVE(pos, next);
			}
			else {//extinction
				REMOVE(pos);
				continue;
			}

			//if (PARTICLE_COUNT > __maxParticles)__maxParticles = PARTICLE_COUNT;
			last_time = time;
		}while (PARTICLE_COUNT); // do..while avalanche alive
		//WRITE_TRACE_SIZE(_n, L, last_time)

			//pad gaps in time for consistent data shape
			// comment that out
		//while (write_time_index < BINS) { WRITE_MOMENTS_TRACE(__trace, write_time_index);write_time_index++; }
		//printf("%f", last_time);
		//BIN_TRACE(__trace);
		//choose N = K*CHUNK_SIZE + 1

		//if ((_n + 1) % CHUNK_SIZE == 0 || (_n + 1) == N) { chunk++; COMMIT_MOMENTS_TRACE(chunk, L, D, bcs); }
		//COMPUTE_AVALANCHE_STATS(__avalanche, N)
		printf("#TRACE %i\n", __trace);
	} //for loop
	//WRITE_TRACE_HIST(L);
	//WRITE_AVALANCHES(L);
	//printf("#ROG(%.3Lf,%d, %d):%.3Lf\n", time, L, N, radius_gyr);
	printf("#okely dokely!");//look for this line int stats out
}


int parse_args(int *bcs, int *D, int *L, int *C, int *Ln, int *N, int *seed, int *min_l, int *max_l, double *h, int argc, char *argv[]) {
	// Define default parameters
	*seed = 5; *N = 5000000; *D = 2; *Ln = -1; *bcs = 0; *L = -1; *C = 1; *h = 0.1;
	*min_l = 2; *max_l = 7;

	//test///////////////
	//*L = 243;//
	//__graph_type__ = 0;
	/////////////////////
	// Check if Windows. Then get parameters as "seed, N, Ln (lattice length), D, bcs"
#ifdef _WIN64
	int temp = time(NULL);
	if (argc > 1 && atoi(argv[1]) != -1)
		temp = atoi(argv[1]);
	if (argc == 3) {
		*seed = temp; *N = atoi(argv[2]);
	}

	else if (argc == 6) {
		*seed = temp; *N = atoi(argv[2]); *Ln = atoi(argv[3]); *D = atoi(argv[4]); *bcs = atoi(argv[5]);
	}
	else {
		printf("#Using default arguments - pass either seed and N or all args (seed,N, Ln, D, Bcs)\n");
	}
#endif

	// Check for linux, if so use getopt
#if defined linux || defined __APPLE__
//	int ch;
//	while ((ch = getopt(argc, argv, "S:N:L:D:B:h")) != -1)  // bw change 'R'
//		switch (ch) {
//		case 'S':// seed
//			*seed = atof(optarg);
//			break;
//		case 'N'://realisations
//			*N = atoi(optarg);
//			break;
//		case 'L':
//			*Ln = atoi(optarg);
//			break;
//		case 'D'://dimensions
//			*D = atoi(optarg);
//			break;
//		case 'B':
//			*bcs = atoi(optarg);
//			break;
//		case 'h':
//			printhelp();
//			_exit(0);
//		default:
//			break;
//		}

	int             c;
	int option_index = 0;
	const char    * short_opt = "C:N:L:D:h:";
	struct option   long_opt[] =
	{
	   {"help",          no_argument,       NULL, 0},
	   {"verbose",          no_argument,       NULL, 0 },
	   {"seed",          required_argument, NULL, 0},
	   {"D",          required_argument, NULL, 0},
	   {"L",          required_argument, NULL, 0},
	   {"graph",          required_argument, NULL, 0 },
	   {"specimen",          required_argument, NULL, 0 },
	   {"Ln",          required_argument, NULL, 0},
	   {"BCs",          required_argument, NULL, 0},
	   {"minLn",          required_argument, NULL, 0},
	   {"maxLn",         required_argument, NULL, 0},
	   {"wROG",          no_argument, NULL, 0 },
	   {"wHist",          no_argument, NULL, 0},
	   {"wAval",          no_argument, NULL, 0},
	   {"wImage",          no_argument, NULL, 0},
	   {"wHull",          no_argument, NULL, 0 },
	   {"write_times",          no_argument, NULL, 0},
	   {"version", 			no_argument,NULL, 0},
	   {NULL,    0,                 NULL, 0  }
	};

	while ((c = getopt_long(argc, argv, short_opt, long_opt, &option_index)) != -1)
	{
		switch (c)
		{
		case -1:       /* no more arguments */
		case 0:
			if (long_opt[option_index].flag != 0)
				break;
			if (strcmp(long_opt[option_index].name, "seed") == 0) { *seed = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "wROG") == 0) { __compute_rog = 1; }
			if (strcmp(long_opt[option_index].name, "wHist") == 0) { write_hist = 1; }
			if (strcmp(long_opt[option_index].name, "wAval") == 0) { write_avalanches = 1; }
			if (strcmp(long_opt[option_index].name, "wImage") == 0) { write_image = 1; }
			if (strcmp(long_opt[option_index].name, "wHull") == 0) { write_hull = 1; }
			if (strcmp(long_opt[option_index].name, "maxLn") == 0) { *max_l = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "minLn") == 0) { *min_l = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "L") == 0) { *L = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "graph") == 0) { __graph_type__ = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "specimen") == 0) { __sample_n__ = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "D") == 0) { *D = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "Ln") == 0) { *Ln = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "BCs") == 0) { *bcs = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "write_times") == 0) { print_write_times(); }
			if (strcmp(long_opt[option_index].name, "verbose") == 0) { __verbose = 1; }
			if (strcmp(long_opt[option_index].name, "version") == 0) { PRINT_VERSION return 0; }
			if (strcmp(long_opt[option_index].name, "help") == 0) {

				printf("__________                             .__    .__                \n");
				printf("\\______   \\____________    ____   ____ |  |__ |__| ____    ____   \n");
				printf(" |    |  _/\\_  __ \\__  \\  /    \\_/ ___\\|  |  \\|  |/    \\  / ___\\    \n");
				printf(" |    |   \\ |  | \\// __ \\|   |  \\  \\___|   Y  \\  |   |  \\/ \\/_\\/>   \n");
				printf(" |______  / |__|  (____  /___|  /\\___  >___|  /__|___|  /\\___  /     \n");
				printf("        \\/             \\/     \\/     \\/     \\/        \\//_____/      \n");
				printf("           __      __.__                                         \n");
				printf("          /  \\    /  \\__| ____   ____   ___________                 \n");
				printf("         \\   \\/\\/   /  |/ __ \\ /    \\_/ __ \\_  __ \\               \n");
				printf("           \\        /|  \\  ___/|   |  \\  ___/|  | \\/             \n");
				printf("            \\__/\\  / |__|\\___  >___|  /\\___  >__|              \n");
				printf("                 \\/          \\/     \\/     \\/          \n");
				printf("\n");

				printf("Simulation of a Branching Wiener Sausage at Critical Point \n");
				printf("Usage: %s [OPTIONS]\n", argv[0]);
				printf("  -h                    set hopping rate\n");
				printf("  -C                    choose the numer of chunks to used. Chunksize will be N/C\n");
				printf("  -N                    choose the numer of realisations to do\n");
				printf("  -L                    choose a specific system size. This takes precedence over other system size parameters\n");
				printf("  -D                    choose the dimension e.g. 1-5. Default is 2\n");

				printf("\n");
				printf("  --seed                choose a sensible seed\n");
				printf("  --Ln                  choose a system size index for L=2^(Ln)-1. This takes precedence over minLn and maxLn\n");
				printf("  --minLn               choose a lower bound system size index for minLn=2^(minLn)-1\n");
				printf("  --maxLn               choose an upper bound system size index for maxLn=2^(maxLn)-1\n");
				printf("  --BCs                 set bitmask to choose specific dimension boundaries. Defaults to 0 i.e. all open boundaries\n");
				//printf("  --chunks              choose chunk size for flush stats (default 1000) \n");

				printf("\n");
				printf("  --wROG                set this flag (without argument) to write radius of gyration data\n");
				printf("  --wHist               set this flag (without argument) to write trace histrogram data\n");
				printf("  --wAval               set this flag (without argument) to write avalanche time integral data\n");
				printf("  --wImage              set this flag (without argument) to write lattice data. WARNING: This can be large!\n");
				printf("  --wHull               set this flag (without argument) to write hull.\n");
				printf("  --verbose             set this flag (without argument) to write additional info e.g (max trace per sample path,..).\n");
				printf("  --specimen            specify sample path N for extended data such as the hull and image. This N value can be determined from sample runs when running with--verbose. Look for #MAX_TRACE(,,index). \n");
				printf("  --graph               choose a particular type of graph. Default to regular lattice or choose sierpinski (1), small-world network(2).\n");

				printf("\n");
				printf("  --write_times         print out the write times used by the program here and now...\n");
				printf("  --version             print out current (git-)version of code\n");
				printf("  --help                print help and exit\n");
				printf("\n");
				return 0;
			}
			break;

		case 'h':
			*h = atof(optarg);
			break;

		case 'L':
			*L = atoi(optarg);
			break;

		case 'N':
			*N = atoi(optarg);
			break;

		case 'C':
			*C = atoi(optarg);
			break;

		case 'D':
			*D = atoi(optarg);
			break;

		case ':':
		case '?':
			fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
			return(-2);

		default:
			fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
			fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
			return(-2);
	};
}
#endif

	//if Ln is specified, use that as the log range - note, if L also specified, L will take precedence in main function
	if (*Ln > 0 && *L < 0) {
		*min_l = *Ln;
		*max_l = *Ln;
	}
	if (*L > 0)MAX_L = *L;
	else MAX_L = pow(2, *max_l) - 1;

	return TRUE;
}

void printhelp() {
	printf("Simulation of a Branching Wiener Sausage at Critical Point \n\
	Correct syntax (UNIX-only): ./bws -S (Seed) -N (#Realisations) -L (L <= 2^input) -D (dim) -B (Boundary-conditions as bitmask) -h [help]\n ");

	//printf("#AVALCOMMENT below are avalanche moment <s^n>. However, read <s^n> = <s^0>^n*<s^n>! Moments have to be inflated by N^n.\n");
}

void print_write_times() {
	INIT_WRITE_TIMES
		int i = 0;
	for (i = 0;i < BINS; i++) {
		printf("%0.1Lf\n", write_times[i]);
	}
}
