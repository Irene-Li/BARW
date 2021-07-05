/* 
 * $Header: /home/ma/p/pruess/.cvsroot/branching_wiener_sausage/branching_wiener_reference.c,v 1.9 2017/04/06 12:38:50 pruess Exp $
 *
 * This is a reference implementation (supposedly) of the branching 
 * Wiener moments. 
 */
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>



/*
 cc -Wall -O3 -o branching_wiener_reference branching_wiener_reference.c mt19937ar_clean_bkp.c
 ./branching_wiener_reference_run1D.sh
./branching_wiener_reference_process.sh branching_wiener_reference_L*1D*.dat > branching_wiener_reference_1D_0.5.dat

p 'branching_wiener_reference_1D_0.05.dat' u 2:9:10 w e, x**0.5
p 'branching_wiener_reference_1D_0.05.dat' u 2:11:12 w e, 1e2*x
p 'branching_wiener_reference_1D_0.5.dat' u 2:11:12 w e, 1e2*x

./branching_wiener_reference_process_tt.sh branching_wiener_reference_L5111D_0.05.dat

p 'branching_wiener_reference_1D_0.05.dat' u 2:($9/$2**0.25) , 'branching_wiener_reference_1D_0.5.dat' u 2:($9/$2**0.25)

creep, creep, creep. It does look like the second moment is dieing off.

I wonder whether the clash with the literature is down to system sizes. I will send some stuff off on the cluster.
 */



#define DIMENSION (2)

#define xxMAX_TIME  (0.2)

#define LENGTH (200)
int length[DIMENSION]={0};

#define EPSILON (0.05)
#define SIGMA (0.05)
#define HOPPING (0.51)

double epsilon=EPSILON;
double sigma=SIGMA;
double hopping=HOPPING;

#undef EPSILON
#define EPSILON epsilon
#undef SIGMA
#define SIGMA sigma
#undef HOPPING
#define HOPPING hopping

int num_chunks=-1;
#define ITERATIONS (1000LL)

long long int iterations=ITERATIONS;
#undef ITERATIONS
#define ITERATIONS (iterations)


#define TT_NUM_SLOTS (200)

long long int total_malloced;
#define MALLOC(a,n) {if ((a=malloc((n)*sizeof(*(a))))==NULL) \
   {fprintf(stderr, "# Fatal error (%s::%i) malloc(3)ing %i bytes for variable %s: %i (%s).\n", \
         __FILE__, __LINE__, (int)(n*sizeof(*(a))), #a, errno, strerror(errno)); exit(EXIT_FAILURE);} else {\
	    total_malloced+=((long long int) ((n)*sizeof((*a)) )); printf("#Info: malloc(3)ing %i items of size %i for %s, %lli bytes, total %lli\n", (int)(n), (int)sizeof(*(a)), #a, ((long long int) ((n)*sizeof((*a)) )), total_malloced); } }

/*******************************************/
/* Moment facilities */
// SNIP: moments
#define MOMENT_TYPE long double   /*@\colabel{moment_type}@*/
#define MOMENT_OUT_FMT "%10.20Lg"

MOMENT_TYPE mom_pow;
int         mom_cnt;
int         mom_tt;
int 	    repeated_up_to=-1;

#define MOMENTS_DECL(n,m) MOMENT_TYPE mom_ ##n[(m) + 1]; const int mom_max_ ##n = m;
#define MOMENTS_INIT(n) {for (mom_cnt=0; mom_cnt<=mom_max_ ##n; mom_cnt++) mom_ ##n [mom_cnt]=0;}
#define MOMENTS(n,x) {mom_ ##n [0]++; for (mom_pow=1, mom_cnt=1; mom_cnt<=mom_max_ ##n; mom_cnt++) {mom_pow*=((MOMENT_TYPE)(x)); mom_ ##n[mom_cnt]+=mom_pow;}}
#define MOMENTS_MAGIC "#M_"
#define MOMENTS_OUT(n)  {printf("%s" MOMENTS_MAGIC "%s %i %i %g %g %g %i %lli %lli %i " MOMENT_OUT_FMT " ", (repeated_up_to>=chunk) ? "#REPEATED" : "", \
     #n, chunk, length[0], hopping, epsilon, sigma, DIMENSION, num_sites, (long long int)ITERATIONS, mom_max_ ##n, mom_ ##n [0]); \
  for (mom_cnt=0; mom_cnt<=mom_max_ ##n; mom_cnt++) \
  printf(" " MOMENT_OUT_FMT, mom_ ##n [mom_cnt]/((mom_ ##n [0]) ? mom_ ##n [0] : -1)); fputc('\n', stdout);}

#define MOMENTS_DECL_TT(n,m,l) MOMENT_TYPE **mom_ ##n; const int mom_max_ ##n = m; const int mom_tt_num=l;
#define MOMENTS_MALLOC_TT(n) {MALLOC(mom_ ##n, mom_tt_num); for (mom_tt=0; mom_tt<mom_tt_num; mom_tt++) {MALLOC(mom_ ##n[mom_tt], mom_max_ ##n+1);}}
#define MOMENTS_INIT_TT(n) {for (mom_tt=0; mom_tt<mom_tt_num; mom_tt++) {for (mom_cnt=0; mom_cnt<=mom_max_ ##n; mom_cnt++) mom_ ##n [mom_tt][mom_cnt]=0;}}
#define MOMENTS_TT(n,x,l) {mom_ ##n [l][0]++; for (mom_pow=1, mom_cnt=1; mom_cnt<=mom_max_ ##n; mom_cnt++) {mom_pow*=((MOMENT_TYPE)(x)); mom_ ##n[l][mom_cnt]+=mom_pow;}}
#define MOMENTS_OUT_TT(n)  {for (mom_tt=0; (tt_slot_time[mom_tt]>=0.); mom_tt++) {if (mom_ ##n [mom_tt][0]) {printf("%s" MOMENTS_MAGIC "%s %i %i %g %i %g %g %g %i %lli %lli %i " MOMENT_OUT_FMT " ", (repeated_up_to>=chunk) ? "#REPEATED" : "", \
     #n, chunk, mom_tt, tt_slot_time[mom_tt], length[0], hopping, epsilon, sigma, DIMENSION, num_sites, (long long int)ITERATIONS, mom_max_ ##n, mom_ ##n [mom_tt][0]); \
  for (mom_cnt=0; mom_cnt<=mom_max_ ##n; mom_cnt++) \
  printf(" " MOMENT_OUT_FMT, mom_ ##n [mom_tt][mom_cnt]/((mom_ ##n [mom_tt][0]) ? mom_ ##n [mom_tt][0] : -1)); fputc('\n', stdout);}}}

MOMENTS_DECL(area, 8); /*@\hfill@*//* Degeneracy of the maximum. */
MOMENTS_DECL_TT(areatt, 8, TT_NUM_SLOTS);



/*******************************************/
/********************************************/
/*Gunnar's boolean version of MT starts here*/
/* RNG (lifted from opgp_cp.c) */
#include "mt19937ar_clean_bkp.h"
#ifndef SEED
#define SEED (5UL)
#endif
#define RNG_MT_BITS (32)
#define RNG_MAX (~0UL)
#define RNG_TYPE unsigned long
#define RNG_MT_MAX (4294967295.)
#define RNG_TYPE_OUT "%lu"

RNG_TYPE seed = SEED;
void init_genrand(RNG_TYPE);
RNG_TYPE genrand_int32(void);
double genrand_real1(void);
double genrand_real2(void);


RNG_TYPE mt_bool_rand=0UL;
int      mt_bool_ptr=RNG_MT_BITS;
RNG_TYPE mt_bool_mask=1UL<<(RNG_MT_BITS-1);
RNG_TYPE mt_8_rand=0UL;
int      mt_8_shifted=0;
RNG_TYPE mt_8_mask=7UL<<(RNG_MT_BITS-2-3);
RNG_TYPE mt_8_randV2=0UL;
int      mt_8_shiftedV2=RNG_MT_BITS-2;
RNG_TYPE mt_4_rand=0UL;
int      mt_4_shifted=0;
RNG_TYPE mt_4_mask=3UL<<(RNG_MT_BITS-2);
RNG_TYPE mt_4_randV2=0UL;
int      mt_4_shiftedV2=RNG_MT_BITS;
RNG_TYPE mt_2_randV2=0UL;
int      mt_2_shiftedV2=RNG_MT_BITS;


#define RANDOM_INTEGER() genrand_int32()
#define OLD_RNG_MT_BOOLEAN ( ( mt_bool_ptr==RNG_MT_BITS ) ?  \
        ((mt_bool_ptr=1, mt_bool_rand=genrand_int32()) & 1) : \
        (mt_bool_rand & (1<<mt_bool_ptr++)) )
#define RNG_MT_BOOLEAN ( ( mt_bool_mask==(1UL<<(RNG_MT_BITS-1)) ) ?  \
        ((mt_bool_mask=1UL, mt_bool_rand=genrand_int32()) & mt_bool_mask) : \
        (mt_bool_rand & (mt_bool_mask+=mt_bool_mask)) )
#define RNG_MT_4 ( ( mt_4_mask==(3UL<<(RNG_MT_BITS-2)) ) ?  \
        ((mt_4_shifted=0, mt_4_mask=3UL, mt_4_rand=genrand_int32()) & mt_4_mask) : \
        (mt_4_shifted+=2, (mt_4_rand & (mt_4_mask<<=2))>>mt_4_shifted) )
/* 30 bits are available */
#define RNG_MT_8 ( ( mt_8_mask==(7UL<<(RNG_MT_BITS-2-3)) ) ?  \
        ((mt_8_shifted=0, mt_8_mask=7UL, mt_8_rand=genrand_int32()) & mt_8_mask) : \
        (mt_8_shifted+=3, (mt_8_rand & (mt_8_mask<<=3))>>mt_8_shifted) )
#define RNG_MT_8V2 ( ( mt_8_shiftedV2==RNG_MT_BITS-2 ) ?  \
        ((mt_8_shiftedV2=3, mt_8_randV2=genrand_int32()) & 7UL) : \
        (mt_8_shiftedV2+=3, (mt_8_randV2>>=3) & 7UL) )
#define RNG_MT_4V2 ( ( mt_4_shiftedV2==RNG_MT_BITS ) ?  \
        ((mt_4_shiftedV2=2, mt_4_randV2=genrand_int32()) & 3UL) : \
        (mt_4_shiftedV2+=2, (mt_4_randV2>>=2) & 3UL) )
#define RNG_MT_2V2 ( ( mt_2_shiftedV2==RNG_MT_BITS ) ?  \
        ((mt_2_shiftedV2=1, mt_2_randV2=genrand_int32()) & 1UL) : \
        (mt_2_shiftedV2++, (mt_2_randV2>>=1) & 1UL) )

/* /RNG */
/*Gunnar's boolean version of MT ends here*/
/********************************************/



int *stack_particle;
int stack_particle_height, stack_particle_height_max, max_stack_particle_height;
int *stack_area;
int stack_area_height, stack_area_height_max, max_stack_area_height;

#define POP_X(a,s)  a=stack_ ## s[--stack_ ## s ## _height]
/* max_stack_ assignment used to be at the beginning, before stack_height was updated. */
#define PUSH_X(a,s) if (stack_## s ##_height<stack_## s ## _height_max) {stack_ ##s[stack_ ##s## _height++]=a; if (stack_## s ##_height>max_stack_## s ## _height) {max_stack_## s ## _height=stack_## s ##_height;} } else {fprintf(stderr, "# Error: stack %s: stack_height %i exceeds maximum %i in line %i.\n", #s, stack_ ##s## _height, stack_ ##s## _height_max, __LINE__); printf("# Error: stack %s: stack_height %i exceeds maximum %i in line %i.\n", #s, stack_ ##s## _height, stack_ ##s## _height_max, __LINE__); exit(EXIT_FAILURE); }




#define NEW(a) PUSH_X(a,particle); MARK(a);
#define DELETE(s) stack_particle[s]=stack_particle[--stack_particle_height]
#define MARK(a) if (lattice[a]==0) {PUSH_X((a),area); lattice[a]=1;}



char *lattice;


#define INDEX2COO(p,index) { int tmp_i; long long int pos; \
     pos=index;\
     p[0]=pos % length[0]; \
     for (tmp_i=1; tmp_i<DIMENSION; tmp_i++) {\
       pos/=length[tmp_i-1]; \
       p[tmp_i]=pos % length[tmp_i];}}

#define COO2INDEX(pos,p) { int tmp_i; \
     pos=0LL;\
     for (tmp_i=DIMENSION-1; tmp_i>0; tmp_i--) {\
       pos+=p[tmp_i];\
       pos*=length[tmp_i-1];\
     } pos+=p[tmp_i];}




int main(int argc, char *argv[])
{
double t;
long long int num_sites=1;
int chunk;
long long int it;
int start_pos;
int i, j, dim, move;
double process;
double hopping_threshold, branching_threshold;
double total_rate;
int t_slot_num;
double tt_slot_time[TT_NUM_SLOTS]={
  0.,  1.,  2.,  5., 
      10., 20., 50.,
     100.,200.,500.,
     1000.,2000.,5000.,
     10000.,20000.,50000.,
     100000.,200000.,500000.,
     1000000.,2000000.,5000000.,
     -1.};
int spos;
int coo[DIMENSION];
int ch;
#ifdef DEBUG
int ohoh_untouched, ohoh_hop;
double ohoh_m0=0., ohoh_m1=0.0;
#endif

setlinebuf(stdout);

printf("# Info: Command: %s", argv[0]);
for (i=1; i<argc; i++) 
  printf(" %s", argv[i]);
printf("\n");

for (i=0; i<DIMENSION; i++) length[i]=LENGTH;

while ((ch=getopt(argc, argv, "c:i:s:L:R:"))!=-1) {
  switch (ch) {
    case 'c':
      num_chunks=atoi(optarg);
      break;
    case 'i':
      iterations=atoll(optarg);
      break;
    case 's':
      seed=strtoul(optarg, NULL, 10);
      break;
    case 'L':
      for (i=0; i<DIMENSION; i++) 
	if ((length[i]=atoi(optarg))<=0) {
	  fprintf(stderr, "length %i not allowed.\n", length[i]);
          exit(EXIT_FAILURE);
        }
       break;
     case 'R':
       if (sscanf(optarg, "%lg:%lg:%lg", &hopping, &epsilon, &sigma)==3) {
         printf("# Info: Read rates hopping=%g, epsilon=%g, sigma=%g\n", hopping, epsilon, sigma);
       } else {
         fprintf(stderr, "Failed to read rates (%i::%s)\n", errno, strerror(errno));
         exit(EXIT_FAILURE);
       }
     default:
       break;
   }
}

hopping_threshold=hopping/(epsilon+sigma+hopping);
branching_threshold=(sigma+hopping)/(epsilon+sigma+hopping);

for (i=0; i<DIMENSION; i++) {
  num_sites*=length[i];
  coo[i]=(length[i]-1)/2;
}
/* Does not work in d>1:
 * start_pos=(num_sites-1)/2;
 *
 * This was a bug I hunted for a long time.
 */

COO2INDEX(start_pos, coo);




tt_slot_time[0]=0.0;
tt_slot_time[1]=0.01;
for (i=2; i<TT_NUM_SLOTS; i++) {
  tt_slot_time[i]=tt_slot_time[i-1]*1.2;
}
tt_slot_time[TT_NUM_SLOTS-1]=-1.;

/* heuristics... XXX this ought to be a command line parameter */
if (DIMENSION==1) stack_particle_height_max=(int)((double)(num_sites*num_sites*2)/hopping_threshold);
else stack_particle_height_max=num_sites*100; 
stack_area_height_max=num_sites;

  printf("# Info: %s\n", "$Header: /home/ma/p/pruess/.cvsroot/branching_wiener_sausage/branching_wiener_reference.c,v 1.9 2017/04/06 12:38:50 pruess Exp $.");
  #define PRINT_PARAM(a,f) fprintf(stdout, "# Info " #a ": " f "\n", a)
  PRINT_PARAM(num_sites, "%lli");
  PRINT_PARAM(DIMENSION, "%i");
  for (i=0; i<DIMENSION; i++) {
    PRINT_PARAM(i, "%i");
    PRINT_PARAM(length[i], "%i");
  }
  PRINT_PARAM(start_pos, "%i");
  INDEX2COO(coo, start_pos);
  for (i=0; i<DIMENSION; i++) {
    PRINT_PARAM(i, "%i");
    PRINT_PARAM(coo[i], "%i");
  }
  PRINT_PARAM(HOPPING, "%g");
  PRINT_PARAM(EPSILON, "%g");
  PRINT_PARAM(SIGMA, "%g");
  PRINT_PARAM(hopping, "%g");
  PRINT_PARAM(epsilon, "%g");
  PRINT_PARAM(sigma, "%g");
  PRINT_PARAM(SEED, "%lu");
  PRINT_PARAM(seed, "%lu");
  PRINT_PARAM(num_chunks, "%i");
  PRINT_PARAM(iterations, "%lli");
  PRINT_PARAM(hopping_threshold, "%g");
  PRINT_PARAM(branching_threshold, "%g");
  PRINT_PARAM(stack_particle_height_max, "%i");
  PRINT_PARAM(stack_area_height_max, "%i");


  //exit(0);

  i=0;
  do {
    printf("# Info: tt_slot_time[%i]=%g\n", i, tt_slot_time[i]);
    } while (tt_slot_time[i++]>=0);

  MALLOC(lattice, num_sites);

  MALLOC(stack_particle, stack_particle_height_max);
  max_stack_particle_height=0;

  MALLOC(stack_area, stack_area_height_max);
  max_stack_area_height=0;

  MOMENTS_MALLOC_TT(areatt);
  init_genrand(seed);

  for (i=0; i<num_sites; i++) lattice[i]=0;

  for (chunk=1; ((chunk<=num_chunks) || (num_chunks<0)); chunk++) {
    MOMENTS_INIT(area);
    MOMENTS_INIT_TT(areatt);

#ifdef DEBUG
ohoh_m0=0;
ohoh_m1=0;
#endif
    for (it=1LL; it<=ITERATIONS; it++) {
#ifdef DEBUG
      ohoh_untouched=1;
      ohoh_hop=1;
#endif
      t=0;
      t_slot_num=0;
      stack_particle_height=0;
      NEW(start_pos);
#ifdef MAX_TIME
      while ((stack_particle_height) && (stack_area_height<num_sites) && (t<=MAX_TIME)) {
#else
      while ((stack_particle_height) && (stack_area_height<num_sites)) {
#endif
	total_rate=stack_particle_height*(EPSILON+SIGMA+HOPPING);
//printf("Time0 is %g max %g\n", t, MAX_TIME);
	t+=(-log(1.-genrand_real2())/total_rate);
//printf("Time1 is %g max %g\n", t, MAX_TIME);

#ifdef DEBUG
if ((t>0.01) && (ohoh_untouched)) {
  ohoh_m0++;
  //ohoh_m1+=stack_area_height;
  ohoh_m1+=ohoh_hop;
  ohoh_untouched=0;
}
#endif


  /* this should really be t>=, so I catch the very first one. Later I realised that I need to update t before taking measurements of the _old_ state. */
	while ((t>=tt_slot_time[t_slot_num]) && (tt_slot_time[t_slot_num]>=0.)) {
  /*
	  if ((t_slot_num==0) && (stack_area_height>1)) {
          printf("t_slot_num=%i but stack_area_height=%i\n", t_slot_num, stack_area_height);
          exit(0);
        }
*/
#ifdef DEBUG
if (tt_slot_time[t_slot_num]==0.01)
printf("For time %g at time %g I have %i\n", tt_slot_time[t_slot_num], t, stack_area_height);
#endif
        MOMENTS_TT(areatt, stack_area_height, t_slot_num);
	t_slot_num++;
      }
      spos=genrand_int32()%stack_particle_height;
      process=genrand_real2();

      if (process<hopping_threshold) {
#ifdef DEBUG
if (t<0.01)
printf("process: At time %g (%i) versus hopping_threshold: %g at rate %g of %g stack_area_height %i\n", t, t_slot_num, process, hopping_threshold, total_rate, stack_area_height);
ohoh_hop++;
//printf("process: %g versus hopping_threshold: %g at rate %g\n", process, hopping_threshold, total_rate);
#endif
        INDEX2COO(coo,stack_particle[spos]);
#ifdef DEBUG
printf("coo: %i %i (%i)\n", coo[0], coo[1], stack_particle[spos]);
#endif
	move=genrand_int32() % (2*DIMENSION);
	coo[dim=move/2]+=((move & 1) ? 1 : -1);
	if ((coo[dim]<0) || (coo[dim]>=length[dim])) {DELETE(spos); }
	else { 
	  COO2INDEX(stack_particle[spos], coo);
          MARK(stack_particle[spos]);
	}
      } else if (process<branching_threshold) {
//printf("process: %g versus branching_threshold: %g\n", process, branching_threshold);
        NEW(stack_particle[spos]);
      } else {
        DELETE(spos);
      }
/*
 * early on forgotten to delete.
for (i=j=0; i<num_sites; i++) 
  if (lattice[i]) j++;
if (j!=stack_area_height) {
  printf("j is %i but stack_area_height is %i\n", j, stack_area_height);
}
*/
    } /* while stack_particle_height */
  while (tt_slot_time[t_slot_num]>=0.) {
    MOMENTS_TT(areatt, stack_area_height, t_slot_num);
    t_slot_num++;
  }
  MOMENTS(area, stack_area_height);
  while (stack_area_height) lattice[stack_area[--stack_area_height]]=0;
  } /* iterations */

#ifdef DEBUG
printf("# ohoh: %g %g %g\n", ohoh_m0, ohoh_m1, ohoh_m1/ohoh_m0);
#endif
printf("# Info: max_stack_particle_height=%i and max_stack_area_height=%i so far.\n", max_stack_particle_height, max_stack_area_height);
MOMENTS_OUT(area);
MOMENTS_OUT_TT(areatt);
} /* chunk */

    


}
