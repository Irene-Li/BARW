//#basic bosonic L^D lattice
#include "lattice_core.h"
#include "rng/mt19937ar_clean_bkp.h"
#include <stdlib.h>
#include <stdio.h>

int _D = 0, _L = 0, _volume = 0, _BCS = 0, _center = 0;
int lattice_type = 0;

int cache[CACHE_CAPACITY];
char* lattice;
int lattice_actions[20];
int dimension_increment_maps[20];
int wrap_increment_maps[20];

#define TRACE_FLAG_LAT (1<<2) //redefined from bws.c

#define FALSE 0
#define TRUE 1
#define BOUNDARY_FLAG(d) (1<<d)
#define IS_BOUNDARY_CLOSED(d,flag) (BOUNDARY_FLAG(d) == (flag & BOUNDARY_FLAG(d)))
#define IS_ALLOWED_TRANSITION(i,inext,dim_scale) (PARTITION_OF(i, dim_scale*_L) == PARTITION_OF(inext, dim_scale*_L))
#define PARTITION_OF(i,j) (int)( floor(i/((double)j)))
#define IS_ACCESSIBLE(p) (HOLE != (lattice[p] & HOLE))


int Xc(int i, int d) {
	int _d = 0, x = 0, fact = 0;
	for (_d = _D - 1; _d >= 0;_d--) {
		fact = pow(_L, _d);
		x = (int)floor(i / fact);
		if (_d == d)return x;
		i = i % fact;
	}
return(-1);
}


//this code is very literal and crude at the moment and should be optimized and done more elegantly
int get_site_default(int i) {
	if (lattice_type == 0)return 0;
	if (lattice_type == 1) {//sierpinski b=3, m=1
		int b = 3;//move out
		int c = 0, p = 0, p_dash, _p = 0;
		int partitions = log(_L) / log(b);//move out

		for (_p = 0; _p < partitions; _p++)
		{
			int invalid = 1;
			p = pow(b, _p + 1);
			p_dash = pow(b, _p);
			for (c = 0; c < _D; c++) {
				if (((Xc(i, c) % p) / p_dash) != 1) {
					invalid = 0;
					break;
				}
			}
			if (invalid == 1) {//we could not prove it is NOT a hole in D
				//printf("%d,%d\n", Xc(i, 0), Xc(i, 1));
				return 0 | HOLE;
			}
		}
		return 0;
	}
	return 0;
}

#define CACHE_EXCEEDED (1)
#define CACHE_DISABLED (2)

int cache_flags = 0; // (pruess 18 Jul 2021) CACHE_DISABLED to disable cache, otherwise set = 0 
//int capacity_exceeded = -1; // to disable sttucture
int cache_size = 0;
void log_trace(int i) {
	//hint_stack.push(i)
	if (cache_flags & ( CACHE_DISABLED | CACHE_EXCEEDED)) return;//disabled
	if ((cache_size) >= CACHE_CAPACITY) {
		//capacity_exceeded = 1;
		cache_flags |= CACHE_EXCEEDED;
		return;
	}
	cache[cache_size] = i;
	cache_size += 1;
}


/*

int capacity_exceeded = -1; // to disable sttucture
int cache_size = 0;
void log_trace(int i) {
        //hint_stack.push(i)
        if (capacity_exceeded == -1) return;//disabled
        if ((cache_size) == CACHE_CAPACITY) {
                capacity_exceeded = 1;
                return;
        }
        cache[cache_size] = i;
        cache_size += 1;
}

*/





int count_full_resets = 0;
int count_cache_resets = 0;

void __init(int L, int D) {
	_L = L; _D = D;
	_volume = (int)pow(_L, D);

	int bisectAt = round(L / 2);
	int center = 0, i = 0;
	for (i = 0; i < _D; i++) {
		center += (bisectAt * (int)(pow(L, i)));
	}
	_center = center;
}

int allocate_lattice(char **buffer, int L, int D, int type) {
	lattice_type = type;

	__init(L, D);

	int choice;
	for (choice = 0; choice < 2 * D; choice++) {
		int dim = (int)floor(choice / 2);
		int sign = 1; 
		if (choice % 2 == 0) sign = -1;
		//these could be computed once and added to a list
		int j = (int)pow(L, dim); //the jump for moving around the lattice structure e.g. 1, L, L^2 etc. see hier. pend.
		dimension_increment_maps[choice] = j;
		lattice_actions[choice] = sign * j;
		wrap_increment_maps[choice] = (-1 * sign) * (j*L - j);
	}

	char* myBuf = (char*)*buffer;
	int size = _volume * sizeof(char), i = 0;
	lattice = myBuf = (char*)malloc(size);
	for (i = 0;i < _volume;i++) { lattice[i] = get_site_default(i); }//first time init, we may not always run this

	//int test1 = Xc(12, 0);
	//int test2 = Xc(12, 1);
	//int nacc = IS_ACCESSIBLE(30);
	//int yacc = IS_ACCESSIBLE(29);

	if (*buffer == NULL) { printf("\nERROR: Memory allocation did not complete successfully! Required size %i bytes", size); exit(EXIT_FAILURE); }
	return 0;
}

/* pruess 18 July 2021
 * I think the previous version contained the following bug.
 * By default capacity_exceeded=0. 
 * Bug: log_trace would be triggered irrespective of whether]
 * or not it actually needs triggering, because of a typo in
 * the ADD(p) macro.
 * As a result, capacity would often be exceeded and capacity_exceeded=1
 * would almost always be triggered. 
 * Bug: In init_lattice the lattice would be fully flushed only when
 * capacity_exceeded!=1, i.e. when capacity_exceeded=0 or when capacity_exceeded=-1.
 * When capacity was exceeded, i.e. capacity_exceeded=1, then the cache would be
 * used.
 * Further, there was a bug in log_trace, that the cache would be exceeded,
 * but log_trace would carry on writing.
 *
 * Further, when used in the mode of increasing system sizes, not
 * run_for_realisations, then all the information for the boundary condition
 * is not recalculated and therefore the next_pos returned by diffusion 
 * may often be wrong. I turned this off.
 *
 * Fixes: ADD(p) needs curly brackets, init_lattice needs fix, log_trace needs fix, 
 *   init_lattice needs either get_site_default or 0 in both cases, switched off
 *   non-run_for_realisations
 */

void init_lattice(int boundary_conditions, int L, int D)
{
	_BCS = boundary_conditions;
	__init(L, D);
	int i = 0;
	//reset logic follows
	//
	

	if (cache_flags & (CACHE_DISABLED | CACHE_EXCEEDED) ) {//turn off this behaviour using cap_exceeded == -1
	//printf("# Info: Init via volume %i.\n", _volume);
		for (i = 0;i < _volume;i++) { lattice[i] = get_site_default(i); }
		count_full_resets++;
		cache_flags &= ~CACHE_EXCEEDED;
		cache_size=0;

	}
	else {
	//printf("# Info: Init via cache.\n");
		for (i = 0;i < cache_size;i++) {
		/* pruess 18 July 2021: get_site_default(i) is what should happen. In any case, it needs to be the same as in cache_flags & (CACHE_DISABLED | CACHE_EXCEEDED) */
			lattice[cache[i]] = get_site_default(cache[i]); 
			//lattice[cache[i]] = 0;//get_site_default(i); <- if it was visited, it should be zero on reset
		}
		//#warning "Debug only."
		//for (i = 0;i < _volume;i++) { if (lattice[i] != get_site_default(i)) { printf("Site %i has %i should be %i\n", i, lattice[i], get_site_default(i)); exit(EXIT_FAILURE);} }
		count_cache_resets++;
		cache_size = 0;
	}
	//mark bit for boundary - todo
	//check BC for closed boundaries
}

/*
void init_lattice(int boundary_conditions, int L, int D)
{
	_BCS = boundary_conditions;
	__init(L, D);
	int i = 0;
	//reset logic follows
	if (capacity_exceeded == 1) {//turn off this behaviour using cap_exceeded == -1
		for (i = 0;i < cache_size;i++) {
			lattice[cache[i]] = 0;//get_site_default(i); <- if it was visited, it should be zero on reset
		}
		count_cache_resets++;
		cache_size = capacity_exceeded = 0;
	}
	else {
		for (i = 0;i < _volume;i++) { lattice[i] = get_site_default(i); }
		count_full_resets++;
	}
	//mark bit for boundary - todo
	//check BC for closed boundaries
}
*/

int diffuse(int pos, int choice) {
	int pos_next = pos + lattice_actions[choice];
	int d_scale = dimension_increment_maps[choice];
	if (IS_ALLOWED_TRANSITION(pos, pos_next, d_scale)==0)
		return -1;
	if (IS_BOUNDARY_CLOSED(choice, _BCS) == TRUE) {
		//return pos;//reflecting
		pos_next= pos + wrap_increment_maps[choice];//wrap?
	}
 	if (IS_ACCESSIBLE(pos_next)==0) {  return -2;  } //reflect and surface the issue, process can go again
	//if we are wrapping?
	#warning "Unnecessary check?"
	if ((pos_next<0) || (pos_next>=_volume)) { fprintf(stderr, "Illegal pos_next=%i not in [0,%i-1]. Parameters were pos=%i choice=%i lattice_actions[choice]=%i wrap_increment_maps[choice]=%i\n", pos_next, _volume, pos, choice, lattice_actions[choice], wrap_increment_maps[choice]); exit(EXIT_FAILURE); }
	return pos_next;
}

int persist_diffuse2d(int pos, int past_pos, double p, double q, double r) {
	// p1: probabiltiy of going forward 
	// p2: probability of going either of the sideway direction 
	// r: random number between 0 and 1 
	int next_pos; 
	int past_step = pos - past_pos; 
	int step_size; 
	if (r <= p) {
		// Go forward；
		next_pos = pos + past_step;
		step_size = (int) abs(past_step); 

	}
	else if (r > (q*2+p)) {
		// Go backward;
		next_pos = past_pos; 
		step_size = (int) abs(past_step); 
	}
	else {
		step_size = get_sideway_step(pos, past_pos);
		if (r <= p+q) { //  p < r < p+q 
			// Go one of the side direction 
			next_pos = pos + step_size; 
		}
		else { // p+q < r < p+2*q 
			next_pos = pos - step_size; 
		}
	}
	// Check if pos_next is accessible 
	if (IS_ALLOWED_TRANSITION(pos,next_pos,step_size)==0){
		return -1;
	}
	else {
		return next_pos; 
	}

}

tuple branch_2d(int pos, int past_pos, double p, double q, double r) {
	// branch according to a random number
	// p: probabiltiy of going forward 
	// q: probability of going either of the sideway direction 
	tuple return_vals; 
	tuple step_sizes; 
	int past_step = pos - past_pos; 
	int past_step_size = (int) abs(past_step); 
	int sideway_step = get_sideway_step(pos, past_pos);
	if (r < p/(2*p+q)) {
		// branch forward and one side direction 
		return_vals.a = pos + past_step;  
		return_vals.b = pos + sideway_step;
		step_sizes.a = past_step_size; 
		step_sizes.b = sideway_step; 
	}
	else if (r < 2*p/(2*p+q)) {
		// branch forward and the other side direction 
		return_vals.a = pos + past_step; 
		return_vals.b = pos - sideway_step;
		step_sizes.a = past_step_size;
		step_sizes.b = sideway_step; 
	}
	else { 
		// branch sideways 
		return_vals.a = pos + sideway_step; 
		return_vals.b = pos - sideway_step;
		step_sizes.a = sideway_step; 
		step_sizes.b = sideway_step;  
	}

	if (IS_ALLOWED_TRANSITION(pos,return_vals.a,step_sizes.a)==0){ // if not allowed transition
		return_vals.a = -1; 
	}
	if (IS_ALLOWED_TRANSITION(pos,return_vals.b,step_sizes.b)==0){ // if not allowed transition
		return_vals.b = -1; 
	}

	return return_vals;
}

int get_sideway_step(int pos, int past_pos) {
	int past_step = pos - past_pos;
	int step_size; 

	if (abs(past_step) == 1) {
		step_size = _L; // set it equal to the new step size 
	}
	else {
		step_size = 1; 
	} 
	return step_size; 
}

int get_center(int L, int D) {
	__init(L, D);

	if (lattice_type == 1) {
		//int random_dir = abs(genrand_int32() % (2 * _D));
		int random_dir = genrand_int32() % (2 * _D);
		int random_dir_offset = lattice_actions[random_dir];
		int sierpinski_tile_offset = random_dir_offset * (int)ceil(L / (double)(3 * 2));//see notes S3,8 always has central tile of side sqrt(L). So we place ourparticle beside it
		//printf("center at %d\n", sierpinski_tile_offset);
		return _center + sierpinski_tile_offset; //must be a random offset from center - not a single one
	}

	return _center;
}

int probe(int p, int orientation, int offsets[]) {
	//consider lattice boundary cond.
	int L = offsets[4];
	if (orientation != -1) {
		//first check that the offset is legal for this position - if not return 0
		int temp = p;
		int offset = offsets[orientation];
		temp += offset;
		if (IS_ALLOWED_TRANSITION(p, temp, L) == FALSE) {//2d only
			return 0;
		}
		p = temp;
	}
	if (TRACE_FLAG_LAT == (lattice[p] & TRACE_FLAG_LAT))
		return 1;
	return 0;
}


void right_turning_walk(int L, int n, long double t) {
	int pos = L / 2, orientation = 0;
	int offsets[8] = { -L,-L + 1,1,L + 1,+L,+L - 1,-1,-L - 1 };
	while (probe(pos, -1, offsets) == FALSE) {
		pos += L;
	} //move down
	int ring_start = pos;
	printf("#HULL(%.3Lf,%d): ", t, n);
	do {
		//turn until we find something
		while (probe(pos, orientation, offsets) == FALSE) {
			orientation = (orientation + 1) % 8;
		}
		printf("%d ", pos);
		//move the cursor to this something
		pos += offsets[orientation];
		//re-orient - by looking "backwards" to where you came from + 1 to the right. There should be nothing there
		orientation = (orientation + 5) % 8;
		//terminate when we come back around
	} while (ring_start != pos);
	printf("\n");
}

double distance_from_center(int i) {
	int _d = 0;
	long double sdist = 0;
	for (_d = 0; _d < _D; _d++) {
		sdist += pow(Xc(_center, _d) - Xc(i, _d), 2);
	}
	return sdist;//sqrt(dist);
}
