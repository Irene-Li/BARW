#pragma once

#include <math.h>

#define HOLE (1<<6)
#define CACHE_CAPACITY 1000000

int inline get_site_default(int i);
void log_trace(int i);
double distance_from_center(int i);
int diffuse(int pos, int choice);
void init_lattice(int boundary_conditions, int L, int D);
int allocate_lattice(char** buffer, int L, int D, int type);
int get_center(int L, int D);
int probe(int loc, int dir, int offsets[]);
void right_turning_walk(int L, int n, long double time);

