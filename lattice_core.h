#ifndef LATTICE_CORE_H
#define LATTICE_CORE_H
#pragma once

#include <math.h>

#define HOLE (1<<6)
#define CACHE_CAPACITY 100000
//#warning "Debug only"
//#define CACHE_CAPACITY 10

int inline get_site_default(int i);
void log_trace(int i);
double distance_from_center(int i);
int diffuse(int pos, int choice);
int persist_diffuse2d(int pos, int past_pos, double p, double q, double r); 
void init_lattice(int boundary_conditions, int L, int D);
int allocate_lattice(char** buffer, int L, int D, int type);
int get_center(int L, int D);
int probe(int loc, int dir, int offsets[]);
void right_turning_walk(int L, int n, long double time);

#ifndef LATTICE_CORE_C
extern int count_cache_resets;
extern int count_full_resets;
#endif

#endif
