/*
 * Created:  23/08/2010 18:17:19
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include "parseinput.h"

#define SWAP(a, b, type) do {type temp = a; a = b; b = temp;} while(0)
#define MALLOC(count, type) Malloc(count * sizeof(type))
#define RANDOM(limit) (lrand48() % (limit))
#define RANDOM_UNIT() (drand48())
#define SEED_REUSE 0
#define SEED_NEW 1

int chit, cmiss;

int cached_dist(data *d, int, int);
void *Malloc(size_t);
int seed_prng(long int, int, long int *);
void show_sol(int *, int, int);
int check_tour_validity(int *, int);
int tour_cost(data *, int *, int);
double current_usertime_secs(void);

#endif /* UTIL_H */
