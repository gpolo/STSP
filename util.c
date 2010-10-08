/*
 * Created:  23/08/2010 18:16:40
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "util.h"

#define M 131072
static int cachesig[M];
static int cacheval[M];

int chit, cmiss;

/* Took from Bentley's paper.
 * I have yet to find a case where this actually decreases the running time.
 * XXX Not used. */
int
cached_dist(data *d, int i, int j)
{
    int temp, ind;

    if (i < j) {
        temp = i; i = j; j = temp;
    }

    ind = i ^ j;
    if (cachesig[ind] != i) {
        cachesig[ind] = i;
        cacheval[ind] = d->dist(d, i, j);
        /*cmiss++;
    } else {
        chit++; */
    }
    return cacheval[ind];
}


void *
Malloc(size_t size)
{
    void *ptr = malloc(size);
    if (!ptr)
        abort();
    return ptr;
}

void
show_sol(int *tour, int length, int cost)
{
    int i;

    for (i = 0; i < length; i++) {
        printf("%d ", tour[i]);
    }
    printf(": %d\n", cost);
}

int
tour_cost(data *d, int *tour, int dim)
{
    int i, len = 0;

    for (i = 1; i < dim; i++) {
        len += d->dist(d, tour[i - 1], tour[i]);
    }
    len += d->dist(d, tour[dim - 1], tour[0]);

    return len;
}

int
check_tour_validity(int *tour, int dim)
{
    int i, status = 0;
    int *visited = MALLOC(dim, int);

    memset(visited, 0, sizeof(int) * dim);

    for (i = 0; i < dim; i++) {
        if (visited[tour[i]]) {
            fprintf(stderr, "Not a hamiltonian cycle, city %d was "
                            "visited at %d.\n", tour[i], visited[tour[i]]);
            status = -1;
            break;
        }
        visited[tour[i]] = i + 1;
    }
    for (i = 0; i < dim; i++) {
        if (!visited[i]) {
            fprintf(stderr, "Not a hamiltonian cycle, city %d was not "
                            "visited.\n", i);
            status = -1;
            break;
        }
    }

    free(visited);
    return status;
}

double
current_usertime_secs(void)
{
    double usertime, systemtime;
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) < 0) {
        perror("getrusage");
        return -1;
    }
    usertime = usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec * 1e-6);
    systemtime = usage.ru_stime.tv_sec + (usage.ru_stime.tv_usec * 1e-6);
    return (usertime + systemtime);
}

int
seed_prng(long int seed, int flag, long int *genseed)
{
    int gseed;
    struct timeval tval;

    if (flag & SEED_NEW) {
        if (gettimeofday(&tval, NULL) < 0) {
            perror("gettimeofday");
            return 1;
        }
        gseed = tval.tv_sec * tval.tv_usec;
    } else {
        gseed = seed;
    }

    srand48(gseed);
    if (genseed != NULL)
        *genseed = gseed;

    return 0;
}
