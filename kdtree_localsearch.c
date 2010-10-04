/*
 * Created:  23/08/2010 09:51:39
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */

/* See "BENTLEY, J. L.: Fast Algorithms for Geometric Traveling Salesman
 * Problems. 1992" */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "parseinput.h"
#include "kdtree.h"
#include "util.h"


struct opt2 {
    int a, b, c;
    int gain;
    int tour_size;
    int *tour, *index;

    /* For multiple considerations. */
    int curr_b, curr_c;

    /* 2Q-opt */
    bool allow_2qopt;
    int c_prev;
    bool opt_2q;
};

static int
attempt_2opt_or_2qopt(data *dat, int c, struct opt2 *o)
{
    int a = o->a;
    int b = o->b;
    int d = o->tour[(o->index[c] + 1) % o->tour_size];
    int c_prev; /* For 2Q-opt. */
    int gain;

    if (d == a)
        return 0;

    gain = dat->dist(dat, a, c) + dat->dist(dat, b, d);
    gain -= dat->dist(dat, a, b) + dat->dist(dat, c, d);

    if (gain < o->gain) {
        if (o->index[a] > o->index[c]) {
            o->curr_b = d;
            o->curr_c = a;
        } else { 
            o->curr_b = b;
            o->curr_c = c;
        }
        o->gain = gain;
        o->opt_2q = false;
        return 1;
    } else {
        if (!o->allow_2qopt)
            return 0;

        /* 2Q-opt. */
        if (o->index[a] > o->index[c]) {
            a = c; b = d;
            c = o->a; d = o->b;
        }

        c_prev = o->tour[(o->index[c] - 1) % o->tour_size];
        if (c_prev == b)
            return 0;

        gain = dat->dist(dat,a,c)+dat->dist(dat,c,b)+dat->dist(dat,c_prev,d);
        gain -= dat->dist(dat,a,b)+dat->dist(dat,c_prev,c)+dat->dist(dat,c,d);

        if (gain < o->gain) {
            o->curr_b = b;
            o->c_prev = c_prev;
            o->curr_c = c;
            o->opt_2q = true;
            o->gain = gain;
            return 1;
        }
    }

    return 0;
}

int nconsidered = 0;
int ncavg = 0;

static int
rfrnn(kdtree *tree, kdnode *p, data *d, int rad, double x, double y,
        struct opt2 *o)
{
    int i, c;
    int target = o->a;
    int thisdist;
    double diff;

    if (p->empty)
        return 0;

    if (p->bucket) {
        for (i = p->lopt; i <= p->hipt; i++) {
            c = tree->perm[i];
            if (c == target)
                continue;
            thisdist = d->dist(d, c, target);
            if (thisdist < rad) {
                nconsidered++;
                if (attempt_2opt_or_2qopt(d, c, o)) {
                    if (nconsidered > 0) {
                        /* XXX Didn't gain much by considering more
                         * neighbors (actually found worst results).
                         * Keeping it at > 0 for now. */
                        return 1;
                    } else { 
                        /* Continue searching. */
                        return 0;
                    }
                }
            }
        }
        return 0;
    } else {
        assert(p->cutdim == 0 || p->cutdim == 1);

        diff = ((p->cutdim == 0) ? x : y) - p->cutval;
        if (diff < 0) {
            if (rfrnn(tree, p->loson, d, rad, x, y, o))
                return 1;
            if (rad >= -diff) {
                if (rfrnn(tree, p->hison, d, rad, x, y, o))
                    return 1;
            }
        } else {
            if (rfrnn(tree, p->hison, d, rad, x, y, o))
                return 1;
            if (rad >= diff) {
                if (rfrnn(tree, p->loson, d, rad, x, y, o))
                    return 1;
            }
        }
    }

    return 0;
}

int
kdtree_fixed_radius_nearneighbour(kdtree *tree, data *d, int rad,
        struct opt2 *o)
{
    int target = o->a;
    double diff;
    double x = d->x[target];
    double y = d->y[target];
    kdnode *lastp, *p = tree->bucket_ptr[target];

    if (rfrnn(tree, p, d, rad, x, y, o))
        return 1;

    while (1) {
        lastp = p;
        p = p->father;
        if (p == NULL)
            return 0;

        assert(p->cutdim == 0 || p->cutdim == 1);

        diff = ((p->cutdim == 0) ? x : y) - p->cutval;
        if (lastp == p->loson) {
            if (rad >= -diff) {
                if (rfrnn(tree, p->hison, d, rad, x, y, o))
                    return 1;
            }
        } else {
            if (rad >= diff) {
                if (rfrnn(tree, p->loson, d, rad, x, y, o))
                    return 1;
            }
        }
        if (p->bnds_available &&
                ball_in_bounds(d, p->bnds_x, p->bnds_y, target, rad)) {
            return 0;
        }
    }
}


static inline void
invert(int *path, int *index, int start, int end)
{
    if (start > end) {
        SWAP(start, end, int);
    }
    for (; start <= end; start++) {
        SWAP(index[path[start]], index[path[end]], int);
        SWAP(path[start], path[end], int);
        end--;
    }
}

static inline void
slide(int *path, int *index, struct opt2 *data)
{
    int i;
    int aux1 = data->curr_c, aux2 = index[data->curr_b];

    for (i = index[data->c_prev]; i >= index[data->curr_b]; i--) {
        index[path[i]] = i + 1;
        path[i + 1] = path[i];
    }
    path[aux2] = aux1;
    index[aux1] = aux2;
}

int
kdtree_2opt(kdtree *tree, data *d, int *tour, int *len, int do_2qopt)
{
    int i, dim = d->numnodes;
    int *neighbor, *tour_index;
    int old_swaps, total_gain = 0, swaps = 0;
    int a, b;
    int rad;
    struct opt2 opt2_data;

    double start_time, inversion_time = 0;

    neighbor = MALLOC(dim, int);
    tour_index = MALLOC(dim, int);

    for (i = 0; i < dim; i++) {
        neighbor[i] = kdtree_nearest_neighbour(tree, d, i);
    }

    for (i = 0; i < dim; i++) {
        tour_index[tour[i]] = i;
    }

    opt2_data.tour_size = dim;
    opt2_data.index = tour_index;
    opt2_data.tour = tour;
    opt2_data.allow_2qopt = do_2qopt;

    do {
        old_swaps = swaps;

        for (i = 0; i < dim; i++) {
            a = tour[i];
            b = tour[(i + 1) % dim];

            if (neighbor[a] == b) {
                continue;
            }

            rad = d->dist(d, a, b);
            assert(rad > 0);
            opt2_data.a = a;
            opt2_data.b = b;
            opt2_data.curr_b = b;
            opt2_data.gain = 0;
            opt2_data.opt_2q = false;

            nconsidered = 0;

            kdtree_fixed_radius_nearneighbour(tree, d, rad, &opt2_data);
            if (opt2_data.gain < 0) {
                ncavg += nconsidered;
                swaps++;
                total_gain += opt2_data.gain;

                start_time = current_usertime_secs();
                if (!opt2_data.opt_2q) {
                    invert(tour, tour_index,
                        tour_index[opt2_data.curr_b],
                        tour_index[opt2_data.curr_c]);
                } else { /* 2Q-opt swap. */
                    slide(tour, tour_index, &opt2_data);
                }
                inversion_time += (current_usertime_secs() - start_time);

                /* Decrement i by 2 so the loop increments it by 1
                 * and we effectively back up one step in the loop.
                 * (Suggestion from the article at top) */
                if (i == 0) i++; /* If 'i' is 0, don't let it go negative. */
                i -= 2;
            }
        }
    } while (old_swaps != swaps);

    /*
    printf("Time used for inversions: %f s\n", inversion_time);
    printf("> %d swaps, gain : %d\n", swaps, -total_gain);
    printf("Average neighbours considered: %f\n", (double)ncavg / swaps);*/
    ncavg = 0;

    *len += total_gain;

    free(neighbor);
    free(tour_index);

    return swaps;
}
