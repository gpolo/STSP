/*
 * Created:  23/08/2010 09:51:39
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include "util.h"
#include "kdtree.h"
#include "heapmax.h"
#include "parseinput.h"

struct rcl {
    int rcl_size, rcl_used;
    struct heapelm *candidates;
};

static void
rmnear(kdtree *tree, data *d, kdnode *p, int target, int *nndist,
        struct rcl *pr)
{
    int i;
    int thisdist;
    double diff;

    if (p->empty)
        return;

    if (p->bucket) {
        for (i = p->lopt; i <= p->hipt; i++) {
            if (tree->perm[i] != target) {
                thisdist = d->dist(d, tree->perm[i], target);
                if (thisdist < *nndist) {
                    /*
                    printf("N(%d): %d [%d]\n", target, tree->perm[i], thisdist);
                    */
                    *nndist = thisdist;

                    /* Possibly add a new element to the RCL. */
                    if (pr->rcl_used < pr->rcl_size) {
                        pr->candidates[pr->rcl_used].num = tree->perm[i];
                        pr->candidates[pr->rcl_used].cost = thisdist;
                        pr->rcl_used++;
                        if (pr->rcl_used == pr->rcl_size) {
                            heap_build_max(pr->candidates, pr->rcl_size);
                        }
                    } else {
                        heap_replace_root(pr->candidates, pr->rcl_size,
                                tree->perm[i], thisdist);
                    }

                }
            }
        }
    } else {
        assert(p->cutdim == 0 || p->cutdim == 1);

        diff = ((p->cutdim == 0) ? d->x[target] : d->y[target]) - p->cutval;
        if (diff < 0) {
            rmnear(tree, d, p->loson, target, nndist, pr);
            if (*nndist >= -diff) {
                rmnear(tree, d, p->hison, target, nndist, pr);
            }
        } else {
            rmnear(tree, d, p->hison, target, nndist, pr);
            if (*nndist >= diff) {
                rmnear(tree, d, p->loson, target, nndist, pr);
            }
        }
    }
}

void
kdtree_m_nearest_neighbours(kdtree *tree, data *d, int target, void *pr)
{
    double thisval;
    int nndist;
    kdnode *p, *lastp;

    nndist = INT_MAX;

    p = tree->bucket_ptr[target];
    rmnear(tree, d, p, target, &nndist, pr);
    while (1) {
        lastp = p;
        p = p->father;
        if (p == NULL)
            break;
       
        assert(p->cutdim == 0 || p->cutdim == 1);
        thisval = (p->cutdim == 0) ? d->x[target] : d->y[target];

        if (lastp == p->loson) {
            if (thisval + nndist > p->cutval) {
                rmnear(tree, d, p->hison, target, &nndist, pr);
            }
        }
        else {
            if (thisval - nndist < p->cutval) {
                rmnear(tree, d, p->loson, target, &nndist, pr);
            }
        }
    }
}


#define BIAS_LOG(r) (1/log((r)+1))
#define BIAS_LINEAR(r) (1./(r))
#define BIAS_EXP(r) (1/exp(r))
#define BIAS_RANDOM(r) (1)
#define BIAS(r) BIAS_LOG(r)

void
kdtree_semigreedy_tour(kdtree *tree, data *d, int start, int *tour, int *len)
{
    int dim = d->numnodes;
    int i, current, next;
    int lenght = 0;
    /* RCL. */
    int j, k;
    double total_bias, prob_acc, prob_pick, prob;
    struct rcl params;

    params.candidates = MALLOC(dim, struct heapelm);

    tour[0] = start;
    current = start;
    for (i = 1; i < dim; i++) {
        kdtree_delete(tree, current);

        /* Define RCL size. */
        params.rcl_used = 0;
        params.rcl_size = RANDOM(3) + 1;

        kdtree_m_nearest_neighbours(tree, d, current, &params);

        /*printf("RCL size: %d, used: %d\n", params.rcl_size, params.rcl_used);
        heap_dump(params.candidates, params.rcl_used);*/

        /* Pick an element from RCL. */
        total_bias = 0;
        for (j = 1; j <= params.rcl_used; j++) {
            total_bias += BIAS(j);
        }

        prob_acc = 0;
        prob_pick = RANDOM_UNIT();
        for (j = params.rcl_used, k = 1; j >= 0; j--, k++) {
            prob = BIAS(k) / total_bias;
            if (prob + prob_acc > prob_pick) {
                break;
            }
            prob_acc += prob;
        }
        /*printf("Picked %d: %d, %d\n", j, params.candidates[j - 1].num,
                params.candidates[j - 1].cost);*/
        next = params.candidates[j - 1].num;

        /* Update tour. */
        tour[i] = next;
        lenght += params.candidates[j - 1].cost;
        current = next;
    }
    lenght += d->dist(d, current, start);
    *len = lenght;

    kdtree_undelete_all(tree, d->numnodes);

    free(params.candidates);
}
