/*
 * Created:  23/08/2010 09:51:39
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include "kdtree.h"
#include "parseinput.h"


static void
rnn(kdtree *tree, data *d, kdnode *p, int target, int *nndist, int *nearest)
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
                    *nndist = thisdist;
                    *nearest = tree->perm[i];
                }
            }
        }
    } else {
        assert(p->cutdim == 0 || p->cutdim == 1);

        diff = ((p->cutdim == 0) ? d->x[target] : d->y[target]) - p->cutval;
        if (diff < 0) {
            rnn(tree, d, p->loson, target, nndist, nearest);
            if (*nndist >= -diff) {
                rnn(tree, d, p->hison, target, nndist, nearest);
            }
        } else {
            rnn(tree, d, p->hison, target, nndist, nearest);
            if (*nndist >= diff) {
                rnn(tree, d, p->loson, target, nndist, nearest);
            }
        }
    }
}

int
kdtree_nearest_neighbour(kdtree *tree, data *d, int target)
{
    double thisval;
    int nndist, nearest;
    kdnode *p, *lastp;

    nndist = INT_MAX;

    p = tree->bucket_ptr[target];
    rnn(tree, d, p, target, &nndist, &nearest);
    while (1) {
        lastp = p;
        p = p->father;
        if (p == NULL)
            break;
       
        assert(p->cutdim == 0 || p->cutdim == 1);
        thisval = (p->cutdim == 0) ? d->x[target] : d->y[target];

        if (lastp == p->loson) {
            if (thisval + nndist > p->cutval) {
                rnn(tree, d, p->hison, target, &nndist, &nearest);
            }
        }
        else {
            if (thisval - nndist < p->cutval) {
                rnn(tree, d, p->loson, target, &nndist, &nearest);
            }
        }
        /* XXX ball_in_bounds not used for nn. */
    }
    return nearest;
}

void /* XXX retornar status (tour == NULL, outros) */
kdtree_nearest_neighbour_tour(kdtree *tree, data *d, int start,
        int *tour, int *len)
{
    int i, current, next;
    int lenght = 0;

    tour[0] = start;
    current = start;
    for (i = 1; i < d->numnodes; i++) {
        kdtree_delete(tree, current);
        next = kdtree_nearest_neighbour(tree, d, current);
        tour[i] = next;
        lenght += d->dist(d, current, next);
        current = next;
    }
    lenght += d->dist(d, current, start);
    *len = lenght;

    kdtree_undelete_all(tree, d->numnodes);
}
