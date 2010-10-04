/*
 * Created:  23/08/2010 09:51:39
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */

/* See "BENTLEY, J. L.: K-d Trees for Semidynamic Point Sets. 1990" */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include "kdtree.h"
#include "parseinput.h"
#include "util.h"

/* The paper mentioned above mentions the following good values for
 * the two following params. */
#define BUCKET_CUTOFF 5 /* > 0. */
#define BNDSLEVEL 3

#define MIN_SORTSIZE 16


static int
findmaxspread(kdtree *tree, int l, int u, double *x, double *y)
{
    int i;
    double xval, xmin, xmax, xspread;
    double yval, ymin, ymax, yspread;

    xmax = xmin = x[tree->perm[l]];
    ymax = ymin = y[tree->perm[l]];

    for (i = l + 1; i <= u; i++) {
        xval = x[tree->perm[i]];
        yval = y[tree->perm[i]];

        if (xval < xmin) xmin = xval;
        else if (xval > xmax) xmax = xval;

        if (yval < ymin) ymin = yval;
        else if (yval > ymax) ymax = yval;
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;
    return ((xspread >= yspread) ? 0 : 1);
}

static void
insertion_sort(int *array, double *val, int n)
{
    int i, j;
    int temp;

    for (i = 1; i < n; i++) {
        temp = array[i];
        for (j = i; j > 0 && val[array[j - 1]] > val[temp]; j--) {
            array[j] = array[j - 1];
        }
        array[j] = temp;
    }
}

/* Robert Sedgewick, Algorithms in C, Chapter 9 - Quicksort (1990).
 * This is a merge of the "partition" and "select" algorithms
 * plus an adaptation to use insertion sort when the size of the
 * subarray gets small enough. */
static void
select_rs(int *array, double *val, int l, int r, int k)
{
    double v;
    int i, j;

    while (r > l && r - l > MIN_SORTSIZE) {
        v = val[array[r]];
        i = l - 1;
        j = r;
        while (1) {
            do i++; while (val[array[i]] < v && i < r);
            do j--; while (val[array[j]] > v && j > l);
            if (i >= j)
                break;
            SWAP(array[i], array[j], int);
        }
        SWAP(array[i], array[r], int);
        if (i >= k) r = i - 1;
        if (i <= k) l = i + 1;
    }

    insertion_sort(array + l, val, r - l + 1);
}

static kdnode *
build(kdtree *tree, int depth, double bx0, double bx1, double by0, double by1,
        int l, int u, double *x, double *y)
{
    int m, i;
    kdnode *p = MALLOC(1, kdnode);

    p->empty = false;
    depth++;

    if (u-l+1 <= BUCKET_CUTOFF) {
        p->bnds_available = false;
        p->bucket = 1;
        p->lopt = l;
        p->hipt = u;
        for (i = l; i <= u; i++) {
            tree->bucket_ptr[tree->perm[i]] = p;
        }
    } else {
        if (depth % BNDSLEVEL) {
            p->bnds_available = false;
        } else {
            p->bnds_available = true;
            p->bnds_x[0] = bx0; p->bnds_x[1] = bx1;
            p->bnds_y[0] = by0; p->bnds_y[1] = by1;
        }
        p->bucket = 0;
        p->cutdim = findmaxspread(tree, l, u, x, y);
        m = (l + u) / 2;
        assert(p->cutdim == 0 || p->cutdim == 1);
        if (p->cutdim == 0) {
            select_rs(tree->perm, x, l, u, m);
            p->cutval = x[tree->perm[m]];

            p->loson = build(tree, depth, bx0, p->cutval, by0, by1,
                    l, m, x, y);
            p->hison = build(tree, depth, p->cutval, bx0, by0, by1,
                    m + 1, u, x, y);
        } else {
            select_rs(tree->perm, y, l, u, m);
            p->cutval = y[tree->perm[m]];

            p->loson = build(tree, depth, bx0, bx1, by0, p->cutval,
                    l, m, x, y);
            p->hison = build(tree, depth, bx0, bx1, p->cutval, by1,
                    m + 1, u, x, y);
        }

        p->loson->father = p;
        p->hison->father = p;
    }
    return p;
}

void
kdtree_init(kdtree *tree, data *d)
{
    int i;
    int numnodes = d->numnodes;

    tree->perm = MALLOC(numnodes, int);
    for (i = 0; i < numnodes; i++)
        tree->perm[i] = i;

    tree->bucket_ptr_size = numnodes;
    tree->bucket_ptr = MALLOC(numnodes, kdnode *);
    tree->root = build(tree, 0, INFINITY, -INFINITY, INFINITY, -INFINITY,
            0, numnodes - 1, d->x, d->y);
    tree->root->father = NULL;
}

bool
ball_in_bounds(data *d, double bx[], double by[], int target, int dist)
{
    bool in_x = d->x[target] >= bx[0] && d->x[target] <= bx[1];
    bool in_y = d->y[target] >= by[0] && d->y[target] <= by[1];

    if (in_x && in_y) {
        if (dist > d->x[target] - bx[0]) return false;
        if (dist > bx[1] - d->x[target]) return false;
        if (dist > d->y[target] - by[0]) return false;
        if (dist > by[1] - d->y[target]) return false;
        return true;
    }
    return false;
}

/* kdtree_delete actually "hides" the point. */
void
kdtree_delete(kdtree *tree, int pointnum)
{
    int j;
    kdnode *p;

    p = tree->bucket_ptr[pointnum];
    if (p->empty) {
        /* Object already removed (the entire bucket is empty).
         * (Something else may be wrong). */
        fprintf(stderr, "ERROR: Point %d has been removed earlier\n",
                pointnum);
        return;
    }
    
    j = p->lopt;
    while (tree->perm[j] != pointnum) {
        j++;
    }
    SWAP(tree->perm[j], tree->perm[p->hipt], int);
    (p->hipt)--;
    if (p->lopt > p->hipt) {
        p->empty = true;
        while ((p = p->father) != NULL && p->loson->empty && p->hison->empty)
            p->empty = true;
    }
}

/* kdtree_undelete "unhides" the point. */
void
kdtree_undelete(kdtree *tree, int pointnum)
{
    int j;
    kdnode *p;

    p = tree->bucket_ptr[pointnum];
    j = p->lopt;
    while (tree->perm[j] != pointnum) {
        j++;
    }

    if (j > p->hipt) {
        (p->hipt)++;
        SWAP(tree->perm[j], tree->perm[p->hipt], int);
        p->empty = false;
        while ((p = p->father) != NULL && p->empty)
            p->empty = false;
    }
}

void
kdtree_undelete_all(kdtree *tree, int numnodes)
{
    int i;

    for (i = 0; i < numnodes; i++) {
        kdtree_undelete(tree, i);
    }
}


void free_tree(kdnode *node)
{
    if (node->bucket) {
        free(node);
        return;
    }
    free_tree(node->loson);
    free_tree(node->hison);
    free(node);
}

void
kdtree_free(kdtree *tree)
{
    free(tree->perm);
    free(tree->bucket_ptr);
    free_tree(tree->root);
    return;
}


static void
dump(kdtree *tree, kdnode *node, data *d, int depth)
{
#define INDENT do { int i; for (i = 0; i < depth; i++) putchar(' '); } while(0)

    int i;

    if (node == NULL || node->empty)
        return;

    INDENT;
    if (node->bucket == 0) {
        assert(node->cutdim == 0 || node->cutdim == 1);
        printf("Split %c at %f\n", node->cutdim ? 'Y' : 'X', node->cutval);
        dump(tree, node->loson, d, depth + 2);
        dump(tree, node->hison, d, depth + 2);
    } else {
        printf("Bucket, ");
        for (i = node->lopt; i <= node->hipt; i++) {
            printf("(%.1f;%.1f) ", d->x[tree->perm[i]], d->y[tree->perm[i]]);
        }
        printf("\n");
    }
}

void
kdtree_dump(kdtree *tree, data *d)
{
    int i;

    if (!tree) {
        return;
    }
    printf("Data:\n  X: ");
    for (i = 0; i < d->numnodes; i++) {
        printf("%.2f ", d->x[i]);
    }
    printf("\n  Y: ");
    for (i = 0; i < d->numnodes; i++) {
        printf("%.2f ", d->y[i]);
    }

    printf("\nData permutation: ");
    for (i = 0; i < tree->bucket_ptr_size; i++) {
        printf("%d ", tree->perm[i]);
    }
    printf("\n\n");

    dump(tree, tree->root, d, 0);
}
