/*
 * Created:  23/08/2010 19:39:26
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef KDTREE_H
#define KDTREE_H

#include <stdbool.h>
#include "parseinput.h"

typedef struct kdnode kdnode;
typedef struct kdtree kdtree;

struct kdnode {
    int bucket;
    int cutdim;
    double cutval;
    kdnode *loson, *hison;
    int lopt, hipt;

    bool empty;
    kdnode *father;
    double bnds_x[2], bnds_y[2];
    bool bnds_available;
};

struct kdtree {
    kdnode *root;
    kdnode **bucket_ptr;
    int *perm;
    int bucket_ptr_size;
};

bool ball_in_bounds(data *, double *, double *, int, int);
void kdtree_init(kdtree *, data *);
void kdtree_delete(kdtree *, int);
void kdtree_undelete(kdtree *, int);
void kdtree_undelete_all(kdtree *, int);
void kdtree_free(kdtree *);
void kdtree_dump(kdtree *, data *d);
int kdtree_nearest_neighbour(kdtree *, data *, int);
void kdtree_nearest_neighbour_tour(kdtree *, data *, int, int *, int *);
void kdtree_m_nearest_neighbours(kdtree *, data *, int, void *);
void kdtree_semigreedy_tour(kdtree *, data *, int, int *, int *);
int kdtree_2opt(kdtree *, data *, int *, int *, int);

#endif /* KDTREE_H */
