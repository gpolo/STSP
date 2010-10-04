/*
 * Created:  05/09/2010 00:47:39
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <assert.h>
#include <stdio.h>
#include "util.h"
#include "dist.h"
#include "kdtree.h"

int
main(void)
{
    data d;
    kdtree tree;

    d.numnodes = 6;
    d.x = MALLOC(d.numnodes, double);
    d.y = MALLOC(d.numnodes, double);

    d.x[0] = 5; d.x[1] = 25; d.x[2] = 20; d.x[3] = 10; d.x[4] = 30; d.x[5] = 10;
    d.y[0] = 5; d.y[1] = 10; d.y[2] = 15; d.y[3] = 20; d.y[4] = 25; d.y[5] = 30;

    kdtree_init(&tree, &d);
    kdtree_dump(&tree, &d);

    d.dist = dist_euc2d;
    assert(kdtree_nearest_neighbour(&tree, &d, 5) == 3);

    free(d.x);
    free(d.y);

    return 0;
}
