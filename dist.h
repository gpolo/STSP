/*
 * Created:  23/08/2010 18:21:01
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef DIST_H
#define DIST_H

#include "parseinput.h"

#define EARTH_RADIUS 6378.388

int matrix_at(data *d, int, int);
int dist_ceil2d(data *, int, int);
int dist_euc2d(data *, int, int);
int dist_att(data *, int, int);
int dist_geo(data *, int, int);

void geo_latlon(double *, double *);

#endif /* DIST_H */
