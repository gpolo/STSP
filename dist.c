/*
 * Created:  23/08/2010 18:22:05
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <math.h>
#include "dist.h"

/* QC stands for Quick & Cheap (take care). */
#define QCROUND(x) ((x) >= 0 ? (int)((x) + 0.5) : (int)((x) - 0.5))
#define QCCEIL(x) (((x) - (int)(x)) == 0 ? (int)(x) : \
        (int)(x + ((x) >= 0 ? 1 : -1)))


int
matrix_at(data *unused, int i, int j)
{
    return small_matrix[i][j];
}

int
dist_ceil2d(data *d, int i, int j)
{
    double xd, yd;
    double result;

    xd = d->x[i] - d->x[j];
    yd = d->y[i] - d->y[j];
    result = sqrt(xd * xd + yd * yd);
    return QCCEIL(result);
}

int
dist_euc2d(data *d, int i, int j)
{
    double xd, yd;
    double result;

    xd = d->x[i] - d->x[j];
    yd = d->y[i] - d->y[j];
    result = sqrt(xd * xd + yd * yd);
    return QCROUND(result);
}

int
dist_att(data *d, int i, int j)
{
    double xd, yd;
    double rij;
    int tij;

    xd = d->x[i] - d->x[j];
    yd = d->y[i] - d->y[j];

    rij = sqrt((xd * xd + yd * yd) / 10);
    tij = QCROUND(rij);
    return ((tij < rij) ? (tij + 1) : tij);
}

int
dist_geo(data *d, int i, int j)
{
    double q1, q2, q3;
    double lat1, lon1, lat2, lon2;

    lat1 = d->x[i]; lon1 = d->y[i];
    lat2 = d->x[j]; lon2 = d->y[j];

    q1 = cos(lon1 - lon2);
    q2 = cos(lat1 - lat2);
    q3 = cos(lat1 + lat2);

    return (EARTH_RADIUS * acos(0.5 * ((1 + q1)*q2 - (1 - q1)*q3)) + 1);
}

void
geo_latlon(double *lat, double *lon)
{
    int deg;

    deg = (int)*lat;
    *lat = M_PI * (deg + 5 * (*lat - deg) / 3) / 180;
    deg = (int)*lon;
    *lon = M_PI * (deg + 5 * (*lon - deg) / 3) / 180;
}

