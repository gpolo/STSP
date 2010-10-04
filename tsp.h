/*
 * Created:  09/09/2010 14:48:33
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef TSP_H
#define TSP_H

#include <stdio.h>

#define DEBUGGING 0
#define DEBUG(str, vaargs...) do { \
    if (DEBUGGING) fprintf(stderr, str, ##vaargs); \
} while(0)

/* Flags used for command line parameters. */
#define SET_SEED      1
#define USE_KDTREE    2
#define USE_MATRIX    4
#define DO_2HHOPT     8
#define DO_GRASP      16
#define DO_3OPT       32
#define ONLY_GREEDYNN 64
#define BE_VERBOSE    128

#endif /* TSP_H */
