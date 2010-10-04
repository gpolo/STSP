/*
 * Created:  24/08/2010 16:52:34
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef NAIVE_H
#define NAIVE_H

#include "parseinput.h"

void local_2opt(data *, int *, int *, int *);
void local_3opt(data *, int *, int *, int *);
int naive(data *, int);

#endif /* NAIVE_H */
