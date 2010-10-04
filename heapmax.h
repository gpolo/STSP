/*
 * Created:  10/09/2010 02:33:26
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef HEAPMAX_H
#define HEAPMAX_H

struct heapelm {
    int cost;
    int num;
};

void heap_build_max(struct heapelm *, int);
void heap_dump(struct heapelm *, int);
void heap_replace_root(struct heapelm *, int, int, int);

#endif /* HEAPMAX_H */
