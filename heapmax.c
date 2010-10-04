/*
 * Created:  24/08/2010 16:53:13
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include "util.h"
#include "tsp.h"
#include "heapmax.h"

#define LEFT(i) (2*i)
#define RIGHT(i) (2*i + 1)

static void
max_heapify(struct heapelm *heap, int size, int i)
{
    struct heapelm temp;
    int l = LEFT(i), r = RIGHT(i);
    int largest;

    if (l < size && heap[l].cost > heap[i].cost) {
        largest = l;
    } else {
        largest = i;
    }
    if (r < size && heap[r].cost > heap[largest].cost) {
        largest = r;
    }

    if (largest != i) {
        temp = heap[i];
        heap[i] = heap[largest];
        heap[largest] = temp;

        max_heapify(heap, size, largest);
    }
}


void
heap_build_max(struct heapelm *heap, int size)
{
    int i;

    for (i = (size / 2); i >= 0; i--) {
        max_heapify(heap, size, i);
    }
}

#if DEBUGGING
void
heap_dump(struct heapelm *heap, int size)
{
    int i;

    for (i = 0; i < size; i++) {
        DEBUG("%d:%d, ", heap[i].num, heap[i].cost);
    }
    DEBUG("\n");
}
#else
void
heap_dump(struct heapelm *heap, int size)
{
    return;
}
#endif

void
heap_replace_root(struct heapelm *heap, int size, int newnum, int newcost)
{
    heap[0].cost = newcost;
    heap[0].num = newnum;
    max_heapify(heap, size, 0);
}
