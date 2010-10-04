/*
 * Created:  24/08/2010 16:53:13
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "parseinput.h"
#include "heapmax.h"
#include "util.h"
#include "naive.h"
#include "tsp.h"

static void
semi_greedy(data *d, int start, int *tour, int *tour_length)
{
    int i, j, k;
    int dim = d->numnodes;
    int current, len, dist;
    int *visited = MALLOC(dim, int);
    /* RCL. */
    struct heapelm *candidates = MALLOC(dim, struct heapelm); 
    int rclsize, rclused;
    double total_bias, rprob, racc, rpick;

    memset(visited, 0, dim * sizeof(int));

    len = 0;
    current = start;
    for (i = 0; i < dim - 1; i++) {
        visited[current] = 1;
        tour[i] = current;

        /* Define RCL size. */
        rclused = 0;
        if (RANDOM_UNIT() < 0.05) rclsize = dim - i - 2;
        else rclsize = RANDOM(((dim - i - 2) / 4) + 1);
        rclsize++;

        DEBUG("RCL size: %d\n", rclsize);

        /* Define RCL. */
        for (j = 0; j < dim; j++) {
            if (!visited[j]) {
                dist = d->dist(d, current, j);
                if (rclused < rclsize) {
                    candidates[rclused].num = j;
                    candidates[rclused].cost = dist;
                    rclused++;

                    if (rclused == rclsize) {
                        heap_build_max(candidates, rclsize);
                    }
                } else if (dist < candidates[0].cost) {
                    heap_replace_root(candidates, rclsize, j, dist);
                }
            }
        }

        DEBUG("RCL used: %d\n", rclused);
        heap_dump(candidates, rclused);

        /* Pick RCL element based on logarithmic bias. */
#define BIAS_LOG(r) (1/log((r)+1))
#define BIAS_LINEAR(r) (1./(r))
#define BIAS_EXP(r) (1/exp(r))
#define BIAS_RANDOM(r) 1
#define BIAS(r) BIAS_LOG(r)
        total_bias = 0;
        for (j = 1; j <= rclused; j++) {
            total_bias += BIAS(j);
        }
        racc = 0;
        rpick = RANDOM_UNIT();
#if DEBUGGING
        for (j = rclused; j >= 1; j--) {
            DEBUG("%f ", BIAS(j) / total_bias);
        }
        DEBUG("\nR: %f\n", rpick);
#endif
        for (j = rclused, k = 1; j >= 0; j--, k++) {
            rprob = BIAS(k) / total_bias;
            if (rprob + racc > rpick) {
                break;
            }
            racc += rprob;
        }
        DEBUG("Picked: j = %d, %d %d\n\n", j,
                candidates[j-1].num, candidates[j-1].cost);


        current = candidates[j - 1].num;
        len += candidates[j - 1].cost;
    }
    tour[i] = current;
    len += d->dist(d, current, start);
    *tour_length = len;

    free(visited);
    free(candidates);

#undef BIAS
#undef BIAS_LOG
}


int
naive_grasp(data *d, int flags)
{
    int dim = d->numnodes;
    int start, cost, swaps;
    int *tour, *best_tour, best_cost, no_improvements;
    int use_three_opt = flags & DO_3OPT;
    long int i, limit;
    bool break_line = false;
    
    tour = MALLOC(dim, int);
    best_tour = MALLOC(dim, int);

    limit = (long int)dim * dim;
    printf("\nGRASP, %ld runs starting.. \n", limit);
    best_cost = INT_MAX;
    no_improvements = 0;
    for (i = 0; i < limit; i++) {
        start = RANDOM(dim);
        semi_greedy(d, start, tour, &cost);
        if (use_three_opt) { 
            local_3opt(d, tour, &cost, &swaps);
        } else {
            local_2opt(d, tour, &cost, &swaps);
        }
        if (cost < best_cost) {
            if (flags & BE_VERBOSE) {
                if (break_line) {
                    break_line = false;
                    printf("\n");
                }
                printf("%ld: %d, %d\n", i, cost, start);
            }
            best_cost = cost;
            memcpy(best_tour, tour, dim * sizeof(int));
            no_improvements = 0;
        } else {
            no_improvements++;
            if (flags & BE_VERBOSE && !(no_improvements % 50)) {
                printf(">");
                break_line = true;
                fflush(stdout);
            }
        }
        if (no_improvements == 5000) {
            if (flags & BE_VERBOSE) {
                if (break_line) {
                    printf("\n");
                    break_line = false;
                }
                printf("No improvements for %d runs, stopping.\n",
                        no_improvements);
            }
            break;
        }
    }
    if (flags & BE_VERBOSE) {
        if (break_line) printf("\n");
        show_sol(best_tour, dim, best_cost);
        printf("Executions: %ld\n", i);
    }
    printf("\nGRASP final cost: %d\n", best_cost);

    free(tour);
    free(best_tour);
    return 0;
}
