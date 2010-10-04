/*
 * Created:  23/08/2010 19:43:26
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <stdbool.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include "kdtree.h"
#include "naive.h"
#include "naive_grasp.h"
#include "tsp.h"

#include <assert.h> /* XXX Used only in commented code. */
#include <limits.h> /* XXX Used by code that should be somewhere else. */

void
usage(const char *prog)
{
    fprintf(stderr, "Usage: %s [OPTIONS] input.tsp\n", prog);
    fprintf(stderr, "\nOPTIONS:\n");
    fprintf(stderr, "\t-S seed\t  Set a seed\n");
    fprintf(stderr, "\t-G\t  Run GRASP\n");
    fprintf(stderr, "\t-K\t  Use a k-d tree to represent the instance\n");
    fprintf(stderr, "\t-M\t  Use a matrix to represent the instance or find a\n"
                    "\t  \t  solution in naïve mode\n");
    fprintf(stderr, "\t-4\t  Perform 2.(1/4)-Opt in local search when ");
    fprintf(stderr, "using a k-d tree\n");
    fprintf(stderr, "\t-N\t  Run Greedy-NN only\n");
    fprintf(stderr, "\t-3\t  Use 3-opt for local search [only matrix]\n");
    fprintf(stderr, "\t-V\t  Be a bit more verbose\n");
    exit(1);
}

static int
try_openinput(int argc, char *argv[], long int *seed)
{
#ifdef LEAVE
#undef LEAVE
#endif
#define LEAVE(msg) fprintf(stderr, msg); exit(1)

    int ch;
    int status = 0;

    if (argc < 2) {
        usage(argv[0]);
    } else {
        while ((ch = getopt(argc, argv, "s:S:gGkKmM4nN3vV")) != -1) {
            switch (ch) {
            case 's':
            case 'S':
                status |= SET_SEED;
                *seed = atol(optarg);
                break;
            case 'g':
            case 'G':
                if (status & ONLY_GREEDYNN) {
                    LEAVE("CONFLICT: Cannot run Grasp when 'Greed only' "
                            "is specified.\n");
                }
                status |= DO_GRASP;
                break;
            case 'k':
            case 'K':
                if (status & USE_MATRIX) {
                    LEAVE("CONFLICT: Cannot use both matrix and k-d tree.\n");
                } else if (status & DO_3OPT) {
                    LEAVE("NOT SUPPORTED: Cannot run 3-opt on k-d tree yet.\n");
                }
                status |= USE_KDTREE;
                break;
            case 'm':
            case 'M':
                if (status & USE_KDTREE) {
                    LEAVE("CONFLICT: Cannot use both k-d tree and matrix.\n");
                }
                status |= USE_MATRIX;
                break;
            case '4':
                if (status & ONLY_GREEDYNN) {
                    LEAVE("CONFLICT: Cannot run 2.1/4-Opt on Greedy.\n");
                }
                status |= DO_2HHOPT;
                break;
            case 'n':
            case 'N':
                if (status & DO_2HHOPT || status & DO_3OPT ||
                        status & DO_GRASP) {
                    LEAVE("CONFLICT: Local search in use.\n");
                }
                status |= ONLY_GREEDYNN;
                break;
            case '3':
                if (status & ONLY_GREEDYNN) {
                    LEAVE("CONFLICT: Cannot run 3-Opt on Greedy.\n");
                } else if (status & USE_KDTREE) {
                    LEAVE("NOT SUPPORTED: "
                            "3-opt not implemented on k-d tree yet.\n");
                }
                status |= DO_3OPT;
                break;
            case 'v':
            case 'V':
                status |= BE_VERBOSE;
                break;

            case '?':
            default:
                usage(argv[0]);
            }
        }
        if (!(status & USE_MATRIX) && !(status & USE_KDTREE)) {
            LEAVE("MISSING: Need to set the representation structure.\n");
        }

        argc -= optind;
        argv += optind;
        if (argc != 1) {
            usage((argv - optind)[0]);
        } else {
            if (freopen(*argv, "rb", stdin) == NULL) {
                perror("freopen");
                exit(1);
            }
        }
    }

    return status;

#undef LEAVE
}


int
main(int argc, char *argv[])
{
    int argc_info, verbose, start, result = 0;
    long int seed;
    data d;
    kdtree tree;
    double start_time, build_time, read_time, greedy_time, show_time, opt2_time;
    struct Info input_info;

    int *rota, tamanho, trocas;


    argc_info = try_openinput(argc, argv, &seed);
    verbose = argc_info & BE_VERBOSE;

    if (seed_prng(seed, ((argc_info & SET_SEED) ? SEED_REUSE : SEED_NEW),
                &seed) < 0) {
        return 1;
    }
    if (get_info(&input_info) < 0)
        return 1;

    d.numnodes = input_info.dimension;
    d.x = MALLOC(d.numnodes, double);
    d.y = MALLOC(d.numnodes, double);

    printf("Seed: %ld\n", seed);
    printf("Working on instance '%s' using %s...\n", input_info.name,
            (argc_info & USE_KDTREE) ? "k-d tree" : "matrix");
    free(input_info.name);
    start_time = current_usertime_secs();
    if (argc_info & USE_MATRIX) { /* Naïve */
        if (read_explicit(&input_info, &d) < 0) {
            result = 1;
            goto exit;
        }
        if (argc_info & DO_GRASP) {
            result = naive_grasp(&d, argc_info);
        } else {
            result = naive(&d, argc_info);
        }
        goto exit;
    } else if (setdist_readdata(&input_info, &d) < 0) {
        result = 1;
        goto exit;
    }
    read_time = current_usertime_secs() - start_time;

    start_time = current_usertime_secs();
    kdtree_init(&tree, &d);
    build_time = current_usertime_secs() - start_time;

    /*kdtree_dump(&tree, &d);*/

    rota = MALLOC(d.numnodes, int);

    /* XXX Move all the following to somewhere else. */
    if (argc_info & DO_GRASP) {
        bool break_line = false;
        long int i;
        long int limit = (long int)d.numnodes * d.numnodes;
        int best_cost = INT_MAX, start;
        int *best_tour = MALLOC(d.numnodes, int);

        int no_improvements = 0;
        printf("\nGRASP, %ld runs starting..\n", limit);
        for (i = 0; i < limit; i++) {
            start = RANDOM(d.numnodes);
            kdtree_semigreedy_tour(&tree, &d, start, rota, &tamanho);
            kdtree_2opt(&tree, &d, rota, &tamanho, argc_info & DO_2HHOPT);
            if (tamanho < best_cost) {
                if (verbose) {
                    if (break_line) {
                        break_line = false;
                        printf("\n");
                    }
                    printf("%ld: %d, %d\n", i, tamanho, start);
                }
                best_cost = tamanho;
                memcpy(best_tour, rota, d.numnodes * sizeof(int));
                no_improvements = 0;
            } else {
                no_improvements++;
                if (verbose && !(no_improvements % 100)) {
                    printf("><");
                    break_line = true;
                    fflush(stdout);
                }
            }
            if (no_improvements == 5000) {
                if (verbose) {
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
        if (verbose) {
            if (break_line) printf("\n");
            show_sol(best_tour, d.numnodes, best_cost);
            /*assert(tour_cost(&d, best_tour, d.numnodes) == best_cost);
            check_tour_validity(best_tour, d.numnodes);*/
            printf("Executions: %ld\n", i);
        }
        printf("\nGRASP final cost: %d\n", best_cost);
        free(best_tour);
        free(rota);
        result = 0;
        goto exit_kdtree;
    }

    start_time = current_usertime_secs();
    start = RANDOM(d.numnodes);
/*
    int i, aa;
    aa = INT_MAX;
    for (i = 0; i < d.numnodes; i++) {
        kdtree_init(&tree, &d);
        kdtree_nearest_neighbour_tour(&tree, &d, i, rota, &tamanho);
        if (tamanho < aa) {
            aa = tamanho;
            printf("Best at %d: %d\n", i, tamanho);
        }
        kdtree_free(&tree);
    }
    exit(1);
*/
    kdtree_nearest_neighbour_tour(&tree, &d, start,
            rota, &tamanho);
    /*check_tour_validity(rota, d.numnodes);*/
    greedy_time = current_usertime_secs() - start_time;

    start_time = current_usertime_secs();
    if (verbose) {
        show_sol(rota, d.numnodes, tamanho);
    }
    printf("\nStart point: %d\n", start);
    printf("Greedy-NN tour cost: %d\n", tamanho);
    show_time = current_usertime_secs() - start_time;


    if (!(argc_info & ONLY_GREEDYNN)) {
        start_time = current_usertime_secs();
        trocas = kdtree_2opt(&tree, &d, rota, &tamanho, argc_info & DO_2HHOPT);
        opt2_time = current_usertime_secs() - start_time;
        printf("Cost after running 2-opt: %d %s[%d swaps]\n", tamanho,
                (argc_info & DO_2HHOPT) ? "(with 2.1/4-Opt) " : "", trocas);
    }
    /*assert(tour_cost(&d, rota, d.numnodes) == tamanho);*/
    if (verbose) {
        printf("\nStructure build time: %f s\n", build_time);
        printf("Greedy tour time: %f s\n", greedy_time);
        printf("Read input time: %f s\n", read_time);
        printf("Display solution time: %f s\n", show_time);
        if (!(argc_info & ONLY_GREEDYNN)) {
            printf("2-Opt time: %f s\n", opt2_time);
        }
    }
    /*check_tour_validity(rota, d.numnodes);*/

    free(rota);
exit_kdtree:
    kdtree_free(&tree);
exit:
    free(d.x);
    free(d.y);

    return result;
}
