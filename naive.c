/*
 * Created:  24/08/2010 16:53:13
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "parseinput.h"
#include "util.h"
#include "naive.h"
#include "tsp.h"

#include <math.h>

/*#define C(t, d, a, b) (cached_dist(d, t[a], t[b]))*/
#define C(t, d, a, b) (d->dist(d, t[a], t[b]))

#define C3(t, dat, a, b, c, d, e, f) (C(t, dat, a, b) + C(t, dat, c, d) + C(t, dat, e, f))
#define C2(t, dat, a, b, c, d) (C(t, dat, a, b) + C(t, dat, c, d))


#define DO_2OPT(tour, a, b) do { \
    invert(tour, a, b); \
} while(0)


static void
invert(int *array, int start, int end)
{
    int temp;

    for (; start <= end; start++) {
        temp = array[start];
        array[start] = array[end];
        array[end] = temp;
        end--;
    }
}


void
local_2opt(data *dat, int *tour, int *cost, int *swaps)
{
    int dim = dat->numnodes;
    int i, j;
    int a, b, c, d;
    int ganho, melhor_ganho, melhor_b = -1, melhor_c = -1;
    int trocas = 0, ganho_total = 0;

    int consid = 0;

    while (1) {
        melhor_ganho = 0;

        for (i = 0; i < dim; i++) {
            a = i;
            b = (i + 1) % dim;

            for (j = i + 2; j < dim; j++) {
                c = j;
                d = (j + 1) % dim;
                if (d == a) {
                    continue;
                }
                consid++;

                /* 2-opt */
                ganho = C(tour, dat, a, c) + C(tour, dat, b, d);
                ganho -= C(tour, dat, a, b) + C(tour, dat, c, d);

                if (ganho < melhor_ganho) {
                    melhor_b = b;
                    melhor_c = c;
                    melhor_ganho = ganho;
                }
            }
        }

        if (melhor_ganho < 0) {
            assert(melhor_b >= 0 && melhor_c >= 0);
            DO_2OPT(tour, melhor_b, melhor_c);
            ganho_total += melhor_ganho;
            trocas++;
        } else {
            break;
        }
    }

    *cost += ganho_total;
    *swaps = trocas;
    DEBUG("Consideraçoes = %d\n", consid);
}

/* XXX Works but it is not used. */
#if 0
#define SAT_SUBPATH_REVERSE(sat, a, b, c, d) do { \
    sat[a] = c ^ 1; sat[c] = a ^ 1; \
    sat[d ^ 1] = b; sat[b ^ 1] = d; \
} while(0)

static void
local_2opt_satellite(data *dat, int *tour, int *outt, int *cost, int *swaps)
{
    int dim = dat->numnodes;
    int i, j;
    int a, b, c, d;
    int ganho, melhor_ganho;
    int trocas = 0, ganho_total = 0;
    int troca[4];

    int consid = 0;

    /* Satellite list init. */
    int sd, val; /* XXX Used below. */
    /*int sat[dim * 2];*/
    int *sat = MALLOC(dim * 2, int);
    sat[0] = 2;
    sat[1] = dim * 2 - 1;
    sat[dim * 2 - 2] = 0;
    sat[dim * 2 - 1] = dim * 2 - 3;
    for (i = 2; i < dim * 2 - 2; i += 2) {
        sat[i] = i + 2;
        sat[i + 1] = i - 1;
    }
    /* End. */

    while (1) {
        melhor_ganho = 0;

        a = 0; /* sat index */
        for (i = 0; i < dim - 2; i++) {
            b = sat[a]; /* b é vizinho de a. */

            c = sat[b]; /* c é vizinho de b. */
            for (j = i + 2; j < dim; j++) {
                d = sat[c]; /* d é vizinho de c. */

                /* XXX Falta alguma verificação para Satellite ? */
                if (d == a)
                    break;

                consid++;

                /* 2-opt */
                ganho = C(tour, dat, a/2, c/2) + C(tour, dat, b/2, d/2);
                ganho -= C(tour, dat, a/2, b/2) + C(tour, dat, c/2, d/2);

                if (ganho < melhor_ganho) {
                    troca[0] = a; troca[1] = b;
                    troca[2] = c; troca[3] = d;
                    melhor_ganho = ganho;
                }

                c = d;
            }
            a = b;
        }

        if (melhor_ganho < 0) {
            SAT_SUBPATH_REVERSE(sat,
                    troca[0], troca[1], troca[2], troca[3]);
            ganho_total += melhor_ganho;
            trocas++;
        } else {
            break;
        }
    }

    /* Atualizar tour agora. */
    val = 0;
    i = 0;
    do {
        sd = sat[val];
        outt[i] = tour[val / 2];
        val = sd;
    } while (++i < dim);

    *cost += ganho_total;
    *swaps = trocas;
    DEBUG("Considerações = %d\n", consid);

    free(sat);
}
#endif

#define DEFINE_TROCA_CIDADES(cidades, a, b, c, d, e, f) do { \
    cidades[0] = a; cidades[1] = b; cidades[2] = c; \
    cidades[3] = d; cidades[4] = e; cidades[5] = f; \
} while (0)

#define TROCA_BC_DE     0 /* 2 2-opt */
#define TROCA_BC        2 /* 1 2-opt */
#define TROCA_DE        3 /* 1 2-opt */
#define TROCA_BE        4 /* 1 2-opt */
#define TROCA_2_2OPT_BE_ED  8
#define TROCA_2_2OPT_BC_CE  9
#define TROCA_3_2OPT_BC_CE_ED 16

void
local_3opt(data *dat, int *tour, int *cost, int *swaps)
{
    int dim = dat->numnodes;
    int i, j, k, temp;
    int a, b, c, d, e, f;
    int ganho, ganho2, melhor_ganho;
    int custo3, cidades[6];
    int tipo_troca, tipo_troca_final; /* XXX Ver flags acima. */
    int trocas = 0, ganho_total = 0;

    while (1) {
        melhor_ganho = 0;
        tipo_troca_final = -1;
        for (i = 0; i < dim; i++) {
            a = i;
            b = (i + 1) % dim;

            for (j = i + 2; j < dim; j++) {
                c = j;
                d = (j + 1) % dim;
                if (c == a || c == b || d == a || d == c)
                    continue;

                for (k = j + 2; k < dim; k++) {
                    e = k;
                    f = (k + 1) % dim;
                    if (e == a || e == b || e == c || e == d ||
                            f == a || f == b || f == c || f == d)
                        continue;

                    /* XXX Papel. */

                    custo3 = C3(tour, dat, a, b, c, d, e, f);

                    ganho = 0;

                    /* 0 */
                    ganho2 = C3(tour, dat, a, c, b, e, d, f) - custo3;
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_BC_DE;
                        ganho = ganho2;
                    }

                    /* 16 */
                    ganho2 = C3(tour, dat, a, d, e, b, c, f) - custo3;
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_3_2OPT_BC_CE_ED;
                        ganho = ganho2;
                    }

                    /* 8 */
                    ganho2 = C3(tour, dat, a, d, e, c, b, f) - custo3;
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_2_2OPT_BE_ED;
                        ganho = ganho2;
                    }

                    /* 9 */
                    ganho2 = C3(tour, dat, a, e, d, b, c, f) - custo3;
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_2_2OPT_BC_CE;
                        ganho = ganho2;
                    }

                    /* 2 */
                    ganho2 = C2(tour, dat, a,c,b,d) - C2(tour, dat, a,b,c,d);
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_BC;
                        ganho = ganho2;
                    }

                    /* 3 */
                    ganho2 = C2(tour, dat, c,e,d,f) - C2(tour, dat, c,d,e,f);
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_DE;
                        ganho = ganho2;
                    }

                    /* 4 */
                    ganho2 = C2(tour, dat, a,e,b,f) - C2(tour, dat, a,b,e,f);
                    if (ganho2 < ganho) {
                        tipo_troca = TROCA_BE;
                        ganho = ganho2;
                    }

                    if (ganho < melhor_ganho) {
                        DEFINE_TROCA_CIDADES(cidades, a, b, c, d, e, f);
                        melhor_ganho = ganho;
                        tipo_troca_final = tipo_troca;
                    }
                }
            }
        }

        if (melhor_ganho < 0) {
            /*printf("melhor_ganho = %d, tipo_troca = %d\n",
                    melhor_ganho, tipo_troca_final);*/
            assert(tipo_troca_final != -1);
            switch (tipo_troca_final) {
            case TROCA_BC: /* 2 */
                DO_2OPT(tour, cidades[1], cidades[2]);
                break;

            case TROCA_BE: /* 4 */
                DO_2OPT(tour, cidades[1], cidades[4]);
                break;

            case TROCA_DE: /* 3 */
                DO_2OPT(tour, cidades[3], cidades[4]);
                break;

            case TROCA_BC_DE: /* 0 */
                /* Equivale a dois 2-opt, B<->C e D<->E. */
                DO_2OPT(tour, cidades[1], cidades[2]);
                DO_2OPT(tour, cidades[3], cidades[4]);
                break;

            case TROCA_2_2OPT_BE_ED: /* 8 */
                /* Primeiro 2-opt, trocar B por E. */
                DO_2OPT(tour, cidades[1], cidades[4]);

                /* Segundo 2-opt, trocar novo E por novo D. */
                temp = cidades[4];
                /* Troquei B por E antes, logo novo E esta em B agora. */
                cidades[4] = cidades[1];
                /* Distância de E a D se mantém. */
                cidades[3] = cidades[4] + (temp - cidades[3]);

                DO_2OPT(tour, cidades[4], cidades[3]);
                break;

            case TROCA_2_2OPT_BC_CE:
                /* Primeiro 2-opt, trocar B por C. */
                DO_2OPT(tour, cidades[1], cidades[2]);

                /* Segundo 2-opt, trocar novo C por E (não se alterou). */
                /* Troquei B por C antes, logo novo C está em B agora. */
                cidades[2] = cidades[1];
                
                DO_2OPT(tour, cidades[2], cidades[4]);
                break;

            case TROCA_3_2OPT_BC_CE_ED:
                /* Primeiro 2-opt, trocar B por C. */
                DO_2OPT(tour, cidades[1], cidades[2]);

                /* Segundo 2-opt, trocar novo C por E (não se alterou). */
                /* Troquei B por C antes, logo novo C está em B agora. */
                cidades[2] = cidades[1];

                DO_2OPT(tour, cidades[2], cidades[4]);

                /* Terceiro 2-opt, trocar novo E por D. */
                temp = cidades[4];
                /* C estava no lugar de B antes da troca anterior,
                 * após trocar E com C, E passa estar no lugar de B. */
                cidades[4] = cidades[1];
                cidades[3] = cidades[4] + (temp - cidades[3]);

                DO_2OPT(tour, cidades[4], cidades[3]);
                break;
            } /* switch */

            ganho_total += melhor_ganho;
            trocas++;
        } else {
            break;
        }
    }

    *cost += ganho_total;
    *swaps = trocas;
}

static void
greedy(data *d, int start, int *tour, int *tour_length)
{
    int i, j;
    int dim = d->numnodes;
    int current, new_curr, min, len;
    int *visited = MALLOC(dim, int);
    int dist;

    memset(visited, 0, dim * sizeof(int));

    len = 0;
    current = start;
    for (i = 0; i < dim - 1; i++) {
        visited[current] = 1;
        tour[i] = current;

        min = INT_MAX;
        new_curr = -1;
        for (j = 0; j < dim; j++) {
            if (!visited[j]) {
                dist = d->dist(d, current, j);
                if (dist < min) {
                    min = dist;
                    new_curr = j;
                }
            }
        }
        assert(new_curr != -1);
        current = new_curr;
        len += min;
    }
    tour[i] = current;
    len += d->dist(d, current, start);
    *tour_length = len;

    free(visited);
}


int
naive(data *d, int flags)
{
    int dim = d->numnodes, start;
    int cost, trocas, twoopt_gain;
    int use_three_opt = flags & DO_3OPT;

    int *tour = MALLOC(dim, int);
    /*int i, outtour[dim];*/

    double start_time, greedy_time, show_time, opt2_time;

    start_time = current_usertime_secs();
    start = RANDOM(dim);
    greedy(d, start, tour, &cost);
    greedy_time = current_usertime_secs() - start_time;
    start_time = current_usertime_secs();
    if (flags & BE_VERBOSE) {
        show_sol(tour, dim, cost);
    }
    printf("\nStart point: %d\n", start);
    printf("Greedy-NN tour cost: %d\n", cost);
    show_time = current_usertime_secs() - start_time;

    if (!(flags & ONLY_GREEDYNN)) {
        twoopt_gain = cost;
        start_time = current_usertime_secs();
        if (use_three_opt) {
            local_3opt(d, tour, &cost, &trocas);
        } else {
            local_2opt(d, tour, &cost, &trocas);
        }
        /*show_sol(tour, dim, cost);*/
        /*local_2opt_satellite(d, tour, outtour, &cost, &trocas);
        show_sol(outtour, dim, cost);*/
        opt2_time = current_usertime_secs() - start_time;
        printf("Cost after running %c-opt: %d [%d swaps]\n",
                (use_three_opt) ? '3':'2', cost, trocas);
    }

    if (flags & BE_VERBOSE) {
        printf("\nGreedy tour time: %f s\n", greedy_time);
        printf("Display solution time: %f s\n", show_time);
        if (!(flags & ONLY_GREEDYNN)) {
            printf("%c-Opt time: %f s\n", (use_three_opt) ? '3':'2', opt2_time);
        }
    }

    free(tour);
    return 0;
}
