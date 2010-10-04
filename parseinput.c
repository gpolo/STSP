/*
 * Created:  19/08/2010 16:10:17
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parseinput.h"
#include "dist.h"
#include "util.h"

int
get_info(struct Info *si)
{
    size_t i;
    char buffer[150];

    si->edge_weight_format = NONE;
    si->problem_type = -1;
    si->name = NULL;

    while (1) {
        if (fgets(buffer, sizeof(buffer), stdin) == NULL) {
            perror("fgets");
            goto err;
        }
        /* Remove '\n'. */
        buffer[strlen(buffer) - 1] = '\0';

        if (strstr(buffer, "_SECTION")) {
            /* Dados seguem. */
            break;
        }

        if (strlen(buffer) == 0) {
            fprintf(stderr, "ERROR: Input file incorrectly formatted.\n");
            goto err;
        }

        for (i = 0; i < strlen(buffer); i++) {
            if (*buffer == 'N' && *(buffer + 1) == 'A') {
                /* NAME: */
                for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                i += 2;
                si->name = MALLOC(strlen(buffer + i) + 1, char);
                strcpy(si->name, buffer + i);
                break;
            } else if (*buffer == 'D') {
                if (*(buffer + 2) == 'M') {
                    /* DIMENSION: */
                    for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                    i += 2;
                    si->dimension = atoi(buffer + i);
                    break;
                } else if (*(buffer + 2) == 'S') {
                    /* DISPLAY_DATA_TYPE: */
                    for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                    i += 2;
                    if (*(buffer + i) == 'C') {
                        si->display_data_type = COORD_DISPLAY;
                        break;
                    } else if (*(buffer + i) == 'T') {
                        si->display_data_type = TWOD_DISPLAY;
                        break;
                    } else {
                        fprintf(stderr, "ERROR: Display data type '%s' is not "
                                        "supported.\n", (buffer + i));
                        goto err;
                    }
                    break;
                }
            } else if (*buffer == 'E' && *(buffer + 12) == 'T') {
                /* EDGE_WEIGHT_TYPE: */
                for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                i += 2;
                if (*(buffer + i) == 'G') {
                    si->edge_weight_type = GEO;
                } else if (*(buffer + i) == 'C') {
                    si->edge_weight_type = CEIL_2D;
                } else if (*(buffer + i) == 'E' && *(buffer + i + 4) == '2') {
                    si->edge_weight_type = EUC_2D;
                } else if (*(buffer + i) == 'E' && *(buffer + i + 1) == 'X') {
                    si->edge_weight_type = EXPLICIT;
                } else if (*(buffer + i) == 'A') {
                    si->edge_weight_type = ATT;
                } else {
                    fprintf(stderr, "ERROR: Edge weight type '%s' is not "
                                    "supported.\n", buffer + i);
                    goto err;
                }
                break;
            } else if (*buffer == 'E' && *(buffer + 12) == 'F') {
                /* EDGE_WEIGHT_FORMAT: */
                for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                i += 2;
                if (*(buffer + i) == 'F' && *(buffer + i + 1) == 'U') {
                    si->edge_weight_format = FUNCTION;
                } else if (*(buffer + i) == 'U' && *(buffer + i + 6) == 'R') {
                    si->edge_weight_format = UPPER_ROW;
                } else if (*(buffer + i) == 'L' && strlen(buffer + i) > 9 &&
                        *(buffer + i + 11) == 'R') {
                    si->edge_weight_format = LOWER_DIAG_ROW;
                } else {
                    fprintf(stderr, "ERROR: Edge weight format '%s' is not "
                                    "supported.\n", buffer + i);
                    goto err;
                }
                break;
            } else if (*buffer == 'T') {
                for (; i < strlen(buffer) && buffer[i] != ':'; i++) ;
                i += 2;
                if (strcmp(buffer + i, "TSP")) {
                    fprintf(stderr, "ERROR: Expected a TSP problem, "
                                    "got '%s'.\n", buffer + i);
                    goto err;
                } else {
                    si->problem_type = TSP;
                    break;
                }
            } else {
                /* Descarta campo. */
                break;
            }
        }
    }

    goto success;

err:
    if (si->name) {
        free(si->name);
    }
    return -1;

success:
    return 0;
}

int
setdist_readdata(struct Info *info, data *d)
{
    int i;
    int status = 0;

    switch (info->edge_weight_type) {
        case CEIL_2D:
            d->dist = dist_ceil2d;
            goto read;
        case EUC_2D:
            d->dist = dist_euc2d;
            goto read;
        case ATT:
            d->dist = dist_att;
read:
            for (i = 0; i < info->dimension; i++) {
                if (scanf("%*d %lf %lf", &d->x[i], &d->y[i]) != 2) {
                    fprintf(stderr, "Incorrect data format.\n");
                    abort();
                }
            }
            break;

        case GEO:
            d->dist = dist_geo;
            for (i = 0; i < info->dimension; i++) {
                if (scanf("%*d %lf %lf", &d->x[i], &d->y[i]) != 2) {
                    fprintf(stderr, "Incorrect data format.\n");
                    abort();
                }

                geo_latlon(&d->x[i], &d->y[i]);
            }
            break;

        case EXPLICIT:
            status = read_explicit(info, d);
            break;
    }

    return status;
}

int
read_explicit(struct Info *info, data *d)
{
    int i, j;
    int status = 0;

    if (info->edge_weight_type != EXPLICIT) {
        /* Assumindo que realmente quer tratar o problema sem
         * utilizar kdtree. */
        return setdist_readdata(info, d);
    }

    if (info->dimension > MAX_SMALLDIM) {
        fprintf(stderr, "ERROR: Dimension > %d is not supported for EXPLICIT "
                        "edge weight.\n", MAX_SMALLDIM);
        status = -1;
    }


    d->dist = matrix_at;

    switch (info->edge_weight_format) {
    case UPPER_ROW:
        for (i = 0; i < info->dimension; i++) {
            for (j = i + 1; j < info->dimension; j++) {
                if (scanf("%d", &small_matrix[i][j]) != 1)
                    abort();
                small_matrix[j][i] = small_matrix[i][j];
            }
        }
        break;

    case LOWER_DIAG_ROW:
        for (i = 0; i < info->dimension; i++) {
            for (j = 0; j <= i; j++) {
                if (scanf("%d", &small_matrix[i][j]) != 1)
                    abort();
                small_matrix[j][i] = small_matrix[i][j];
            }
        }
        break;

    default:
        fprintf(stderr, "ERROR: The declared edge weight format is not "
                        "supported.\n");
        status = -1;
        break;
    }

    return status;
}
