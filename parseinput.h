/*
 * Created:  23/08/2010 15:43:26
 *
 * Author:  Guilherme Polo, ggpolo@gmail.com
 *
 */
#ifndef PARSEINPUT_H
#define PARSEINPUT_H

typedef struct data data;

struct Info {
    char *name;
    int dimension;
    enum { TSP } problem_type;
    enum { CEIL_2D, EUC_2D, GEO, EXPLICIT, ATT } edge_weight_type;
    enum { NONE, FUNCTION, UPPER_ROW, LOWER_DIAG_ROW } edge_weight_format;
    enum { COORD_DISPLAY, TWOD_DISPLAY } display_data_type; /* XXX Não utilizado. */
};

struct data {
    double *x, *y;
    int numnodes;
    int (*dist)(data *, int, int);
};

/* Para problemas com edge_weight_type == EXPLICIT. */
#define MAX_SMALLDIM 128
extern int small_matrix[MAX_SMALLDIM][MAX_SMALLDIM];

int get_info(struct Info *);
int setdist_readdata(struct Info *, data *);
int read_explicit(struct Info *, data *);

#endif /* PARSEINPUT_H */
