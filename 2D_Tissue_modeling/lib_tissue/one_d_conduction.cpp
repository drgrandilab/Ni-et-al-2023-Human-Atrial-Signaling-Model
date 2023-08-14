/*
 *  one_d_conduction.c
 *
 *  generate lap and nbh information for purely one dimensional
 *  simulations.
 *  assume that FE were to be used.
 *
 *   Haibo Ni <qiangzini@gmail.com>
 *   Dec 09 2014
 *   Last Update :Fri 30 Jan 2015 15:38:55 GMT
 *   To include
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "util.h"

/* cell map info */
/*
    0 - self
    1 - x + 1
    2 - x - 1

    just to make it consistent with the 3D version code conduction.c

*/

int **generate_one_D_neighboors(int cell_num) {
    int **nbd;
    int i, j;

    // nbd = (int **) malloc(cell_num * sizeof(int *));
    nbd = new int* [cell_num];
    if (!nbd) {
        std::cerr << "Laplacian array malloc \n";
        std::exit(0);
    }
    // allocate the neighbours.
    for (i = 0; i < cell_num; ++i) {
        // nbd[i] = (int *) malloc(3 * sizeof(int));
        nbd[i] = new int [3];
        if (!nbd[i]) {
            std::cerr << "Laplacian array malloc \n";
            std::exit(0);
        }
        nbd[i][0] = i;
        nbd[i][1] = i + 1;
        nbd[i][2] = i - 1;
    }
    nbd[cell_num - 1][1] = cell_num - 1;
    nbd[0][2] = 0;
    return nbd;
}

double *get_ddx(double *diff, int num, double dx) {
    double *d_diff;
    int i, j;

    // allocate a long enough spine.
    // d_diff = (double *) malloc(num * sizeof(double));
    d_diff = new double [num];
    if (!d_diff) {
        std::cerr << "d_diff array malloc wrong!\n";
        std::exit(0);
    }

    for (i = 1; i < num - 1; i++) {
        d_diff[i] = (diff[i + 1] - diff[i - 1]) / (2 * dx);
    }

    d_diff[0] = (diff[1] - diff[0]) / (2 * dx);
    d_diff[num - 1] = (diff[num - 1] - diff[num - 2]) / (2 * dx);
    return d_diff;
}


double **generate_one_D_laplacian(int cell_num, double *diff, double *d_diff, double dx) {
    // allocates a 2d laplacian array.

    double **lap;
    int i, j;

    // allocate a long enough spine.
    // lap = (double **) malloc(cell_num * sizeof(double *));
    lap = new double* [cell_num];
    if (!lap) {
        std::cerr << "Laplacian array malloc\n";
        std::exit(0);
    }

    // allocate the neighbours.
    for (i = 0; i < cell_num; ++i) {
        // lap[i] = (double *) malloc(3 * sizeof(double));
        lap[i] = new double [3];
        if (!lap[i]) {
            std::cerr << "Laplacian array malloc\n";
            std::exit(0);
        }
        lap[i][0] = 0.0;
        lap[i][1] = 0.0;
        lap[i][2] = 0.0;

    }
    double factor;
    for (i = 0; i < cell_num; ++i) {
        factor = diff[i] / (dx * dx);
        lap[i][0] += -2.0 * factor;
        lap[i][1] += 1.0 * factor;
        lap[i][2] += 1.0 * factor;
    }

    for (i = 0; i < cell_num; ++i) {
        factor = d_diff[i] / (2 * dx);

        // upwind scheme
        if (factor > 0 ) {
            lap[i][0] += -2.0 * factor;
            lap[i][1] += 2.0 * factor;
            lap[i][2] += 0.0 * factor;
        } else {
            lap[i][0] += 2.0 * factor;
            lap[i][1] += 0.0 * factor;
            lap[i][2] += -2.0 * factor;
        }

        // centered difference
        /*
            lap[i][0] += 0.0 * factor;
            lap[i][1] += 1.0 * factor;
            lap[i][2] += -1.0 * factor;
        */
    }

    return lap;
}


double ** generate_one_D_laplacian_test(int cell_num, unsigned char *tissue, double diff_coef, double dx) {
    double *diff = new double [cell_num];
    if (!diff) {
        std::cerr << "diff array malloc failed!\n";
        std::exit(0);
    }
    int i;
    for (i = 0; i < cell_num; i++) {
        diff[i] = diff_coef;
        if (fabs(i - 60) <= 2 ) {
            // diff[i] = diff_coef / 5.0;
            // std::cout << "yese\n";
        }

    }

    double *d_diff = get_ddx(diff, cell_num, dx);
    double ** lap = generate_one_D_laplacian(cell_num, diff, d_diff, dx);
    delete [] diff;
    delete [] d_diff;
    return lap;
}


double ** generate_one_D_laplacian_test(int cell_num, unsigned char *tissue, double * diff, double dx) {
    if (!diff) {
        std::cerr << "diff array malloc failed!\n";
        std::exit(0);
    }
    double *d_diff = get_ddx(diff, cell_num, dx);
    double ** lap = generate_one_D_laplacian(cell_num, diff, d_diff, dx);
    delete [] d_diff;
    return lap;
}


void update_one_D_laplacian(double * lap, int num_neighboor, int cell_num, unsigned char *tissue, double * diff, double dx) {

    double** laplacian = generate_one_D_laplacian_test(cell_num, tissue, diff, dx);
    for (int n = 0; n < cell_num; n++) {
        for (int i = 0; i < num_neighboor; i++) {
            lap[(n * num_neighboor) + i] = laplacian[n][i];
        }
    }
    deallocate_and_zero_2d_matrix(laplacian, cell_num, num_neighboor);
}