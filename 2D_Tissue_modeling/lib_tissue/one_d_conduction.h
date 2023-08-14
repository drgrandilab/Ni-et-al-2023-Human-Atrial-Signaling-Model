/*
 *  one_d_conduction.c
 *
 *  generate lap and nbh information for purely one dimensional
 *  simulations.
 *  assume that FE were to be used.
 *
 *   Haibo Ni <qiangzini@gmail.com>
 *   Dec 09 2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>

 #ifndef ONE_D_CONDUCTOIN_H
 #define ONE_D_CONDUCTOIN_H


int **generate_one_D_neighboors(int cell_num);
double *get_ddx(double *diff, int num, double dx);
double **generate_one_D_laplacian(int cell_num, double *diff, double *d_diff, double dx);
double ** generate_one_D_laplacian_test(int cell_num, unsigned char *tissue, double diff_coef, double dx);
double ** generate_one_D_laplacian_test(int cell_num, unsigned char *tissue, double * diff_coef, double dx);
void update_one_D_laplacian(double * lap, int num_neighboor, int cell_num, unsigned char *tissue, double * diff, double dx);
#endif /*end of one_d_conduction.h*/