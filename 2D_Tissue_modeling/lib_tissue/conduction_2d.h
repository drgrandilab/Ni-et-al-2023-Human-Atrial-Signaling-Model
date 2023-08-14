#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <iostream>

// #include "conduction.h"
#include "util.h"
#include <assert.h>
#ifndef CONDUCTION_2D_H
#define CONDUCTION_2D_H


double **generate_laplacian_heterogeneity_2d ( unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double d2);

double **generate_laplacian_heterogeneity_2d_test(unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double* d2);

double **generate_laplacian_with_diffusion_heterogeneity_2d(unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy,
        double **dd, double **d2);



void calculate_fibre_unit_vectors_2d(unsigned char ** phi,
                                     int ny, int nx, float **y, float **x);



void write_fibre_vtk_2d(float ** y, float ** x, int ny, int nx, double dx = 0.3);




/* compute the diffusion tensor */
void calculate_diffusion_tensor_components_with_diffusion_heterogeneity_2d(unsigned char ** atrium,
        int ny, int nx,
        float **y, float **x,
        double dx, double dy,
        double **dd, double **d2,
        double  * **dc, double  * **df,
        double  * **d_dd, double  * **d_d2);

double * * *generate_2d_laplacian_with_diffusion_heterogeneity(unsigned char ** atrium,
        int ny, int nx,
        const char *phi_name, 
        double dx, double dy,
        double **dd, double **d2);


int ** generate_numeric_map_2d(unsigned char ** atrium, int ny, int nx, int * count) ;


int ** allocate_numeric_map_2d(int count);


int ** generate_neighbours_from_numeric_map_2d(int ** map, int ny, int nx, int count);


int ** generate_neighbours_map_2d(unsigned char ** atrium, int ny, int nx, int * count) ;



void calculate_laplacian_components_with_diffusion_heterogeneity_2d ( unsigned char ** atrium,
        int ny, int nx,
        float **y, float **x,
        double dx, double dy,
        double **dd, double **d2,
        double  ***dc, double ***df,
        double  ***d_dd, double ***d_d2,
        double ***lap);
unsigned char * generate_tissue_map_2d(unsigned char ** atrium,
                                 int ny, int nx, int count);
unsigned char * create_stimulation_map_vec_2D_geo(unsigned char ** atrium, const char *filename, int ny, int nx, int count);

void update_2D_laplacian(double * lap, int num_neighboor, unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double* d2);
#endif