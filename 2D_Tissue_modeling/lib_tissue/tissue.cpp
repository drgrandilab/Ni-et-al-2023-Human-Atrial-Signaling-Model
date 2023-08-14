/*
 * tissue.c
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 *


 *  Update by Haibo Ni
 *  qiangzi.ni@gmail.com
 *  Dec 07 2014
 */
#include "util.h"
#include "tissue.h"
#include <stdlib.h>
#include <stdio.h>
#include <string>

int ** generate_neighbours_map(unsigned char *** atrium, int nz, int ny, int nx, int * count) {
// generates a numeric map of the tissue for all nodes which have an active
// cell.  The result is an ncell*nneighbours (19), sized array, where all cells
// have the indicies of all their neighbours.

    int ***map;
    int ** nbd;

    // generate the numeric map we'll need.
    map = generate_numeric_map(atrium, nz, ny, nx, count);

    // generate the neighbours from it.
    nbd = generate_neighbours_from_numeric_map(map, nz, ny, nx, *count);

    // free the map.
    deallocate_and_zero_3d_int(map, nz + 2, ny + 2, nx + 2);

    // return the neighbours
    return nbd;
}

unsigned char * generate_tissue_map(unsigned char *** atrium,
                                    int nz, int ny, int nx, int count) {
// makes a linear map of the tissue.
    unsigned char * tissue;
    int i, j, k, c;

    // allocate space for the tissue.
    // tissue = (unsigned char*) malloc(count*sizeof(unsigned char));
    tissue = new unsigned char [count];
    if (!tissue) {
        perror("malloc: tissue map arrray");
        exit(EXIT_FAILURE);
    }

    c = 0;
    for (k = 1; k < nz + 1; k++)
        for (j = 1; j < ny + 1; j++)
            for (i = 1; i < nx + 1; i++)
                if (atrium[k][j][i] > 0) {
                    tissue[c] = atrium[k][j][i];
                    c++;
                }
    if (c != count) {
        perror("in generating tissue map, the total num of cells does not equal to the counter of input");
        exit(EXIT_FAILURE);
    }
    return tissue;
}

int ** generate_neighbours_from_numeric_map(int *** map, int nz, int ny, int nx, int count) {
    // Generates a neighbours map from the given numeric map.
    int ** nbd;
    int n, i, j, k;

    // allocate the map.
    nbd = allocate_numeric_map(count);

    // iterate over the map, assigning neighbours.
    // self         = 0
    // xplus        = 1
    // xminus       = 2
    // yplus        = 3
    // yminus       = 4
    // zplus        = 5
    // zminus       = 6
    // xplusyplus   = 7
    // xminusyminus = 8
    // xplusyminus  = 9
    // xminusyplus  = 10
    // xpluszplus   = 11
    // xminuszminus = 12
    // xpluszminus  = 13
    // xminuszplus  = 14
    // ypluszplus   = 15
    // yminuszminus = 16
    // ypluszminus  = 17
    // yminuszplus  = 18
    n = 0;
    for (k = 1; k < nz + 1; k++)
        for (j = 1; j < ny + 1; j++)
            for (i = 1; i < nx + 1; i++) {
                if (map[k][j][i] > -1) {
                    nbd[n][ 0] = map[k][j][i];

                    // simple directions.
                    // all get assigned the real answer, or self.
                    nbd[n][ 1] = (map[k][j][i + 1] > -1) ? map[k][j][i + 1] : map[k][j][i];
                    nbd[n][ 2] = (map[k][j][i - 1] > -1) ? map[k][j][i - 1] : map[k][j][i];
                    nbd[n][ 3] = (map[k][j + 1][i] > -1) ? map[k][j + 1][i] : map[k][j][i];
                    nbd[n][ 4] = (map[k][j - 1][i] > -1) ? map[k][j - 1][i] : map[k][j][i];
                    nbd[n][ 5] = (map[k + 1][j][i] > -1) ? map[k + 1][j][i] : map[k][j][i];
                    nbd[n][ 6] = (map[k - 1][j][i] > -1) ? map[k - 1][j][i] : map[k][j][i];

                    // now for more complex directions.
                    nbd[n][ 7] = (map[k][j + 1][i + 1] > -1) ? map[k][j + 1][i + 1] : map[k][j][i];
                    nbd[n][ 8] = (map[k][j - 1][i - 1] > -1) ? map[k][j - 1][i - 1] : map[k][j][i];
                    nbd[n][ 9] = (map[k][j - 1][i + 1] > -1) ? map[k][j - 1][i + 1] : map[k][j][i];
                    nbd[n][10] = (map[k][j + 1][i - 1] > -1) ? map[k][j + 1][i - 1] : map[k][j][i];
                    nbd[n][11] = (map[k + 1][j][i + 1] > -1) ? map[k + 1][j][i + 1] : map[k][j][i];
                    nbd[n][12] = (map[k - 1][j][i - 1] > -1) ? map[k - 1][j][i - 1] : map[k][j][i];
                    nbd[n][13] = (map[k - 1][j][i + 1] > -1) ? map[k - 1][j][i + 1] : map[k][j][i];
                    nbd[n][14] = (map[k + 1][j][i - 1] > -1) ? map[k + 1][j][i - 1] : map[k][j][i];
                    nbd[n][15] = (map[k + 1][j + 1][i] > -1) ? map[k + 1][j + 1][i] : map[k][j][i];
                    nbd[n][16] = (map[k - 1][j - 1][i] > -1) ? map[k - 1][j - 1][i] : map[k][j][i];
                    nbd[n][17] = (map[k - 1][j + 1][i] > -1) ? map[k - 1][j + 1][i] : map[k][j][i];
                    nbd[n][18] = (map[k + 1][j - 1][i] > -1) ? map[k + 1][j - 1][i] : map[k][j][i];

                    n++;
                }
            }

    return nbd;
}

int ** allocate_numeric_map(int count) {
// allocates a numeric map, assigning indexes as -1.

    int ** nbd;
    int i, j;

    // allocate a long enough spine.
    // nbd = (int**) malloc(count*sizeof(int*));
    nbd = new int * [count];
    if (!nbd) {
        perror("Neighbours array malloc");
        exit(EXIT_FAILURE);
    }

    // allocate the neighbours.
    for (i = 0; i < count; ++i) {
        // nbd[i] = (int*) malloc(19*sizeof(int));
        nbd[i] = new int [19];
        if (!nbd[i]) {
            perror("Neighbours array malloc");
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < 19; ++j) {
            nbd[i][j] = -1;
        }
    }

    return nbd;
}



int *** generate_numeric_map(unsigned char *** atrium, int nz, int ny, int nx, int * count) {
// generates a numeric map of the tissue, with all tissue nodes assigned a
// number, and all non-tissue nodes assigned -1.  Also returns a count of the
// cells.
    int *** map = allocate_and_set_3d_int(nz + 2, ny + 2, nx + 2, -1);
    int i, j, k;
    int c = 0;

    for (k = 1; k < nz + 1; k++)
        for (j = 1; j < ny + 1; j++)
            for (i = 1; i < nx + 1; i++) {
                if (atrium[k][j][i] > 0) {
                    map[k][j][i] = c;
                    c++;
                }
            }

    *count = c;
    return map;
}



int * generate_patch_pattern_2d(unsigned char ** atrium, int ny, int nx, int count, int block_size) {

    if (block_size <= 0) {

        std::cerr << "block_size cannot be zero or negative!!!! generate_patch_pattern_2d, tissue.h" << std::endl;
        std::exit(0);

    }
    int * block_id = new int[count];
    int i, j, k;
    int c = 0;

    for (j = 1; j < ny + 1; j++)
        for (i = 1; i < nx + 1; i++) {
            if (atrium[j][i] > 0) {
                // map[j][i] = c;
                c++;
            }
        }

    if (c != count) {
        std::cerr << "Error, c != count, generate_patch_pattern_2d, tissue.h" << std::endl;
        std::exit(0);
    }
    c = 0;
    for (j = 1; j < ny + 1; j++)
        for (i = 1; i < nx + 1; i++) {
            if (atrium[j][i] > 0) {
                // map[j][i] = c;

                block_id[c] = ((i - 1) / block_size) + (nx / block_size) * ((j - 1) / block_size);
                c++;
            }
        }

    return block_id;

}




unsigned char * create_stimulation_map_vec_band(unsigned char ** atrium, int ny, int nx, int count, std::string mode, int band_width ) {

    // unsigned char ** stim = read_and_embed_geometry_2d(filename, ny, nx);

    // makes a linear map of the tissue.

    int i, j, c;

    // allocate space for the stim_vec.
    // stim_vec = (unsigned char*) malloc(count*sizeof(unsigned char));
    unsigned char * stim_vec = new unsigned char [count];
    if (!stim_vec) {
        perror("new: stim_vec map arrray failed!\n");
        exit(EXIT_FAILURE);
    }

    c = 0;
    for (j = 1; j < ny + 1; j++)
        for (i = 1; i < nx + 1; i++)
            if (atrium[j][i] > 0) {

                stim_vec[c] = 0;
                // stim_vec[c] = stim[j][i];
                if (mode == "left" or mode == "LEFT") {
                    if (i <= band_width)
                        stim_vec[c] = 1;
                } else if (mode == "right" or mode == "RIGHT") {
                    if (i >= nx + 1 - band_width)
                        stim_vec[c] = 1;
                } else if (mode == "bottom" or mode == "BOTTOM") {
                    if (j <= band_width)
                        stim_vec[c] = 1;
                } else if (mode == "top" or mode == "TOP") {
                    if (j >= ny + 1 - band_width)
                        stim_vec[c] = 1;
                } else {
                    std::cerr << " wrong mode detected !!! create_stimulation_map_vec_band" << std::endl;
                    // delete [] stim_vec;
                    std::exit(0);
                }

                c++;
            }
    if (c != count) {
        perror("in generating stim_vec map, the total num of cells does not equal to the counter of input! \n");
        exit(EXIT_FAILURE);
    }


    // deallocate_and_zero_2d_matrix(stim, ny + 2, nx + 2);
    return stim_vec;

}



unsigned char * create_stimulation_map_vec_circle_2D(unsigned char ** atrium, int ny, int nx, int count, int ori_y, int ori_x, int radius) {

    // unsigned char ** stim = read_and_embed_geometry_2d(filename, ny, nx);

    // makes a linear map of the tissue.

    int i, j, c;

    // allocate space for the stim_vec.
    // stim_vec = (unsigned char*) malloc(count*sizeof(unsigned char));
    unsigned char * stim_vec = new unsigned char [count];
    if (!stim_vec) {
        perror("new: stim_vec map arrray failed!\n");
        exit(EXIT_FAILURE);
    }

    c = 0;
    for (j = 1; j < ny + 1; j++)
        for (i = 1; i < nx + 1; i++)
            if (atrium[j][i] > 0) {

                stim_vec[c] = 0;


                if ((j - ori_y) * (j - ori_y) + (i - ori_x) * (i - ori_x) <= radius * radius)
                    stim_vec[c] = 1;


                // // stim_vec[c] = stim[j][i];
                // if (mode == "left" or mode == "LEFT") {
                //     if (i <= band_width)
                //         stim_vec[c] = 1;
                // } else if (mode == "right" or mode == "RIGHT") {
                //     if (i >= nx + 1 - band_width)
                //         stim_vec[c] = 1;
                // } else if (mode == "bottom" or mode == "BOTTOM") {
                //     if (j <= band_width)
                //         stim_vec[c] = 1;
                // } else if (mode == "top" or mode == "TOP") {
                //     if (j >= ny + 1 - band_width)
                //         stim_vec[c] = 1;
                // } else {
                //     std::cerr << " wrong mode detected !!! create_stimulation_map_vec_band" << std::endl;
                //     // delete [] stim_vec;
                //     std::exit(0);
                // }

                c++;
            }
    if (c != count) {
        perror("in generating stim_vec map, the total num of cells does not equal to the counter of input! \n");
        exit(EXIT_FAILURE);
    }


    // deallocate_and_zero_2d_matrix(stim, ny + 2, nx + 2);
    return stim_vec;

}