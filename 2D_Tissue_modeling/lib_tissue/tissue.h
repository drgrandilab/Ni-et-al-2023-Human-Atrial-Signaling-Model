/*
 * tissue.h
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 *
 * Functions dealing with mapping or interacting with the tissue.
 */

#ifndef TISSUE_H
#define TISSUE_H
#include <string>

int ** generate_neighbours_map(unsigned char *** atrium, int nz, int ny, int nx, int * count);
// generates a numeric map of the tissue for all nodes which have an active
// cell.  The result is an ncell*nneighbours (19), sized array, where all cells
// have the indicies of all their neighbours.

unsigned char * generate_tissue_map(unsigned char *** atrium, int nz, int ny, int nx, int count);
// makes a linear map of the tissue.

int *** generate_numeric_map(unsigned char *** atrium, int nz, int ny, int nx, int * count);
// generates a numeric map of the tissue, with all tissue nodes assigned a
// number, and all non-tissue nodes assigned -1.  Also returns a count of the
// cells.

int ** generate_neighbours_from_numeric_map(int *** map, int nz, int ny, int nx, int count);
// Generates a neighbours map from the given numeric map.

int ** allocate_numeric_map(int count);
// allocates a numeric map, assigning indexes as -1.

// create 2d patch of cells of repeated pattern
int * generate_patch_pattern_2d(unsigned char ** atrium, int ny, int nx, int count, int block_size);


// stimulation flag
unsigned char * create_stimulation_map_vec_band(unsigned char ** atrium, int ny, int nx, int count, std::string mode, int band_width);

unsigned char *  create_stimulation_map_vec_circle_2D(unsigned char ** atrium, int ny, int nx, int count, int ori_y, int ori_x, int radius);

#endif /* TISSUE_H */

