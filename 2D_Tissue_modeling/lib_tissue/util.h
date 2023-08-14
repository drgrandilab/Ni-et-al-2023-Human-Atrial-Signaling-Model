/*
 * util.h
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 *
 * Collection of 'utility' functions, i.e. those that are pretty generic.
 */

#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>

#include <iostream>
#include <cstdlib>
#include "matrix_mem.hpp"
#include <sstream>
#include <string.h>

   
template <class T>
int File_Not_Found(T name);

// template <typename T>
std::string get_str_from_float(double data);

float *** read_and_embed_float_data(const char * filename, int nz, int ny, int nx);

unsigned char *** read_and_embed_geometry(const char * filename, int nz, int ny, int nx);
unsigned char ** read_and_embed_geometry_2d(const char *filename,  int ny, int nx);
// reads an unsigned character geometry from a (possibly gzip'd) file called
// filename and embeds it in a geometry that is one larger in all directions.

unsigned char *** allocate_and_zero_3d_unsigned_char(int nz, int ny, int nx);
// allocates and zeros a 3D array of unsigned chars of the given dimensions.
// nz is the slowest dimension, nx the fastest.

void deallocate_and_zero_3d_unsigned_char(unsigned char *** array, int nz, int ny, int nx);
// deallocates and zeros a 3D array of unsigned chars of the given dimensions.
// nz is the slowest dimension, nx the fastest.
double *** allocate_and_zero_3d_double(int nz, int ny, int nx);
void deallocate_and_zero_3d_double(double *** array, int nz, int ny, int nx);

float *** allocate_and_zero_3d_float(int nz, int ny, int nx);
// allocates and zeros a 3D array of floats of the given dimensions.
// nz is the slowest dimension, nx the fastest.

void deallocate_and_zero_3d_float(float *** array, int nz, int ny, int nx);
// deallocates and zeros a 3D array of floats of the given dimensions.
// nz is the slowest dimension, nx the fastest.

double **** allocate_and_zero_3d_double_pointer(int nz, int ny, int nx);
// allocates and zeros a 3d array of pointers to doubles of the given dimensions.
// nz is the slowest dimension, nx the fastest.

void deallocate_and_zero_3d_double_pointer(double **** array, int nz, int ny, int nx);
// deallocates and zeros a 3d array of floats of the given dimensions.
// nz is the slowest dimension, nx the fastest.

int *** allocate_and_zero_3d_int(int nz, int ny, int nx);
// allocates and zeros a 3d array of ints of the given dimensions.
// nz is the slowest dimension, nx the fastest.

int *** allocate_and_set_3d_int(int nz, int ny, int nx, int set);
// allocates and sets a 3d array of ints of the given dimensions.
// nz is the slowest dimension, nx the fastest.

void deallocate_and_zero_3d_int(int *** array, int nz, int ny, int nx);
// deallocates and zeros a 3d array of ints of the given dimensions.
// nz is the slowest dimension, nx the fastest.

void output_voltage_array(int count, float *voltages, const char *format, ...);
// outputs voltage array of the given size.

void output_state_array(int count, float *voltages, char *format, ...);
// outputs voltage array of the given size.


#endif /* UTIL_H */

