/*
 * util.c
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 *

Update by Haibo Ni <haibo.ni02@gmail.com>
convert to c++ (new/delete[])
update functions.

 */
#include "util.h"
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <iostream>
#include <cstdlib>



template <class T>
int File_Not_Found(T name) {
    std::cerr << "--------------- Fatal Errorr --------------" << std::endl;
    std::cerr << "\t\tFile " << name << " Note Found!!!" << std::endl;
    std::cerr << "--------------- Fatal Errorr --------------" << std::endl;
    return 0;
}


unsigned char * * *read_and_embed_geometry(const char *filename, int nz, int ny, int nx) {
    // reads an unsigned character geometry from a (possibly gzip'd) file called
    // filename and embeds it in a geometry that is one larger in all directions.
    int j = 0, k = 0, rw = 0;
    gzFile gz;

    unsigned char *** embedded_array;
    allocate_and_zero_3d_maxtrix(embedded_array, nz + 2, ny + 2, nx + 2);

    // open file and test
    gz = gzopen(filename, "r");
    if (!gz) {
        // perror(filename);
        File_Not_Found(filename);
        deallocate_and_zero_3d_matrix(embedded_array, nz + 2, ny + 2, nx + 2);
        std::exit(EXIT_FAILURE);
    }

    // read in the file, one x row at a time.
    for (k = 1; k < nz + 1; k++) {
        for (j = 1; j < ny + 1; j++) {
            rw = gzread(gz, embedded_array[k][j] + 1, nx);
            if (rw != nx) {
                fprintf(stderr, "%d/%d bytes read from '%s' at (%d,%d), exiting!\n",
                        rw, nx, filename, k, j);
                deallocate_and_zero_3d_matrix(embedded_array, nz + 2, ny + 2, nx + 2);
                std::exit(EXIT_FAILURE);
            }
        }
    }

    // close gz file and free memory
    gzclose(gz);

    return embedded_array;
}


unsigned char ** read_and_embed_geometry_2d(const char *filename,  int ny, int nx) {
    // reads an unsigned character geometry from a (possibly gzip'd) file called
    // filename and embeds it in a geometry that is one larger in all directions.
    int j, k, rw;
    gzFile gz;

    unsigned char ** embedded_array;
    allocate_and_zero_2d_maxtrix(embedded_array, ny + 2, nx + 2);

    // open file and test
    gz = gzopen(filename, "r");
    if (!gz) {
        // perror(filename);
        File_Not_Found(filename);
        deallocate_and_zero_2d_matrix(embedded_array,  ny + 2, nx + 2);
        std::exit(EXIT_FAILURE);
    }

    // read in the file, one x row at a time.

    for (j = 1; j < ny + 1; j++) {
        rw = gzread(gz, embedded_array[j] + 1, nx);
        if (rw != nx) {
            fprintf(stderr, "%d/%d bytes read from '%s' at (%d), exiting!\n",
                    rw, nx, filename, j);
            deallocate_and_zero_2d_matrix(embedded_array, ny + 2, nx + 2);
            std::exit(EXIT_FAILURE);
        }
    }


    // close gz file and free memory
    gzclose(gz);

    return embedded_array;
}



float * * *read_and_embed_float_data(const char *filename, int nz, int ny, int nx) {
    // reads an unsigned character geometry from a (possibly gzip'd) file called
    // filename and embeds it in a geometry that is one larger in all directions.
    int j, k, rw;
    gzFile gz;

    float *** embedded_array = allocate_and_zero_3d_float(nz + 2, ny + 2, nx + 2);

    // open file and test
    gz = gzopen(filename, "r");
    if (!gz) {
        perror(filename);
        deallocate_and_zero_3d_float(embedded_array, nz + 2, ny + 2, nx + 2);
        std::exit(EXIT_FAILURE);
    }

    // read in the file, one x row at a time.
    for (k = 1; k < nz + 1; k++) {
        for (j = 1; j < ny + 1; j++) {
            rw = gzread(gz, embedded_array[k][j] + 1, nx * 4);
            if (rw != nx * 4) {
                fprintf(stderr, "%d/%d bytes read from '%s' at (%d,%d), exiting!\n",
                        rw, nx, filename, k, j);
                deallocate_and_zero_3d_float(embedded_array, nz + 2, ny + 2, nx + 2);
                std::exit(EXIT_FAILURE);
            }
        }
    }

    // close gz file and free memory
    gzclose(gz);

    return embedded_array;
}


unsigned char * * *allocate_and_zero_3d_unsigned_char(int nz, int ny, int nx) {
    // allocates and zeros a 3D array of unsigned chars of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    unsigned char *** array;
    char error_message[256];
    int k, j, i;

    // allocate the backbone of our array
    // array = (unsigned char***)malloc(nz*sizeof(unsigned char**));
    array = new unsigned char **[nz];
    if (!array) {
        sprintf(error_message, "malloc: error allocating z spine!");
        std::cerr << error_message << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (k = 0; k < nz; ++k) {
        // null it first.
        array[k] = NULL;

        // allocate a pointer to hold the y rows
        // array[k] = (unsigned char**)malloc(ny*sizeof(unsigned char*));
        array[k] = new unsigned char *[ny];
        if (!array[k]) {
            sprintf(error_message, "malloc: error allocating y spur %d", k);
            std::cerr << error_message << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (j = 0; j < ny; ++j) {
            array[k][j] = NULL;

            // allocate a pointer to hold the y rows
            // array[k][j] = (unsigned char*)malloc(nx*sizeof(unsigned char));
            array[k][j] = new unsigned char [nx];
            if (!array[k][j]) {
                sprintf(error_message, "malloc: error allocating x rows %d, %d", k, j);
                std::cerr << error_message << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (i = 0; i < nx; ++i) {
                array[k][j][i] = 0;
            }
        }
    }

    return array;
}


void deallocate_and_zero_3d_unsigned_char(unsigned char *** array, int nz, int ny, int nx) {
    // deallocates and zeros a 3D array of unsigned chars of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int k, j, i;
    if (array) {
        for (k = 0; k < nz; ++k) {
            if (array[k]) {
                for (j = 0; j < ny; ++j) {
                    if (array[k][j]) {
                        for (i = 0; i < nx; ++i) {
                            array[k][j][i] = 0;
                        }
                        delete [] array[k][j];
                        array[k][j] = NULL;
                    }
                }
                delete [] array[k];
                array[k] = NULL;
            }
        }
        delete [] array;
        array = NULL;
    }
}



double * * *allocate_and_zero_3d_double(int nz, int ny, int nx) {
    // allocates and zeros a 3D array of doubles of the given dimensions.
    // nz is the slowest dimension, nx the fastest.

    double *** array;
    char error_message[256];
    int k, j, i;

    // allocate the backbone of our array
    // array =(double ***) malloc(nz*sizeof(double**));
    array = new double **[nz];
    if (!array) {
        sprintf(error_message, "malloc: error allocating z spine!");
        std::cerr << error_message << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (k = 0; k < nz; ++k) {
        // null it first.
        array[k] = NULL;

        // allocate a pointer to hold the y rows
        // array[k] = (double **) malloc(ny*sizeof(double*));
        array[k] = new double* [ny];
        if (!array[k]) {
            sprintf(error_message, "malloc: error allocating y spur %d", k);
            std::cerr << error_message << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (j = 0; j < ny; ++j) {
            array[k][j] = NULL;

            // allocate a pointer to hold the y rows
            // array[k][j] = (double *)malloc(nx * sizeof(double));
            array[k][j] = new double [nx];
            if (!array[k][j]) {
                sprintf(error_message, "malloc: error allocating x rows %d, %d", k, j);
                std::cerr << error_message << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (i = 0; i < nx; ++i) {
                array[k][j][i] = 0;
            }
        }
    }
    return array;
}
void deallocate_and_zero_3d_double(double *** array, int nz, int ny, int nx) {
    // deallocates and zeros a 3D array of floats of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int k, j, i;
    if (array) {
        for (k = 0; k < nz; ++k) {
            if (array[k]) {
                for (j = 0; j < ny; ++j) {
                    if (array[k][j]) {
                        for (i = 0; i < nx; ++i) {
                            array[k][j][i] = 0;
                        }
                        delete [] array[k][j];
                        // array[k][j] = NULL;
                    }
                }
                delete [] array[k];
                // array[k] = NULL;
            }
        }
        delete [] array[k];
        // array = NULL;
    }
}


float * * *allocate_and_zero_3d_float(int nz, int ny, int nx) {
    // allocates and zeros a 3D array of floats of the given dimensions.
    // nz is the slowest dimension, nx the fastest.

    float *** array;
    char error_message[256];
    int k, j, i;

    // allocate the backbone of our array
    // array = (float ** *)malloc(nz * sizeof(float **));
    array = new float **[nz];
    if (!array) {
        sprintf(error_message, "malloc: error allocating z spine!");
        std::cerr << error_message << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (k = 0; k < nz; ++k) {
        // null it first.
        array[k] = NULL;

        // allocate a pointer to hold the y rows
        // array[k] = (float **)malloc(ny * sizeof(float *));
        array[k] = new float * [ny];
        if (!array[k]) {
            sprintf(error_message, "malloc: error allocating y spur %d", k);
            std::cerr << error_message << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (j = 0; j < ny; ++j) {
            array[k][j] = NULL;

            // allocate a pointer to hold the y rows
            // array[k][j] = (float *)malloc(nx * sizeof(float));
            array[k][j] = new float [nx];
            if (!array[k][j]) {
                sprintf(error_message, "malloc: error allocating x rows %d, %d", k, j);
                std::cerr << error_message << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (i = 0; i < nx; ++i) {
                array[k][j][i] = 0;
            }
        }
    }

    return array;
}

void deallocate_and_zero_3d_float(float *** array, int nz, int ny, int nx) {
    // deallocates and zeros a 3D array of floats of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int k, j, i;
    if (array) {
        for (k = 0; k < nz; ++k) {
            if (array[k]) {
                for (j = 0; j < ny; ++j) {
                    if (array[k][j]) {
                        for (i = 0; i < nx; ++i) {
                            array[k][j][i] = 0;
                        }
                        // free(array[k][j]);
                        delete [] array[k][j];
                        array[k][j] = NULL;
                    }
                }
                // free(array[k]);
                delete [] array[k];
                array[k] = NULL;
            }
        }
        // free(array);
        delete [] array;
        array = NULL;
    }
}

double * * **allocate_and_zero_3d_double_pointer(int nz, int ny, int nx) {
    // allocates and zeros a 3D array of pointers to doubles of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    double * * **array;
    char error_message[256];
    int k, j, i;

    // allocate the backbone of our array
    // array = (double ** **)malloc(nz * sizeof(double ** *));
    array = new double *** [nz];
    if (!array) {
        sprintf(error_message, "malloc: error allocating z spine!");
        std::cerr << error_message << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (k = 0; k < nz; ++k) {
        // null it first.
        array[k] = NULL;

        // allocate a pointer to hold the y rows
        // array[k] = (double ** *)malloc(ny * sizeof(double **));
        array[k] = new double **[ny];
        if (!array[k]) {
            sprintf(error_message, "malloc: error allocating y spur %d", k);
            std::cerr << error_message << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (j = 0; j < ny; ++j) {
            array[k][j] = NULL;

            // allocate a pointer to hold the y rows
            // array[k][j] = (double **)malloc(nx * sizeof(double *));
            array[k][j] = new double* [nx];
            if (!array[k][j]) {
                sprintf(error_message, "malloc: error allocating x rows %d, %d", k, j);
                std::cerr << error_message << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (i = 0; i < nx; ++i) {
                array[k][j][i] = NULL;
            }
        }
    }

    return array;
}

void deallocate_and_zero_3d_double_pointer(double * * **array, int nz, int ny, int nx) {
    // deallocates and zeros a 3D array of floats of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int k, j, i;
    if (array) {
        for (k = 0; k < nz; ++k) {
            if (array[k]) {
                for (j = 0; j < ny; ++j) {
                    if (array[k][j]) {
                        for (i = 0; i < nx; ++i) {
                            if (array[k][j][i]) {
                                delete [] array[k][j][i];
                                array[k][j][i] = NULL;
                            }
                        }
                        delete [] array[k][j];
                        array[k][j] = NULL;
                    }
                }
                delete [] array[k];
                array[k] = NULL;
            }
        }
        delete [] array;
        array = NULL;
    }

}

int * * *allocate_and_zero_3d_int(int nz, int ny, int nx) {
    // allocates and zeros a 3d array of ints of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    return allocate_and_set_3d_int(nz, ny, nx, 0);
}

int * * *allocate_and_set_3d_int(int nz, int ny, int nx, int set) {
    // allocates and sets a 3d array of ints of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int *** array;
    char error_message[256];
    int k, j, i;

    // allocate the backbone of our array
    // array = (int ** *) malloc(nz * sizeof(int **));
    array = new int **[nz];
    if (!array) {
        sprintf(error_message, "malloc: error allocating z spine!");
        std::cerr << error_message << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (k = 0; k < nz; ++k) {
        // null it first.
        array[k] = NULL;

        // allocate a pointer to hold the y rows
        // array[k] = (int **) malloc(ny * sizeof(int *));
        array[k] = new int *[ny];
        if (!array[k]) {
            sprintf(error_message, "malloc: error allocating y spur %d", k);
            std::cerr << error_message << std::endl;
            std::exit(EXIT_FAILURE);
        }

        for (j = 0; j < ny; ++j) {
            array[k][j] = NULL;

            // allocate a pointer to hold the y rows
            // array[k][j] = (int *) malloc(nx * sizeof(int));
            array[k][j] = new int [nx];
            if (!array[k][j]) {
                sprintf(error_message, "malloc: error allocating x rows %d, %d", k, j);
                std::cerr << error_message << std::endl;
                std::exit(EXIT_FAILURE);
            }

            for (i = 0; i < nx; ++i) {
                array[k][j][i] = set;
            }
        }
    }

    return array;
}

void deallocate_and_zero_3d_int(int *** array, int nz, int ny, int nx) {
    // deallocates and zeros a 3d array of ints of the given dimensions.
    // nz is the slowest dimension, nx the fastest.
    int k, j, i;
    if (array) {
        for (k = 0; k < nz; ++k) {
            if (array[k]) {
                for (j = 0; j < ny; ++j) {
                    if (array[k][j]) {
                        for (i = 0; i < nx; ++i) {
                            array[k][j][i] = 0;
                        }
                        // free(array[k][j]);
                        delete [] array[k][j];
                        array[k][j] = NULL;
                    }
                }
                // free(array[k]);
                delete [] array[k];
                array[k] = NULL;
            }
        }
        // free(array);

        array = NULL;
    }
}

void output_voltage_array(int count, float *voltages, const char *format, ...) {
    char filename[1024];
    va_list args;
    FILE *file;
    int rw;


    // start the args list
    va_start(args, format);

    // print to the filename
    vsnprintf(filename, 1024, format, args);

    // open the file
    file = fopen(filename, "wb");
    if (!file) {
        perror(filename);
        std::exit(EXIT_FAILURE);
    }

    // output
    rw = fwrite(voltages, sizeof(float), count, file);
    if (rw != count) {
        fprintf(stderr, "%d/%d floats written to '%s', exiting!\n", rw, count, filename);
        std::exit(EXIT_FAILURE);
    }

    fclose(file);
}

void output_state_array(int count, float *voltages, char *format, ...) {
    char filename[256];
    va_list args;
    FILE *file;
    int rw;


    // start the args list
    va_start(args, format);

    // print to the filename
    vsnprintf(filename, 256, format, args);

    // open the file
    file = fopen("State.bin", "wb");
    if (!file) {
        perror(filename);
        std::exit(EXIT_FAILURE);
    }

    // output
    rw = fwrite(voltages, sizeof(float), count, file);
    if (rw != count) {
        fprintf(stderr, "%d/%d floats written to '%s', exiting!\n", rw, count, filename);
        std::exit(EXIT_FAILURE);
    }

    fclose(file);
}

