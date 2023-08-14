#ifndef MATRIX_MEM_HPP
#define MATRIX_MEM_HPP

#include <stdlib.h>

#include <iostream>
#include <cstdlib>



template <class T>
void allocate_and_zero_4d_matrix(T * * **(&A), int NZ, int NY, int NX, int BREADTH);
template <class T>
void deallocate_and_zero_4d_matrix(T * * **(&A), int NZ, int NY, int NX, int BREADTH);
template <class T>
void allocate_and_zero_3d_matrix_pointer(T * * **(&A), int NZ, int NY, int NX);


template <class T>
void allocate_and_zero_2d_matrix_pointer(T * * *(&A),  int NY, int NX);

template <class T>
bool allocate_and_zero_3d_maxtrix(T *** (&array), int nz, int ny, int nx);

template <class T>
bool allocate_and_zero_2d_maxtrix(T ** (&array), int ny, int nx);

template <class T>
bool deallocate_and_zero_3d_matrix(T  *** (& array), int nz, int ny, int nx);

template <class T>
bool deallocate_and_zero_2d_matrix(T  ** (& array), int ny, int nx);
template <class T>
void deallocate_and_zero_2d_matrix_pointer(T * * *(&A), int NY, int NX);

template <class T>
void allocate_and_zero_4d_matrix(T * * **(&A), int NZ, int NY, int NX, int BREADTH) {
	A = new T ***[NZ];
	for (int i = 0; i < NZ; ++i) {
		A[i] = new T **[NY];

		for (int j = 0; j < NY; ++j) {
			A[i][j] = new T*[NX];

			for (int k = 0; k < NX; ++k) {
				A[i][j][k] = new T[BREADTH];
				for (int l = 0; l < BREADTH; ++l) {
					A[i][j][k][l] = 0;
				}
			}
		}
	}
}






template <class T>
void deallocate_and_zero_4d_matrix(T * * **(&A), int NZ, int NY, int NX, int BREADTH) {
	for (int i = 0; i < NZ; ++i) {
		for (int j = 0; j < NY; ++j) {
			for (int k = 0; k < NX; ++k) {
				delete [] A[i][j][k];
			}
			delete [] A[i][j];
		}
		delete [] A[i];
	}
	delete [] A;
}

template <class T>
void allocate_and_zero_3d_matrix_pointer(T * * **(&A), int NZ, int NY, int NX) {
	A = new T ***[NZ];
	for (int i = 0; i < NZ; ++i) {
		A[i] = new T **[NY];

		for (int j = 0; j < NY; ++j) {
			A[i][j] = new T* [NX];
			for (int k = 0; k < NX; ++k) {
				A[i][j][k] = NULL;
			}
		}
	}

}


template <class T>
void allocate_and_zero_2d_matrix_pointer(T ***(&A),  int NY, int NX) {
	A = new T **[NY];
	for (int i = 0; i < NY; ++i) {
		A[i] = new T *[NX];

		for (int j = 0; j < NX; ++j) {
			A[i][j] = NULL;
		}
	}

}


template <class T>
void deallocate_and_zero_3d_matrix_pointer(T * * **(&A), int NZ, int NY, int NX) {
	for (int i = 0; i < NZ; ++i) {
		for (int j = 0; j < NY; ++j) {
			for (int k = 0; k < NX; ++k) {
				delete [] A[i][j][k];
			}
			delete [] A[i][j];
		}
		delete [] A[i];
	}
	delete [] A;
}

template <class T>
void deallocate_and_zero_2d_matrix_pointer(T * * *(&A), int ny, int nx) {

	if (A) {
		for (int j = 0; j < ny; ++j) {
			if (A[j]) {
				for (int i = 0; i < nx; ++i) {
					if (A[j][i]) {
						delete [] A[j][i];
						A[j][i] = NULL;
					}
				}
				delete [] A[j];
				A[j] = NULL;
			}
		}
		delete [] A;
		A = NULL;
	}
}





template <class T>
bool allocate_and_zero_3d_maxtrix(T *** (&array), int nz, int ny, int nx) {
	// allocates and zeros a 3D array of unsigned chars of the given dimensions.
	// nz is the slowest dimension, nx the fastest.
	char error_message[256];
	int k, j, i;

	// allocate the backbone of our array
	// array = (unsigned char***)malloc(nz*sizeof(unsigned char**));
	array = new T **[nz];
	if (!array) {
		sprintf(error_message, "malloc: error allocating z spine!");
		std::cerr << error_message << std::endl;
		std::exit(EXIT_FAILURE);
	}

	for (k = 0; k < nz; ++k) {
		// null it first.
		array[k] = NULL;

		// allocate a pointer to hold the y rows
		// array[k] = (T**)malloc(ny*sizeof(T*));
		array[k] = new T* [ny];
		if (!array[k]) {
			sprintf(error_message, "malloc: error allocating y spur %d", k);
			std::cerr << error_message << std::endl;
			std::exit(EXIT_FAILURE);
		}

		for (j = 0; j < ny; ++j) {
			array[k][j] = NULL;

			// allocate a pointer to hold the y rows
			// array[k][j] = (T*)malloc(nx*sizeof(T));
			array[k][j] = new T [nx];
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
	return true;
}

template <class T>
bool deallocate_and_zero_3d_matrix(T  *** (& array), int nz, int ny, int nx) {
	// deallocates and zeros a 3D array of T s of the given dimensions.
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
	return true;
}



template <class T>
bool allocate_and_zero_2d_maxtrix(T ** (&array), int ny, int nx) {
	// allocates and zeros a 3D array of unsigned chars of the given dimensions.
	// nz is the slowest dimension, nx the fastest.
	char error_message[256];
	int j, i;

	// allocate a pointer to hold the y rows
	// array[k] = (T**)malloc(ny*sizeof(T*));
	array = new T *[ny];
	if (!array) {
		sprintf(error_message, "malloc: error allocating z spine!");
		std::cerr << error_message << std::endl;
		std::exit(EXIT_FAILURE);
	}

	for (j = 0; j < ny; ++j) {
		array[j] = NULL;

		// allocate a pointer to hold the y rows
		// array[k][j] = (T*)malloc(nx*sizeof(T));
		array[j] = new T [nx];
		if (!array[j]) {
			sprintf(error_message, "malloc: error allocating x rows %d", j);
			std::cerr << error_message << std::endl;
			std::exit(EXIT_FAILURE);
		}

		for (i = 0; i < nx; ++i) {
			array[j][i] = 0;
		}
	}
	return true;
}


template <class T>
bool deallocate_and_zero_2d_matrix(T  ** (& array), int ny, int nx) {
	// deallocates and zeros a 3D array of T s of the given dimensions.
	// nz is the slowest dimension, nx the fastest.
	int k, j, i;
	if (array) {
		for (j = 0; j < ny; ++j) {
			if (array[j]) {
				for (i = 0; i < nx; ++i) {
					array[j][i] = 0;
				}
				// free(array[k][j]);
				delete [] array[j];
				array[j] = NULL;
			}
		}
		// free(array[k]);
		delete [] array;
		array = NULL;
	}
	return true;
}



#endif