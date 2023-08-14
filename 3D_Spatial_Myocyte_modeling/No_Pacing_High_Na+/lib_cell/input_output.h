/* 
 *  input_output.c
 *   
 *  read in .gz .txt .bin file into arrays.
 *  .gz read into  char arrays.
 *  
 *   Haibo Ni <qiangzini@gmail.com>
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>


/* read from gz txt bin file into char/ float arrays, with error handling */
/* output file into gz txt bin file into char/ float arrays, with error handling */
/* 
* 	Haibo Ni 
*   12/12/2014
*   qiangzini@gmail.com
*
*/



#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

char read_char_from_gz(const char *filename, char *array, long int num);
char read_float_from_txt(const char *filename, float *array, long int num);
char read_double_from_txt(const char *filename, double *array, long int num);
char read_float_from_bin(const char *filename, float *array, long int num);
char read_double_from_bin(const char *filename, double *array, long int num);
void print_error_info_file_not_found(const char *filename);
void print_error_info_array_not_malloced(const char *filename);
void output_matrix(const char * filename, const double **matrix, int row, int col);
char output_float_array_bin(const char *filename, const float *voltages, long int num);
char output_double_array_bin(const char *filename, const double *voltages, long int num);
void print_error_info_filename_empty(const char *filename);
char output_double_array_txt(FILE *file, const double *array, long int num);
char output_double_array_txt(const char *filename, const double *array, long int num);
bool OutPutArrayToTxt(std::ofstream &output, const double *array, long int num );
// template <typename T>
// void print_error_info_file_open_failure(T filename) ;

// template <typename T>
// int read_till_end_file(std::string filename, std::vector<T> & data_vec);

// template <typename T>
// int read_num_data_file(std::string filename, T * data_vec, int num);
template <typename T>
void print_error_info_file_open_failure(T filename) {

    fprintf(stderr, "\n\n\n");
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "can not open file %s\n", filename);
    fprintf(stderr, "**************************************\n");
    fprintf(stderr, "assign array to be 0 \n");
    fprintf(stderr, "\n\n\n");
}


template <typename T>
int read_till_end_file(std::string filename, std::vector<T> & data_vec) {

    std::ifstream iFile(filename.c_str());    // input.txt has integers, one per line
// iFile.open("test.dat1", std::ios::in);

    if (iFile.is_open()) {
        int count = 0;
        while (iFile.good()) {
            T x;
            iFile >> x;//data_vec[count];
            if ( iFile.eof() ) break;
            data_vec.push_back(x);
            count++;
            // std::cerr << x << endl;
        }
        return count;
    } else {
        std::cerr << filename << " open failed, source code read_till_end_file" << std::endl;
        std::exit(0);
    }
    return 0;
}

template <typename T>
int read_num_data_file(std::string filename, T * data_vec, int num) {
    std::ifstream inputStream;
    inputStream.open(filename.c_str(), std::ios::in);

    if (inputStream.is_open()) {
        // double myArray[3][5];
        int count = 0;
        for (int i = 0; i < num; i++) {
            if (inputStream.good()) {
                T x;
                inputStream >> x;//;
                if ( inputStream.eof() ) break;
                data_vec[i] = x;
                count ++;

            }
        }
        if(count!= num) {
            std::cerr << " Error in  read_num_data_file " << std::endl;
            std::cerr << " Inconsistent length of data read: actual in >> " << count << " assigned >> " << num << std::endl;
            std::exit(0);
        }
        return count;
    }
    else {
        std::cerr << "****\n";
        std::cerr << "****\n";
        std::cerr << "****\n";
        std::cerr << "Warning!!!! " << filename << " file not opened successfully, maybe file directory invalid!!" << std::endl;
        std::cerr << "Please make sure this does not affect your simulations!!!!" << std::endl;
        std::cerr << "****\n";
        std::cerr << "****\n";
        std::cerr << "****\n";

    }
    return -1;
}
#endif