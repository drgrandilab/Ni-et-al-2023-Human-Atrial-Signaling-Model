/*
 * image_data_writer.cc
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 * Converted to output binary by Timothy D. Butters
 * Updated by Haibo Ni <haibo.ni.0822@gmail.com>
 *
 */
#include <stdio.h>
#include <zlib.h>
#include <endian.h>
#include <stdlib.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <string>
/*#define NX (3)
#define NY (102)
#define NZ (3)*/

// /* New Atria Geometry */
// #define NX (239)
// #define NY (271)
// #define NZ (300)
/* Downsampled IS ventricle dimension */
/*#define NX (459)
#define NY (325)
#define NZ (357)
*/

#define NX (462)
#define NY (325)
#define NZ (358)


#define MIN (-80.0)
#define MAX (20.0)

unsigned char c[NZ][NY][NX];
float *s;
float out_voltage[NZ][NY][NX];
float ReverseFloat(const float inFloat)
{
    // reverse endian of float 32;
    // from Erick
    // from link:
    // http://stackoverflow.com/questions/2782725/converting-float-values-from-big-endian-to-little-endian
    float retVal;
    char *floatToConvert = ( char *) & inFloat;
    char *returnFloat = ( char *) & retVal;
    // swap the bytes into a temporary buffer
    returnFloat[0] = floatToConvert[3];
    returnFloat[1] = floatToConvert[2];
    returnFloat[2] = floatToConvert[1];
    returnFloat[3] = floatToConvert[0];
    return retVal;
}



int main (int argc, char **argv) {
    int g, h, i, offset, n;
    char command[256];
    gzFile gz;
    int count;
    FILE *p, *out;

    if(argc != 3)  {
    	std::cerr << " ------------ input argument number must be 2! -------------" << std::endl;
    	std::cerr << " please run the code as ./vtk input.bin output.vtk" << std::endl<< std::endl;
    	std::exit(0);
    }

    // gz = gzopen("Last_Atria_Geo_Full.geo.gz", "r");
    // gz = gzopen("Vent_seeman13.mat.gz", "r");
    gz = gzopen("Vent_seeman_Add_RV_Stim_RG.mat.gz", "r");
    if(!gz) {
        printf("Geometry file not opened!!!\n");
        exit(0);
    }

    gzread(gz, c, NZ * NY * NX);
    gzclose(gz);

    count = 0;
    for ( g = 0; g < NZ; g++) {
        for ( h = 0; h < NY; h++) {
            for ( i = 0 ; i < NX; i++) {
                if (c[g][h][i] != 0 ) count ++;
            }
        }
    }

    printf(" total num of cells: %d\n", count);


    s = new float [count];
    FILE *in;
    in = fopen(argv[1], "rb");
    fread(s, sizeof(float), count, in);
    fclose(in);

    float max = -80.0;
    float min = 0.0;
    float space = ReverseFloat(-100);
    float S;
    int gg, hh, ii, nn;
    n = 0;
    gg = hh = ii = nn = -1;


    for (g = 0; g < NZ; ++g) {
        for (h = 0; h < NY; ++h) {
            for (i = 0; i < NX; ++i) {
                if (c[g][h][i] == 0) {
                    // fwrite (&space, sizeof(int), 1, out);
                    out_voltage[g][h][i] = space;
                } else {
                    // assert(s[n] == s[n]);  // failure would indicate NaNs
                    if (s[n] != s[n])
                    {
                        printf("NaNs error, program exiting ... ...\n\n");
                        exit(0);
                    }
                    if (s[n] < min)
                    {
                        min = s[n];
                        gg = g;
                        hh = h;
                        ii = i;
                        nn = n;
                    }
                    if (s[n] > max)
                        max = s[n];
                    out_voltage[g][h][i] = ReverseFloat(s[n]);
                    n++;
                }
            }
        }
    }
    std::string str (argv[2]);
    str += ".vtk";
    out = fopen (str.c_str(), "wb");

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "BINARY\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
    fprintf(out, "SPACING 0.33 0.33 0.33\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", NX * NY * NZ);
    fprintf(out, "SCALARS HumanAtrium float 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");
    
    fwrite(out_voltage, sizeof(float), NX * NY * NZ, out);
    fclose(out);



    out = fopen("AP_max_min_summery", "a+");
    printf("gg = %d, hh = %d, ii = %d, nn = %d\n", gg, hh, ii, nn);
    printf("file %s, max = %f, min = %f\n\n", argv[1], max, min );

    fprintf(out, "file %s, max = %f, min = %f\n", argv[1], max, min );
    fclose(out);
    delete [] s;
    return 0;
} /* end of main() */
