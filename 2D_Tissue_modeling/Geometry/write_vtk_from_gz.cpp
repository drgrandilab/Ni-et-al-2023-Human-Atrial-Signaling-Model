/*
 * image_data_writer.cc
 * Jonathan D. Stott <jonathan.stott@gmail.com>
 * Converted to output binary by Timothy D. Butters
 *
 */
#include <stdio.h>
#include <zlib.h>
#include <endian.h>
#include <stdlib.h>

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

// #define NX 325//462
// #define NY 325 //325
// #define NZ 425 //358



#define NX (125)
#define NY (120)
#define NZ (1)

#define MIN (-80.0)
#define MAX (20.0)

unsigned char c[NZ][NY][NX];
float *s;
int out_voltage[NZ][NY][NX];



int main (int argc, char **argv) {
    int g, h, i, offset, n;
    char command[256];
    gzFile gz;
    int count;
    FILE *p, *out;
    char *str;

    str = (char *) malloc(100 * sizeof(char));

    out = fopen ("vtk_from_gz.vtk", "wt");

    fprintf(out, "# vtk DataFile Version 3.0\n");
    fprintf(out, "vtk output\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET STRUCTURED_POINTS\n");
    fprintf(out, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
    fprintf(out, "SPACING 1 1 1\n");
    fprintf(out, "ORIGIN 0 0 0\n");
    fprintf(out, "POINT_DATA %d\n", NX * NY * NZ);
    fprintf(out, "SCALARS HumanAtrium int 1\n");
    fprintf(out, "LOOKUP_TABLE default\n");

    // gz = gzopen("Last_Atria_Geo_Full.geo.gz", "r");
    gz = gzopen(argv[1], "r");
    gzread(gz, c, NZ * NY * NX);
    gzclose(gz);


    int max = -90;
    int min = -90;
    int space = htobe32(0);
    int S;
    int gg, hh, ii, nn;
    n = 0;
    for (g = 0; g < NZ; ++g) {
        for (h = 0; h < NY; ++h) {
            for (i = 0; i < NX; ++i) {
                fprintf(out, "%d ", (int) (c[g][h][i]));
            }
            fprintf(out, "\n");
        }
    }

    fclose(out);

    free(s);
    free(str);
    return 0;
} /* end of main() */
