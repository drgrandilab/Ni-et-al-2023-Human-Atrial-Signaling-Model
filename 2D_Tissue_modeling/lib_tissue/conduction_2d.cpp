#include "conduction_2d.h"



double **generate_laplacian_heterogeneity_2d (unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double d2) {

	double ** M_dd; // = (double ** *)allocate_and_zero_3d_double( nz + 2,  ny + 2,  nx + 2);
	double ** M_d2; // = (double ** *)allocate_and_zero_3d_double( nz + 2,  ny + 2,  nx + 2);

	allocate_and_zero_2d_maxtrix(M_dd, ny + 2,  nx + 2);
	allocate_and_zero_2d_maxtrix(M_d2, ny + 2,  nx + 2);

	int i, j;

	for (j = 1; j < ny + 1; ++j)
	{

		for (i = 1; i < nx + 1; ++i)
		{
			M_dd[j][i] = dd;
			M_d2[j][i] = d2;

			/* 10 percent decrease in diffusion coefficient from ENDO to EPI*/
			// if (atrium[k][j][i] <= 10 && atrium[k][j][i] >= 1)
			// {
			//     M_dd[k][j][i] = dd * ( 1.0 - 0.1 * ((double)atrium[k][j][i] - 1.0) / 9.0);
			//     M_d2[k][j][i] = d2 * ( 1.0 - 0.1 * ((double)atrium[k][j][i] - 1.0) / 9.0);
			// }
			// if (atrium[k][j][i] <= 20 && atrium[k][j][i] >= 11)
			// {
			//     M_dd[k][j][i] = dd * ( 1.0 - 0.1 * ((double)atrium[k][j][i] - 11.0) / 9.0);
			//     M_d2[k][j][i] = d2 * ( 1.0 - 0.1 * ((double)atrium[k][j][i] - 11.0) / 9.0);
			// }
		}
	}



	double **lap;

	lap = generate_laplacian_with_diffusion_heterogeneity_2d(atrium,
	        ny, nx, count,
	        phi_name,
	        dx, dy,
	        M_dd, M_d2);

	deallocate_and_zero_2d_matrix(M_dd,  ny + 2,  nx + 2);
	deallocate_and_zero_2d_matrix(M_d2,  ny + 2,  nx + 2);
	return lap;
}



double **generate_laplacian_heterogeneity_2d_test(unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double* d2) {

	double ** M_dd; // = (double ** *)allocate_and_zero_3d_double( nz + 2,  ny + 2,  nx + 2);
	double ** M_d2; // = (double ** *)allocate_and_zero_3d_double( nz + 2,  ny + 2,  nx + 2);

	allocate_and_zero_2d_maxtrix(M_dd, ny + 2,  nx + 2);
	allocate_and_zero_2d_maxtrix(M_d2, ny + 2,  nx + 2);

	int i, j;
	int counter = 0;

	for (j = 1; j < ny + 1; ++j)
	{

		for (i = 1; i < nx + 1; ++i)
		{
			M_dd[j][i] = dd;
			M_d2[j][i] = d2[counter];

			counter ++;
		}
	}

	if(counter != count) {
		std::cerr << " counter != count in generate_laplacian_heterogeneity_2d_test" << std::endl;
		std::exit(0);
	}

	double **lap;

	lap = generate_laplacian_with_diffusion_heterogeneity_2d(atrium,
	        ny, nx, count,
	        phi_name,
	        dx, dy,
	        M_dd, M_d2);

	deallocate_and_zero_2d_matrix(M_dd,  ny + 2,  nx + 2);
	deallocate_and_zero_2d_matrix(M_d2,  ny + 2,  nx + 2);
	return lap;
}




double **generate_laplacian_with_diffusion_heterogeneity_2d(unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy,
        double **dd, double **d2) {
	// generates the laplacian operator as a linear array.
	double **laplacian;
	double * * *lap;
	int i, j, n, c;

	if (!dd) {
		perror("dd should a matrix");
		exit(EXIT_FAILURE);
	}


	if (!d2) {
		perror("d2 should a matrix");
		exit(EXIT_FAILURE);
	}

	// first, generate your laplacian.
	// this is the old version, diffusion was universally the same
	// lap = generate_3d_laplacian(atrium, nz, ny, nx, theta_name, phi_name, dx, dd, d2);

	// and the new version:
	lap = generate_2d_laplacian_with_diffusion_heterogeneity(atrium,
	        ny, nx,
	        phi_name,
	        dx, dy,
	        dd, d2);

	// allocate the 2d array.
	// laplacian = allocate_2d_laplacian(count);
	allocate_and_zero_2d_maxtrix(laplacian, count, 9);

	c = 0;

	for (j = 1; j < ny + 1; ++j) {
		for (i = 1; i < nx + 1; ++i) {
			if (lap[j][i] != NULL) {
				for (n = 0; n < 9; ++n) {
					laplacian[c][n] = lap[j][i][n];
				}
				c++;
			}
		}
	}


	deallocate_and_zero_2d_matrix_pointer(lap, ny + 2, nx + 2);

	return laplacian;
}



void calculate_fibre_unit_vectors_2d(unsigned char ** phi,
                                     int ny, int nx, float **y, float **x) {
	int i, j;
	float t, p;


	for (j = 1; j < ny + 1; ++j) {
		for (i = 1; i < nx + 1; ++i) {
			t = (float) ((phi[j][i] / 255.0) * M_PI);
			// p = (float) ((phi[k][j][i] / 254.0) * M_PI);

			x[j][i] = (float) (cos(t));
			y[j][i] = (float) (sin(t));
			// z[j][i] = (float) cos(t);
		} /*else {
				x[k][j][i] = 0.0;
				y[k][j][i] = 0.0;
			}*/
	}
}



void write_fibre_vtk_2d(float ** y, float ** x, int ny, int nx, double dx) {


	FILE *out;
	out = fopen ("atrium_fibre.vtk", "wt");

	fprintf (out, "# vtk DataFile Version 3.0\n");
	fprintf (out, "vtk output\n");
	fprintf (out, "ASCII\n");
	fprintf (out, "DATASET STRUCTURED_POINTS\n");
	fprintf (out, "DIMENSIONS %d %d %d\n", nx, ny, 1);
	fprintf (out, "SPACING %f %f %f \n", dx, dx, dx);
	fprintf (out, "ORIGIN 0 0 0\n");
	fprintf (out, "POINT_DATA %d\n", (ny ) * (nx));
	fprintf (out, "SCALARS HumanAtrium float 3\n");
	fprintf (out, "LOOKUP_TABLE default\n");

	int i, j;
	float t, p;

	// for (k = 1; k < nz + 1; ++k) {
	for (j = 1; j < ny + 1; ++j) {
		for (i = 1; i < nx + 1; ++i) {

			fprintf(out, "%f %f 0.0 ", x[j][i], y[j][i]);

		}
		fprintf(out, "\n");
	}
	// }
	fclose(out);

}





/* compute the diffusion tensor */
void calculate_diffusion_tensor_components_with_diffusion_heterogeneity_2d(unsigned char ** atrium,
        int ny, int nx,
        float **y, float **x,
        double dx, double dy,
        double **dd, double **d2,
        double  * **dc, double  * **df,
        double  * **d_dd, double  * **d_d2) {
	int i, j;

	for (j = 1; j < ny + 1; ++j) {
		for (i = 1; i < nx + 1; ++i) {
			if (atrium[j][i] > 0) {
				// allocate diffusion tensor at this point
				dc[j][i] = new double [4]; // (double *)  malloc(4 * sizeof(double));
				if (!dc[j][i]) {
					fprintf(stderr, "%d, %d failed to allocate!\n", j, i);
				}
				df[j][i] = new double [4]; // (double *) malloc(4 * sizeof(double));
				if (!df[j][i]) {
					fprintf(stderr, "%d, %d failed to allocate!\n", j, i);
				}

				d_dd[j][i] = new double [2]; // (double *) malloc(3 * sizeof(double));
				if (!d_dd[j][i]) {
					fprintf(stderr, "d_dd %d, %d failed to allocate!\n", j, i);
				}

				d_d2[j][i] = new double [2]; // (double *) malloc(3 * sizeof(double));
				if (!d_d2[j][i]) {
					fprintf(stderr, "d_d2 %d, %d failed to allocate!\n", j, i);
				}


				// calculate dc diffusion tensor projected with Fibre Orientation, the material is axially anisotropic
				dc[j][i][0] = d2[j][i] + (dd[j][i] * x[j][i] * x[j][i]);
				dc[j][i][1] = dd[j][i] * x[j][i] * y[j][i];
				dc[j][i][2] = dd[j][i] * y[j][i] * x[j][i];
				dc[j][i][3] = d2[j][i] + (dd[j][i] * y[j][i] * y[j][i]);
				// dc[j][i][5] = dd[j][i] * y[j][i] * z[j][i];
				// dc[j][i][6] = dd[j][i] * z[j][i] * x[j][i];
				// dc[j][i][7] = dd[j][i] * z[j][i] * y[j][i];
				// dc[j][i][8] = d2[j][i] + (dd[j][i] * z[j][i] * z[j][i]);

				// calculate df
				// dfx/dx
				df[j][i][0] = (x[j][i + 1] - x[j][i - 1]) / (2 * dx);
				if (atrium[j][i + 1] == 0)
					df[j][i][0] = (x[j][i] - x[j][i - 1]) / (dx);
				if (atrium[j][i - 1] == 0)
					df[j][i][0] = (x[j][i + 1] - x[j][i]) / (dx);
				if (atrium[j][i + 1] == 0 && atrium[j][i + 1] == 0)
					df[j][i][0] = 0.0;

				//dfx/dy
				df[j][i][1] = (x[j + 1][i] - x[j - 1][i]) / (2 * dy);
				if (atrium[j + 1][i] == 0)
					df[j][i][1] = (x[j][i] - x[j - 1][i]) / (dy);
				if (atrium[j - 1][i] == 0)
					df[j][i][1] = (x[j + 1][i] - x[j][i]) / (dy);
				if (atrium[j + 1][i] == 0 && atrium[j - 1][i] == 0)
					df[j][i][1] = 0.0;

				// dfx/dz
				/*df[k][j][i][2] = (x[k + 1][j][i] - x[k - 1][j][i]) / (2 * dz);
				if (atrium[k + 1][j][i] == 0)
					df[k][j][i][2] = (x[k][j][i] - x[k - 1][j][i]) / (dz);
				if (atrium[k - 1][j][i] == 0)
					df[k][j][i][2] = (x[k + 1][j][i] - x[k][j][i]) / (dz);
				if (atrium[k + 1][j][i] == 0 && atrium[k - 1][j][i] == 0)
					df[k][j][i][2] = 0.0;*/


				// dfy/dx
				df[j][i][2] = (y[j][i + 1] - y[j][i - 1]) / (2 * dx);
				if (atrium[j][i + 1] == 0)
					df[j][i][2] = (y[j][i] - y[j][i - 1]) / (dx);
				else if (atrium[j][i - 1] == 0)
					df[j][i][2] = (y[j][i + 1] - y[j][i]) / (dx);
				else if (atrium[j][i + 1] == 0 && atrium[j][i + 1] == 0)
					df[j][i][2] = 0.0;

				//dfy/dy
				df[j][i][3] = (y[j + 1][i] - y[j - 1][i]) / (2 * dy);
				if (atrium[j + 1][i] == 0)
					df[j][i][3] = (y[j][i] - y[j - 1][i]) / (dy);
				else if (atrium[j - 1][i] == 0)
					df[j][i][3] = (y[j + 1][i] - y[j][i]) / (dy);
				else if (atrium[j + 1][i] == 0 && atrium[j - 1][i] == 0)
					df[j][i][3] = 0.0;

				// dfy/dz
				/*df[k][j][i][5] = (y[k + 1][j][i] - y[k - 1][j][i]) / (2 * dz);
				if (atrium[k + 1][j][i] == 0)
					df[k][j][i][5] = (y[k][j][i] - y[k - 1][j][i]) / (dz);
				if (atrium[k - 1][j][i] == 0)
					df[k][j][i][5] = (y[k + 1][j][i] - y[k][j][i]) / (dz);
				if (atrium[k + 1][j][i] == 0 && atrium[k - 1][j][i] == 0)
					df[k][j][i][5] = 0.0;*/

				// dfz/dx
				// df[j][i][6] = (z[j][i + 1] - z[j][i - 1]) / (2 * dx);
				// if (atrium[j][i + 1] == 0)
				// 	df[j][i][6] = (z[j][i] - z[j][i + 1]) / (dx);
				// if (atrium[j][i - 1] == 0)
				// 	df[j][i][6] = (z[j][i + 1] - z[j][i]) / (dx);
				// if (atrium[j][i + 1] == 0 && atrium[j][i - 1] == 0) // + - changed
				// 	df[j][i][6] = 0.0;

				//dfz/dy
				// df[j][i][7] = (z[j + 1][i] - z[j - 1][i]) / (2 * dy);
				// if (atrium[j + 1][i] == 0)
				// 	df[j][i][7] = (z[j][i] - z[j - 1][i]) / (dy);
				// if (atrium[j - 1][i] == 0)
				// 	df[j][i][7] = (z[j + 1][i] - z[j][i]) / (dy);
				// if (atrium[j + 1][i] == 0 && atrium[j - 1][i] == 0)
				// 	df[j][i][7] = 0.0;

				// dfz/dz
				/*df[j][i][8] = (z[k + 1][j][i] - z[k - 1][j][i]) / (2 * dz);
				if (atrium[k + 1][j][i] == 0)
					df[j][i][8] = (z[j][i] - z[k - 1][j][i]) / (dz);
				if (atrium[k - 1][j][i] == 0)
					df[j][i][8] = (z[k + 1][j][i] - z[j][i]) / (dz);
				if (atrium[k + 1][j][i] == 0 && atrium[k - 1][j][i] == 0)
					df[j][i][8] = 0.0;*/

				// calculate d_dd
				// d_dd/dx
				d_dd[j][i][0] = (dd[j][i + 1] - dd[j][i - 1]) / (2 * dx);
				if (atrium[j][i + 1] == 0)
					d_dd[j][i][0] = (dd[j][i] - dd[j][i - 1]) / (dx);
				else if (atrium[j][i - 1] == 0)
					d_dd[j][i][0] = (dd[j][i + 1] - dd[j][i]) / (dx);
				else if (atrium[j][i + 1] == 0 && atrium[j][i + 1] == 0)
					d_dd[j][i][0] = 0.0;

				// d_dd/dy
				d_dd[j][i][1] = (dd[j + 1][i] - dd[j - 1][i]) / (2 * dy);
				if (atrium[j + 1][i] == 0)
					d_dd[j][i][1] = (dd[j][i] - dd[j - 1][i]) / (dy);
				else if (atrium[j - 1][i] == 0)
					d_dd[j][i][1] = (dd[j + 1][i] - dd[j][i]) / (dy);
				else if (atrium[j + 1][i] == 0 && atrium[j - 1][i] == 0)
					d_dd[j][i][1] = 0.0;

				// d_dd/dz
				/*d_dd[k][j][i][2] = (dd[k + 1][j][i] - dd[k - 1][j][i]) / (2 * dz);
				if (atrium[k + 1][j][i] == 0)
					d_dd[k][j][i][2] = (dd[k][j][i] - dd[k - 1][j][i]) / (dz);
				if (atrium[k - 1][j][i] == 0)
					d_dd[k][j][i][2] = (dd[k + 1][j][i] - dd[k][j][i]) / (dz);
				if (atrium[k + 1][j][i] == 0 && atrium[k - 1][j][i] == 0)
					d_dd[k][j][i][2] = 0.0;*/

				// calculate d_d2
				// d_d2/dx
				d_d2[j][i][0] = (d2[j][i + 1] - d2[j][i - 1]) / (2 * dx);
				if (atrium[j][i + 1] == 0)
					d_d2[j][i][0] = (d2[j][i] - d2[j][i - 1]) / (dx);
				else if (atrium[j][i - 1] == 0)
					d_d2[j][i][0] = (d2[j][i + 1] - d2[j][i]) / (dx);
				else if (atrium[j][i + 1] == 0 && atrium[j][i + 1] == 0)
					d_d2[j][i][0] = 0.0;

				// d_d2/dy
				d_d2[j][i][1] = (d2[j + 1][i] - d2[j - 1][i]) / (2 * dy);
				if (atrium[j + 1][i] == 0)
					d_d2[j][i][1] = (d2[j][i] - d2[j - 1][i]) / (dy);
				else if (atrium[j - 1][i] == 0)
					d_d2[j][i][1] = (d2[j + 1][i] - d2[j][i]) / (dy);
				else if (atrium[j + 1][i] == 0 && atrium[j - 1][i] == 0)
					d_d2[j][i][1] = 0.0;


				// if(d_d2[j][i][1] !=0 or d_d2[j][i][0] != 0 or d_dd[j][i][1] !=0 or d_dd[j][i][0] != 0)
				/*if (df[j][i][3] !=0 or df[j][i][2] !=0 or df[j][i][1] !=0 or df[j][i][3] !=0)
				 {
					std::cout << "EEE" << std::endl;
				}
				*/
				// d_d2/dz
				/*d_d2[k][j][i][2] = (d2[k + 1][j][i] - d2[k - 1][j][i]) / (2 * dz);
				if (atrium[k + 1][j][i] == 0)
					d_d2[k][j][i][2] = (d2[k][j][i] - d2[k - 1][j][i]) / (dz);
				if (atrium[k - 1][j][i] == 0)
					d_d2[k][j][i][2] = (d2[k + 1][j][i] - d2[k][j][i]) / (dz);
				if (atrium[k + 1][j][i] == 0 && atrium[k - 1][j][i] == 0)
					d_d2[k][j][i][2] = 0.0;*/


				/*if (d_d2[k][j][i][2] != 0.0 || d_d2[k][j][i][1] != 0.0 || d_d2[k][j][i][0] != 0.0)
				    printf("k = %d, j = %d, i = %d\n", k, j, i);
				if (d_dd[k][j][i][2] != 0.0 || d_dd[k][j][i][1] != 0.0 || d_dd[k][j][i][0] != 0.0)
				    printf("k = %d, j = %d, i = %d\n", k, j, i);*/
			}
		}
	}
}


double * * *generate_2d_laplacian_with_diffusion_heterogeneity(unsigned char ** atrium,
        int ny, int nx,
        const char *phi_name,
        double dx, double dy,
        double **dd, double **d2) {

	unsigned char ** phi  = read_and_embed_geometry_2d(phi_name, ny, nx);
	// unsigned char ** phi    = read_and_embed_geometry_2d(phi_name,   ny, nx);

	float **x; //  = allocate_and_zero_3d_float(nz + 2, ny + 2, nx + 2);
	float **y; //  = allocate_and_zero_3d_float(nz + 2, ny + 2, nx + 2);
	// float **z; //  = allocate_and_zero_3d_float(nz + 2, ny + 2, nx + 2);
	allocate_and_zero_2d_maxtrix( x, ny + 2, nx + 2);
	allocate_and_zero_2d_maxtrix( y, ny + 2, nx + 2);
	// allocate_and_zero_2d_maxtrix( z, ny, nx);

	// calculate the unit vectors from  the angles.
	calculate_fibre_unit_vectors_2d(phi, ny, nx, y, x);

	write_fibre_vtk_2d(y, x, ny, nx, dx);


	deallocate_and_zero_2d_matrix(phi,  ny + 2, nx + 2);
	// deallocate_and_zero_3d_unsigned_char(phi,    nz + 2, ny + 2, nx + 2);

	double * * *dc  ; // = allocate_and_zero_2d_matrix_pointer( ny + 2, nx + 2);
	double * * *df ; //  = allocate_and_zero_2d_matrix_pointer( ny + 2, nx + 2);
	double * * *d_dd  ; // = allocate_and_zero_2d_matrix_pointer( ny + 2, nx + 2);
	double * * *d_d2  ; // = allocate_and_zero_2d_matrix_pointer( ny + 2, nx + 2);

	allocate_and_zero_2d_matrix_pointer(dc, ny + 2, nx + 2);
	allocate_and_zero_2d_matrix_pointer(df, ny + 2, nx + 2);
	allocate_and_zero_2d_matrix_pointer(d_dd, ny + 2, nx + 2);
	allocate_and_zero_2d_matrix_pointer(d_d2, ny + 2, nx + 2);


	calculate_diffusion_tensor_components_with_diffusion_heterogeneity_2d(atrium,
	        ny, nx,
	        y, x,
	        dx, dy,
	        dd, d2,
	        dc, df,
	        d_dd, d_d2);

	double *** lap;
	allocate_and_zero_2d_matrix_pointer(lap, ny + 2, nx + 2);
	/*calculate_laplacian_components(atrium,
	                               nz, ny, nx,
	                               z, y, x,
	                               dx, dd,
	                               dc, df, lap);*/

	calculate_laplacian_components_with_diffusion_heterogeneity_2d (atrium,
	        ny, nx, \
	        y, x, \
	        dx, dy, \
	        dd, d2, \
	        dc, df, \
	        d_dd, d_d2, \
	        lap);

	// deallocate x, y, z
	deallocate_and_zero_2d_matrix(x, ny + 2, nx + 2);
	deallocate_and_zero_2d_matrix(y, ny + 2, nx + 2);
	// deallocate_and_zero_2d_matrix(z, ny + 2, nx + 2);

	/* deallocate the tempoery memery */
	deallocate_and_zero_2d_matrix_pointer(dc,  ny + 2, nx + 2);
	deallocate_and_zero_2d_matrix_pointer(df,  ny + 2, nx + 2);
	deallocate_and_zero_2d_matrix_pointer(d_dd,  ny + 2, nx + 2);
	deallocate_and_zero_2d_matrix_pointer(d_d2,  ny + 2, nx + 2);

	return lap;
}



int ** generate_numeric_map_2d(unsigned char ** atrium, int ny, int nx, int * count) {
// generates a numeric map of the tissue, with all tissue nodes assigned a
// number, and all non-tissue nodes assigned -1.  Also returns a count of the
// cells.
	int ** map;
	allocate_and_zero_2d_maxtrix(map, ny + 2, nx + 2);

	int i, j;
	int c = 0;

	// important, map id starts from 0, to index Vm arrays. (c is zero indexed)

	// 18:45:55, Wed, 16-December-2020, By Haibo
	// found a bug here. Should loop from 0 to nx/ny + 2; otherwise map id of two-layer grids outside the
	// tissue would be assigned to zero, which would cause all kinds of problems
	for (j = 0; j < ny + 2; j++)
		for (i = 0; i < nx + 2; i++) {

			map[j][i] =  - 1;  // initialise map to -1;
			if (atrium[j][i] > 0) {
				map[j][i] = c;
				c++;
			}
		}

	*count = c;
	return map;
}





int ** allocate_numeric_map_2d(int count) {
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
		nbd[i] = new int [9];
		if (!nbd[i]) {
			perror("Neighbours array malloc");
			exit(EXIT_FAILURE);
		}

		for (j = 0; j < 9; ++j) {
			nbd[i][j] = -1;
		}
	}

	return nbd;
}




int ** generate_neighbours_from_numeric_map_2d(int ** map, int ny, int nx, int count) {
	// Generates a neighbours map from the given numeric map.
	int ** nbd;
	int n, i, j, k;

	// allocate the map.
	nbd = allocate_numeric_map_2d(count);

	// 13:30:21, Tue, 31-December-2019, By Haibo
	// iterate over the map, assigning neighbours.
	// self         = 0
	// xplus        = 1
	// xminus       = 2
	// yplus        = 3
	// yminus       = 4
	// xplusyplus   = 7 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
	// xminusyminus = 8 -2  // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
	// xplusyminus  = 9 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
	// xminusyplus  = 10 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo



	n = 0;
	for (j = 1; j < ny + 1; j++)
		for (i = 1; i < nx + 1; i++) {
			if (map[j][i] > -1) {
				nbd[n][ 0] = map[j][i];

				// simple directions.
				// all get assigned the real answer, or self.
				nbd[n][ 1] = (map[j][i + 1] > -1) ? map[j][i + 1] : map[j][i];
				nbd[n][ 2] = (map[j][i - 1] > -1) ? map[j][i - 1] : map[j][i];
				nbd[n][ 3] = (map[j + 1][i] > -1) ? map[j + 1][i] : map[j][i];
				nbd[n][ 4] = (map[j - 1][i] > -1) ? map[j - 1][i] : map[j][i];
				// nbd[n][ 5] = (map[k + 1][j][i] > -1) ? map[k + 1][j][i] : map[k][j][i];
				// nbd[n][ 6] = (map[k - 1][j][i] > -1) ? map[k - 1][j][i] : map[k][j][i];

				// now for more complex directions.
				// minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				nbd[n][ 7 - 2] = (map[j + 1][i + 1] > -1) ? map[j + 1][i + 1] : map[j][i];
				nbd[n][ 8 - 2] = (map[j - 1][i - 1] > -1) ? map[j - 1][i - 1] : map[j][i];
				nbd[n][ 9 - 2] = (map[j - 1][i + 1] > -1) ? map[j - 1][i + 1] : map[j][i];
				nbd[n][10 - 2] = (map[j + 1][i - 1] > -1) ? map[j + 1][i - 1] : map[j][i];
				// nbd[n][11] = (map[k + 1][j][i + 1] > -1) ? map[k + 1][j][i + 1] : map[k][j][i];
				// nbd[n][12] = (map[k - 1][j][i - 1] > -1) ? map[k - 1][j][i - 1] : map[k][j][i];
				// nbd[n][13] = (map[k - 1][j][i + 1] > -1) ? map[k - 1][j][i + 1] : map[k][j][i];
				// nbd[n][14] = (map[k + 1][j][i - 1] > -1) ? map[k + 1][j][i - 1] : map[k][j][i];
				// nbd[n][15] = (map[k + 1][j + 1][i] > -1) ? map[k + 1][j + 1][i] : map[k][j][i];
				// nbd[n][16] = (map[k - 1][j - 1][i] > -1) ? map[k - 1][j - 1][i] : map[k][j][i];
				// nbd[n][17] = (map[k - 1][j + 1][i] > -1) ? map[k - 1][j + 1][i] : map[k][j][i];
				// nbd[n][18] = (map[k + 1][j - 1][i] > -1) ? map[k + 1][j - 1][i] : map[k][j][i];
				/*if (j == 120 && i == 125) {
					for (int xx = 0; xx < 9; ++xx)
					{
						std::cout << nbd[n][xx]  << " " << map[j][i] << std::endl;
					}
				}*/

				n++;


			}
		}


	assert (n == count);  // error if not equal

	return nbd;
}



int ** generate_neighbours_map_2d(unsigned char ** atrium, int ny, int nx, int * count) {
// generates a numeric map of the tissue for all nodes which have an active
// cell.  The result is an ncell*nneighbours (19), sized array, where all cells
// have the indicies of all their neighbours.

	int ** map;
	int ** nbd;

	// generate the numeric map we'll need.
	map = generate_numeric_map_2d(atrium, ny, nx, count);

	// generate the neighbours from it.
	nbd = generate_neighbours_from_numeric_map_2d(map, ny, nx, *count);

	// free the map.
	deallocate_and_zero_2d_matrix(map,  ny + 2, nx + 2); // +2 because this is the way it is defined

	// return the neighbours
	return nbd;
}




void calculate_laplacian_components_with_diffusion_heterogeneity_2d ( unsigned char ** atrium,
        int ny, int nx,
        float **y, float **x,
        double dx, double dy,
        double **dd, double **d2,
        double  ***dc, double ***df,
        double  ***d_dd, double ***d_d2,
        double ***lap) {

	int i, j, n;
	double factor;

	for (j = 1; j < ny + 1; ++j) {
		for (i = 1; i < nx + 1; ++i) {
			if (atrium[j][i] > 0) {
				lap[j][i] = new double [9];// (double *) malloc(19 * sizeof(double));
				if (!lap[j][i]) {
					fprintf(stderr, " %d, %d failed to allocate!\n", j, i);
				}
				for (n = 0; n < 9; ++n) {
					lap[j][i][n] = 0;
				}




				// indexes
				// self         = 0
				// xplus        = 1
				// xminus       = 2
				// yplus        = 3
				// yminus       = 4
				// xplusyplus   = 5; //7 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xminusyminus = 6;//8 -2  // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xplusyminus  = 7;//9 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xminusyplus  = 8;//10 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo


				// dvdx2



				if ( (atrium[j][i + 1] > 0) && (atrium[j][i - 1] > 0) &&
				        (atrium[j - 1][i + 1] > 0) && (atrium[j - 1][i - 1] > 0) &&
				        (atrium[j + 1][i + 1] > 0) && (atrium[j + 1][i - 1] > 0) &&
				        (atrium[j + 1][i] > 0) && (atrium[j - 1][i] > 0)
				   ) {

					// dvdx2 more robust version
					// double factor = 1.0 / 4.0 * (1 / dx) * (1 / dx) * dc[j][i][0];
					// lap[j][i][0] += -4.0 * factor;
					// lap[j][i][1] +=  2.0 * factor;
					// lap[j][i][2] +=  2.0 * factor;
					// lap[j][i][3] += -2.0 * factor;
					// lap[j][i][4] += -2.0 * factor;
					// lap[j][i][5] +=  1.0 * factor;
					// lap[j][i][6] +=  1.0 * factor;
					// lap[j][i][7] +=  1.0 * factor;
					// lap[j][i][8] +=  1.0 * factor;

					// // dvdy2
					// factor = 1.0 / 4.0 * (1 / dy) * (1 / dy) * dc[j][i][3];
					// lap[j][i][0] += -4.0 * factor;
					// lap[j][i][1] += -2.0 * factor;
					// lap[j][i][2] += -2.0 * factor;
					// lap[j][i][3] +=  2.0 * factor;
					// lap[j][i][4] +=  2.0 * factor;
					// lap[j][i][5] +=  1.0 * factor;
					// lap[j][i][6] +=  1.0 * factor;
					// lap[j][i][7] +=  1.0 * factor;
					// lap[j][i][8] +=  1.0 * factor;

										// dvdx2
					factor = (1 / dx) * (1 / dx) * dc[j][i][0];
					// if ((atrium[j][i + 1] > 0) && (atrium[j][i - 1] > 0)) {
						lap[j][i][0] += -2.0 * factor;
						lap[j][i][1] +=  1.0 * factor;
						lap[j][i][2] +=  1.0 * factor;
					// dvdy2
					factor = (1 / dy) * (1 / dy) * dc[j][i][3];
						lap[j][i][0] += -2.0 * factor;
						lap[j][i][3] +=  1.0 * factor;
						lap[j][i][4] +=  1.0 * factor;

				} else {

					factor = (1 / dx) * (1 / dx) * dc[j][i][0];
					if ((atrium[j][i + 1] > 0) && (atrium[j][i - 1] > 0)) {

						lap[j][i][0] += -2.0 * factor;
						lap[j][i][1] +=  1.0 * factor;
						lap[j][i][2] +=  1.0 * factor;

					} else if ((atrium[j][i + 1] > 0) && (atrium[j][i - 1] == 0)) {
						lap[j][i][0] += -1.0 * factor;
						lap[j][i][1] +=  1.0 * factor;
						lap[j][i][2] +=  0.0 * factor;
					} else if ((atrium[j][i + 1] == 0) && (atrium[j][i - 1] > 0)) {
						lap[j][i][0] += -1.0 * factor;
						lap[j][i][1] +=  0.0 * factor;
						lap[j][i][2] +=  1.0 * factor;
					} else {
						lap[j][i][0] += 0;
						lap[j][i][1] += 0;
						lap[j][i][2] += 0;
					}

					// dvdy2
					factor = (1 / dy) * (1 / dy) * dc[j][i][3];
					if ((atrium[j + 1][i] > 0) && (atrium[j - 1][i] > 0)) {
						lap[j][i][0] += -2.0 * factor;
						lap[j][i][3] +=  1.0 * factor;
						lap[j][i][4] +=  1.0 * factor;




					} else if ((atrium[j + 1][i] > 0) && (atrium[j - 1][i] == 0)) {
						lap[j][i][0] += -1.0 * factor;
						lap[j][i][3] +=  1.0 * factor;
						lap[j][i][4] +=  0.0 * factor;
					} else if ((atrium[j + 1][i] == 0) && (atrium[j - 1][i] > 0)) {
						lap[j][i][0] += -1.0 * factor;
						lap[j][i][3] +=  0.0 * factor;
						lap[j][i][4] +=  1.0 * factor;
					} else {
						lap[j][i][0] += 0.0;
						lap[j][i][3] += 0.0;
						lap[j][i][4] += 0.0;
					}
				}


				/*factor = (1 / dx) * (1 / dx) * dc[j][i][0];
				lap[j][i][0] += -2.0 * factor;
				lap[j][i][1] +=  1.0 * factor;
				lap[j][i][2] +=  1.0 * factor;

				factor = (1 / dy) * (1 / dy) * dc[j][i][3];
				lap[j][i][0] += -2.0 * factor;
				lap[j][i][3] +=  1.0 * factor;
				lap[j][i][4] +=  1.0 * factor;*/



				// // dvdz2  - removed // 15:23:35, Tue, 31-December-2019, By Haibo

				// xplusyplus   = 5; //7 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xminusyminus = 6;//8 -2  // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xplusyminus  = 7;//9 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo
				// xminusyplus  = 8;//10 -2   // minus 2 because this is modified from the 3d code. Haibo. // 13:30:49, Tue, 31-December-2019, By Haibo

				/* comment for second order derivatives */
				// dudxdy
				if ( (atrium[j + 1][i + 1] > 0) &&
				        (atrium[j - 1][i - 1] > 0) &&
				        (atrium[j - 1][i + 1] > 0) &&
				        (atrium[j + 1][i - 1] > 0)) {
					factor = (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  +=  factor;
					lap[j][i][8 - 2]  +=  factor;
					lap[j][i][9 - 2]  += -factor;
					lap[j][i][10 - 2] += -factor;
				}
				// dudxdy
				else if ( (atrium[j + 1][i + 1] == 0) &&
				          (atrium[j - 1][i - 1] > 0) &&
				          (atrium[j - 1][i + 1] > 0) &&
				          (atrium[j + 1][i - 1] > 0)) {
					factor = (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  +=  0.0;
					lap[j][i][8 - 2]  +=  factor;
					lap[j][i][9 - 2]  += -0.5 * factor;
					lap[j][i][10 - 2] += -0.5 * factor;
				}
				// dudxdy
				else if ( (atrium[j + 1][i + 1] > 0) &&
				          (atrium[j - 1][i - 1] == 0) &&
				          (atrium[j - 1][i + 1] > 0) &&
				          (atrium[j + 1][i - 1] > 0)) {
					factor = (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  +=  factor;
					lap[j][i][8 - 2]  +=  0.0;
					lap[j][i][9 - 2]  += -0.5 * factor;
					lap[j][i][10 - 2] += -0.5 * factor;
				}
				// dudxdy
				else if ( (atrium[j + 1][i + 1] > 0) &&
				          (atrium[j - 1][i - 1] > 0) &&
				          (atrium[j - 1][i + 1] == 0) &&
				          (atrium[j + 1][i - 1] > 0)) {
					factor = (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  +=  0.5 * factor;
					lap[j][i][8 - 2]  +=  0.5 * factor;
					lap[j][i][9 - 2]  +=  0.0;
					lap[j][i][10 - 2] += -factor;
				}
				// dudxdy
				else if ( (atrium[j + 1][i + 1] > 0) &&
				          (atrium[j - 1][i - 1] > 0) &&
				          (atrium[j - 1][i + 1] > 0) &&
				          (atrium[j + 1][i - 1] == 0)) {
					factor = (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  += 0.5 * factor;
					lap[j][i][8 - 2]  += 0.5 * factor;
					lap[j][i][9 - 2]  += -factor;
					lap[j][i][10 - 2] += 0.0;
				}
				/* else if (
				    (atrium[j + 1][i + 1] == 0) &&
				    (atrium[j + 1][i - 1] == 0) &&
				    (atrium[j - 1][i - 1] > 0) &&
				    (atrium[j - 1][i + 1] > 0)
				)
				{
					factor =  (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][1]  +=  factor;
					lap[j][i][8 - 2]  +=  factor;
					lap[j][i][9 - 2]  += -factor;
					lap[j][i][2] += -factor;
				}
				else if (
				    (atrium[j + 1][i + 1] == 0) &&
				    (atrium[j + 1][i - 1] > 0) &&
				    (atrium[j - 1][i - 1] > 0) &&
				    (atrium[j - 1][i + 1] == 0)
				)
				{
					factor =  (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][0]  +=  factor;
					lap[j][i][6]  +=  factor;
					lap[j][i][0]  += -factor;
					lap[j][i][8] += -factor;
				}

				else if (
				    (atrium[j + 1][i + 1] > 0) &&
				    (atrium[j + 1][i - 1] > 0) &&
				    (atrium[j - 1][i - 1] == 0) &&
				    (atrium[j - 1][i + 1] == 0)
				)
				{
					factor =  (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][5]  +=  factor;
					lap[j][i][0]  +=  factor;
					lap[j][i][0]  += -factor;
					lap[j][i][8] += -factor;
				}
				else if (
				    (atrium[j + 1][i + 1] > 0) &&
				    (atrium[j + 1][i - 1] == 0) &&
				    (atrium[j - 1][i - 1] == 0) &&
				    (atrium[j - 1][i + 1] > 0)
				)
				{
					factor =  (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][5]  +=  factor;
					lap[j][i][4]  +=  factor;
					lap[j][i][7]  += -factor;
					lap[j][i][3] += -factor;
				}
				*/
				else {

					/*std::cout << j << " " << i << " " << int(atrium[j + 1][i + 1]) << " "  << int(atrium[j + 1][i - 1])
					          << " "  << int(atrium[j - 1][i + 1])
					          << " "  << int(atrium[j - 1][i - 1])
					          << " EEE" << std::endl;*/
					factor =  (1.0 / 4.0) * (1.0 / dx) * (1.0 / dy) * (dc[j][i][1] + dc[j][i][2]);
					lap[j][i][7 - 2]  +=  factor;
					lap[j][i][8 - 2]  +=  factor;
					lap[j][i][9 - 2]  += -factor;
					lap[j][i][10 - 2] += -factor;
				}
				// dudxdz  -- removed // 15:25:10, Tue, 31-December-2019, By Haibo

				// dudydz -- removed // 15:25:10, Tue, 31-December-2019, By Haibo


				// right, that's the easy bit done.  Now for the more
				// complex bits!
				// dvdx
				factor = (1.0 / 2.0) * (1 / dx) * (dd[j][i] * (x[j][i] * df[j][i][0] + x[j][i] * df[j][i][0]) +
				                                   x[j][i] * x[j][i] * d_dd[j][i][0] + d_d2[j][i][0] +
				                                   dd[j][i] * (y[j][i] * df[j][i][2] + x[j][i] * df[j][i][3]) +  // [4] -> [3]
				                                   x[j][i] * y[j][i] * d_dd[j][i][1] /*+
				                                   dd[j][i] * (z[j][i] * df[j][i][2] + x[j][i] * df[j][i][8]) +
				                                   x[j][i] * z[j][i] * d_dd[j][i][2]*/  // remove z dimension here // 15:39:28, Tue, 31-December-2019, By Haibo
				                                  );
				/*factor = (1 / 2) * (1 / dx) * dd[j][i] * ((x[j][i] * df[j][i][0] + x[j][i] * df[j][i][0]) +
				         (y[j][i] * df[j][i][1] + x[j][i] * df[j][i][4]) +
				         (z[j][i] * df[j][i][2] + x[j][i] * df[j][i][8]));*/

				// implement upwind scheme here! haibo 2014.11.06
				if ((atrium[j][i + 1] > 0) && (atrium[j][i - 1] > 0)) {

					/*centred difference*/
					/*lap[j][i][0] += 0;
					lap[j][i][1] +=  factor;
					lap[j][i][2] += -factor;*/

					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += -2 * factor;
						lap[j][i][1] += 2 * factor;
						lap[j][i][2] += 0.0;
					} else {
						lap[j][i][0] += 2 * factor;
						lap[j][i][1] += 0.0;
						lap[j][i][2] += -2 * factor;
					}


				} else if ((atrium[j][i + 1] > 0) && (atrium[j][i - 1] == 0)) {
					/*lap[j][i][0] += -2 * factor;
					lap[j][i][1] +=  2 * factor;
					lap[j][i][2] += 0;*/

					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += -2 * factor;
						lap[j][i][1] += 2 * factor;
						lap[j][i][2] += 0.0;
					} else {
						lap[j][i][0] += 0.0;
						lap[j][i][1] += 0.0;
						lap[j][i][2] += 0.0;
					}
				} else if ((atrium[j][i + 1] == 0) && (atrium[j][i - 1] > 0)) {
					/*lap[j][i][0] += 2 * factor;
					lap[j][i][1] += 0;
					lap[j][i][2] += -2 * factor;*/

					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += 0.0;
						lap[j][i][1] += 0.0;
						lap[j][i][2] += 0.0;
					} else {
						lap[j][i][0] += 2 * factor;
						lap[j][i][1] += 0.0;
						lap[j][i][2] += -2 * factor;
					}

				} else {
					lap[j][i][0] += 0;
					lap[j][i][1] += 0;
					lap[j][i][2] += 0;
				}


				// dvdy
				factor = (1.0 / 2.0) * (1 / dy) * (dd[j][i] * (x[j][i] * df[j][i][1] + y[j][i] * df[j][i][0]) +
				                                   x[j][i] * y[j][i] * d_dd[j][i][0] +
				                                   dd[j][i] * (y[j][i] * df[j][i][3] + y[j][i] * df[j][i][3]) +
				                                   y[j][i] * y[j][i] * d_dd[j][i][1] + d_d2[j][i][1] /*+
				                                   dd[j][i] * (z[j][i] * df[j][i][5] + y[j][i] * df[j][i][8]) +
				                                   z[j][i] * y[j][i] * d_dd[j][i][2]*/
				                                  );
				/*factor = (1 / 2) * (1 / dx) * dd[j][i] * ((x[j][i] * df[j][i][3] + y[j][i] * df[j][i][0]) +
				         (y[j][i] * df[j][i][4] + y[j][i] * df[j][i][4]) +
				         (z[j][i] * df[j][i][5] + y[j][i] * df[j][i][8]));*/

				if ((atrium[j + 1][i] > 0) && (atrium[j - 1][i] > 0)) {

					/*centred difference*/
					/*lap[j][i][0] += 0;
					lap[j][i][3] +=  factor;
					lap[j][i][4] += -factor;*/

					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += -2 * factor;
						lap[j][i][3] += 2 * factor;
						lap[j][i][4] += 0.0;
					} else {
						lap[j][i][0] += 2 * factor;
						lap[j][i][3] += 0.0;
						lap[j][i][4] += -2 * factor;
					}

				} else if ((atrium[j + 1][i] > 0) && (atrium[j - 1][i] == 0)) {


					/*lap[j][i][0] += -2 * factor;
					lap[j][i][3] +=  2 * factor;
					lap[j][i][4] += 0;*/
					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += -2 * factor;
						lap[j][i][3] += 2 * factor;
						lap[j][i][4] += 0.0;
					} else {
						lap[j][i][0] += 0.0;
						lap[j][i][3] += 0.0;
						lap[j][i][4] += 0.0;
					}


				} else if ((atrium[j + 1][i] == 0) && (atrium[j - 1][i] > 0)) {
					/*lap[j][i][0] += 2 * factor;
					lap[j][i][3] += 0;
					lap[j][i][4] += -2 * factor;*/
					/* upwind scheme */
					if (factor >= 0.0)
					{

						lap[j][i][0] += 0.0;
						lap[j][i][3] += 0.0;
						lap[j][i][4] += 0.0;
					} else {
						lap[j][i][0] += 2 * factor;
						lap[j][i][3] += 0.0;
						lap[j][i][4] += -2 * factor;
					}

				} else {
					lap[j][i][0] += 0;
					lap[j][i][3] += 0;
					lap[j][i][4] += 0;
				}

				// dvdz  -- removed // 15:27:42, Tue, 31-December-2019, By Haibo

			}
		}
	}
}



unsigned char * generate_tissue_map_2d(unsigned char ** atrium,
                                       int ny, int nx, int count) {
// makes a linear map of the tissue.
	unsigned char * tissue;
	int i, j, c;

	// allocate space for the tissue.
	// tissue = (unsigned char*) malloc(count*sizeof(unsigned char));
	tissue = new unsigned char [count];
	if (!tissue) {
		perror("malloc: tissue map arrray");
		exit(EXIT_FAILURE);
	}

	c = 0;
	for (j = 1; j < ny + 1; j++)
		for (i = 1; i < nx + 1; i++)
			if (atrium[j][i] > 0) {
				tissue[c] = atrium[j][i];
				c++;
			}
	if (c != count) {
		perror("in generating tissue map, the total num of cells does not equal to the counter of input");
		exit(EXIT_FAILURE);
	}
	return tissue;
}



unsigned char * create_stimulation_map_vec_2D_geo(unsigned char ** atrium, const char *filename, int ny, int nx, int count) {

	unsigned char ** stim = read_and_embed_geometry_2d(filename, ny, nx);

	// makes a linear map of the tissue.
	unsigned char * stim_vec;
	int i, j, c;

	// allocate space for the stim_vec.
	// stim_vec = (unsigned char*) malloc(count*sizeof(unsigned char));
	stim_vec = new unsigned char [count];
	if (!stim_vec) {
		perror("new: stim_vec map arrray failed!\n");
		exit(EXIT_FAILURE);
	}

	c = 0;
	for (j = 1; j < ny + 1; j++)
		for (i = 1; i < nx + 1; i++)
			if (atrium[j][i] > 0) {
				stim_vec[c] = stim[j][i];
				c++;
			}
	if (c != count) {
		perror("in generating stim_vec map, the total num of cells does not equal to the counter of input! \n");
		exit(EXIT_FAILURE);
	}


	deallocate_and_zero_2d_matrix(stim, ny + 2, nx + 2);
	return stim_vec;

}


void update_2D_laplacian(double * lap, int num_neighboor, unsigned char ** atrium,
        int ny, int nx, int count,
        const char *phi_name,
        double dx, double dy, double dd, double* d2) {

// generate_laplacian_heterogeneity_2d_test
    double** laplacian = generate_laplacian_heterogeneity_2d_test(atrium,
        ny, nx, count,
        phi_name,
        dx, dy, dd, d2);

    for (int n = 0; n < count; n++) {
        for (int i = 0; i < num_neighboor; i++) {
            lap[(n * num_neighboor) + i] = laplacian[n][i];
        }
    }
    deallocate_and_zero_2d_matrix(laplacian, count, num_neighboor);
}