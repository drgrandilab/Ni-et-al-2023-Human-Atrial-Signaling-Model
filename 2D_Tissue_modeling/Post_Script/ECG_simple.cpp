/*
* A very simple code to calculate pseudo ecgs.
* The conductivity of tissue was not considered.

* Haibo Ni <haibo.ni0822@gmail.com>

Thu 02 Jul 2015 10:33:56 BST

* compile: g++ PsedoEcg.cpp -lz -std=c++11 -o ECG
* run:    ./ECG xloc yloc zloc potential_file_folder filestart fileend
 * icpc ECG_simple.cpp -std=c++11 -lz


 * this code does not consider the border of the tissue, so if you want to add the contribution of the border tissue to ecg, please add
 * empty tissue around the border, just a single layer should be fine.
*/



// #include "PsudoEcg.h"
#include <vector>
#include <zlib.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <tuple>
#include <iomanip>
#include <sstream>

/*#define NX (462)
#define NY (325)
#define NZ (358)*/

/*
#define NX (235)
#define NY (269)
#define NZ (298)
*/


// #define NX (301)
// #define NY (301)
// #define NZ (1)

#define NX 125//203//;(462)
#define NY 120//203//;(325)
#define NZ 1//3//;(358)

typedef std::tuple<int, int, int> cell_map;


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

double compute_ecg(const std::vector<cell_map> & map, float * potential, float x_loc, float y_loc, float z_loc) ;

/*#define NX1 (100)
#define NY1 (100)
#define NZ1 (100)*/


// #define D2 (0.21)
// #define DD (8*0.21)
/* to be determined later */
/* The diffusion parameter */
#define D1 (0.18)
#define D2 (0.20*0.33/1.5)
#define DD (D1-D2)
#define dx 0.25  // change the spatial step accordingly. haibo
float potential_mat [NZ][NY][NX];
unsigned char geo[NZ][NY][NX];
int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		std::cerr <<  "please run the code with ->  ./ECG xloc yloc zloc potential_file_folder filestart fileend" << std::endl;
		std::exit(1);
	}


	float x_loc, y_loc, z_loc;
	x_loc = atof(argv[1]), y_loc = atof(argv[2]), z_loc = atof(argv[3]);

	gzFile gz = NULL;
	// std::string filename = "ATRIUM_SAN.geo.gz";// "Geometry/Vent_seeman_Add_RV_Stim_RG.mat.gz";
	// std::string filename = "./Geometry/TWO_D_Plat_251X251X1.geo.gz";//TWO_D_Plat_301X301X1.geo.gz";// "Geometry/Vent_seeman_Add_RV_Stim_RG.mat.gz";
	std::string filename = "../Geometry/TWO_D_Plat_125X120.geo.gz";
	gz = gzopen(filename.c_str(), "r");

	if (!gz) {
		printf("Geometry file not opened!!!\n");
		std::exit(0);
	}
	gzread(gz, geo, NZ * NY * NX);
	gzclose(gz);
	std::vector <cell_map> map;
	int count = 0;
	for ( int g = 0; g < NZ; g++) {
		for ( int h = 0; h < NY; h++) {
			for ( int i = 0 ; i < NX; i++) {
				if (geo[g][h][i] != 0 ) {
					map.push_back(std::make_tuple(i, h, g));
					count ++;
				}
			}
		}
	}

	std::cout << "total Num cell " << count << std::endl;
	std::string foldername = argv[4];

	std::string output_name =foldername+ to_string_with_precision(x_loc) + "_" + to_string_with_precision(y_loc) + "_" + to_string_with_precision(z_loc) + "_" + ".dat";
	// output_name.replace(output_name.find("/"), std::string("/").size(), "_");
	// std::cout << output_name << std::
	std::ofstream ecg_out (output_name.c_str(), std::ofstream::out);
	int start = atoi(argv[5]);
	int end = atoi(argv[6]);
	for (int i = start; i <= end; i++ ) {

		std::cout << i << std::endl;
		char int_string [10];
		snprintf(int_string, sizeof(int_string), "%04d", i);
		std::cout << int_string << std::endl;

		std::string potential_filename = foldername + "/v_" + std::string(int_string) + ".bin";
		std::cout << potential_filename << std::endl;

		float * potential = NULL;
		if (potential == NULL)
			potential = new float[count];
		std::ifstream data_in;
		data_in.open(potential_filename.c_str(), std::ios::binary);
		if (data_in.is_open())
			data_in.read( reinterpret_cast<char*>( potential ), sizeof(float) * count );
		else {
			std::cerr << potential_filename << " doesnot exist!!! " << std::endl;
			continue;
		}
		// std::cout << potential_filename << std::endl;

		data_in.clear();
		data_in.close();
		std::cout << 3 << std::endl;

		double ECG = compute_ecg(map, potential, x_loc,  y_loc,  z_loc);
		std::cout << 2 << std::endl;

		std::cout << potential_filename << "  Produced Ecg_singal ->    " << ECG << std::endl;
		ecg_out << i <<  "\t" << ECG << std::endl;

		if (potential != NULL)
		{
			/* code */
			delete [] potential;
		}
	}
	ecg_out.close();
	return 0;
}


double compute_ecg(const std::vector<cell_map> &map, float * potential, float x_loc, float y_loc, float z_loc) {
	int count = 0;
	int i, h, g;
	double ECG = 0.0;

	for (int j = 0; j <= map.size(); j++) {
		i = std::get<0>(map[j]);
		h = std::get<1>(map[j]);
		g = std::get<2>(map[j]);
		potential_mat[g][h][i] = potential[j];
	}
	for (int j = 0; j <= map.size(); j++) {
		i = std::get<0>(map[j]);
		h = std::get<1>(map[j]);
		g = std::get<2>(map[j]);
		// potential_mat[g][h][i] = potential[j];
		if (NZ > 1) {
			if (i > 0 and i < NX - 1 and h > 0 and h < NY - 1 and g > 0 and g < NZ - 1)
			{	double gradVx = 0.0;
				double gradVy = 0.0;
				double gradVz = 0.0;
				// gradVx
				if (geo[g][h][i + 1] != 0 and geo[g][h][i - 1] != 0)
					gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i - 1]) / (2 * dx);
				else if (geo[g][h][i + 1] == 0 and geo[g][h][i - 1] == 0)
					gradVx = 0.0;
				else if (geo[g][h][i + 1] == 0)
					gradVx = (potential_mat[g][h][i] - potential_mat[g][h][i - 1]) / (dx);
				else
					gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i]) / (dx);
				// gradVy
				if (geo[g][h + 1][i] != 0 and geo[g][h - 1][i] != 0)
					gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h - 1][i]) / (2 * dx);
				else if (geo[g][h + 1][i] == 0 and geo[g][h - 1][i] == 0)
					gradVy = 0.0;
				else if (geo[g][h + 1][i] == 0)
					gradVy = (potential_mat[g][h][i] - potential_mat[g][h - 1][i]) / (dx);
				else
					gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h][i]) / (dx);



				// gradVz
				if (geo[g + 1][h][i] != 0 and geo[g - 1][h][i] != 0)
					gradVz = (potential_mat[g + 1][h][i] - potential_mat[g - 1][h][i]) / (2 * dx);
				else if (geo[g + 1][h][i] == 0 and geo[g - 1][h][i] == 0)
					gradVz = 0.0;
				else if (geo[g + 1][h][i] == 0)
					gradVz = (potential_mat[g][h][i] - potential_mat[g - 1][h][i]) / (dx);
				else
					gradVz = (potential_mat[g + 1][h][i] - potential_mat[g][h][i]) / (dx);

				double x = i * dx;
				double y = h * dx;
				double z = g * dx;

				double drx = 0;
				double dry = 0;
				double drz = 0;
				double r = std::sqrt((x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc) + (z - z_loc) * (z - z_loc) );
				if (r == 0.0)
					r = 1e-10;

				drx = -(x - x_loc) / (r * r * r);
				dry = -(y - y_loc) / (r * r * r);
				drz = -(z - z_loc) / (r * r * r);

				ECG -= (gradVx * drx + gradVy * dry + gradVz * drz) * dx * dx * dx;
			}
		} else {
			if (i > 0 and i < NX - 1 and h > 0 and h < NY - 1)
			{	double gradVx = 0.0;
				double gradVy = 0.0;
				double gradVz = 0.0;
				// gradVx
				if (geo[g][h][i + 1] != 0 and geo[g][h][i - 1] != 0)
					gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i - 1]) / (2 * dx);
				else if (geo[g][h][i + 1] == 0 and geo[g][h][i - 1] == 0)
					gradVx = 0.0;
				else if (geo[g][h][i + 1] == 0)
					gradVx = (potential_mat[g][h][i] - potential_mat[g][h][i - 1]) / (dx);
				else
					gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i]) / (dx);
				// gradVy
				if (geo[g][h + 1][i] != 0 and geo[g][h - 1][i] != 0)
					gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h - 1][i]) / (2 * dx);
				else if (geo[g][h + 1][i] == 0 and geo[g][h - 1][i] == 0)
					gradVy = 0.0;
				else if (geo[g][h + 1][i] == 0)
					gradVy = (potential_mat[g][h][i] - potential_mat[g][h - 1][i]) / (dx);
				else
					gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h][i]) / (dx);

				// gradVz
				/*if (geo[g + 1][h][i] != 0 and geo[g - 1][h][i] != 0)
					gradVz = (potential_mat[g + 1][h][i] - potential_mat[g - 1][h][i]) / (2 * dx);
				else if (geo[g + 1][h][i] == 0 and geo[g - 1][h][i] == 0)
					gradVz = 0.0;
				else if (geo[g + 1][h][i] == 0)
					gradVz = (potential_mat[g][h][i] - potential_mat[g - 1][h][i]) / (dx);
				else
					gradVz = (potential_mat[g + 1][h][i] - potential_mat[g][h][i]) / (dx);*/

				gradVz = 0.0; // assign to zero.

				double x = i * dx;
				double y = h * dx;
				double z = g * dx;

				double drx = 0;
				double dry = 0;
				double drz = 0;
				double r = std::sqrt((x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc) + (z - z_loc) * (z - z_loc) );
				if (r == 0.0)
					r = 1e-10;

				drx = -(x - x_loc) / (r * r * r);
				dry = -(y - y_loc) / (r * r * r);
				drz = -(z - z_loc) / (r * r * r);

				ECG -= (gradVx * drx + gradVy * dry + gradVz * drz) * dx * dx * dx;
			}
		}

	}

	/*for ( int g = 0; g < NZ; g++) {
		for ( int h = 0; h < NY; h++) {
			for ( int i = 0 ; i < NX; i++) {
				if (geo[g][h][i] != 0 ) {
					potential_mat[g][h][i] = potential[count];
					count ++;
				}
			}
		}
	}*/



	/*for ( int g = 1; g < NZ - 1; g++) {
		for ( int h = 1; h < NY - 1; h++) {
			for ( int i = 1 ; i < NX - 1; i++) {
				if (geo[g][h][i] != 0 ) {
					double gradVx = 0.0;
					double gradVy = 0.0;
					double gradVz = 0.0;
					// gradVx
					if (geo[g][h][i + 1] != 0 and geo[g][h][i - 1] != 0)
						gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i - 1]) / (2 * dx);
					else if (geo[g][h][i + 1] == 0 and geo[g][h][i - 1] == 0)
						gradVx = 0.0;
					else if (geo[g][h][i + 1] == 0)
						gradVx = (potential_mat[g][h][i] - potential_mat[g][h][i - 1]) / (dx);
					else
						gradVx = (potential_mat[g][h][i + 1] - potential_mat[g][h][i]) / (dx);
					// gradVy
					if (geo[g][h + 1][i] != 0 and geo[g][h - 1][i] != 0)
						gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h - 1][i]) / (2 * dx);
					else if (geo[g][h + 1][i] == 0 and geo[g][h - 1][i] == 0)
						gradVy = 0.0;
					else if (geo[g][h + 1][i] == 0)
						gradVy = (potential_mat[g][h][i] - potential_mat[g][h - 1][i]) / (dx);
					else
						gradVy = (potential_mat[g][h + 1][i] - potential_mat[g][h][i]) / (dx);



					// gradVz
					if (geo[g + 1][h][i] != 0 and geo[g - 1][h][i] != 0)
						gradVz = (potential_mat[g + 1][h][i] - potential_mat[g - 1][h][i]) / (2 * dx);
					else if (geo[g + 1][h][i] == 0 and geo[g - 1][h][i] == 0)
						gradVz = 0.0;
					else if (geo[g + 1][h][i] == 0)
						gradVz = (potential_mat[g][h][i] - potential_mat[g - 1][h][i]) / (dx);
					else
						gradVz = (potential_mat[g + 1][h][i] - potential_mat[g][h][i]) / (dx);

					double x = i * dx;
					double y = h * dx;
					double z = g * dx;

					double drx = 0;
					double dry = 0;
					double drz = 0;
					double r = std::sqrt((x - x_loc) * (x - x_loc) + (y - y_loc) * (y - y_loc) + (z - z_loc) * (z - z_loc) );
					if (r == 0.0)
						r = 1e-10;

					drx = -(x - x_loc) / (r * r * r);
					dry = -(y - y_loc) / (r * r * r);
					drz = -(z - z_loc) / (r * r * r);

					ECG -= (gradVx * drx + gradVy * dry + gradVz * drz) * dx * dx * dx;
				}
			}
		}
	}*/

	return ECG;
}
