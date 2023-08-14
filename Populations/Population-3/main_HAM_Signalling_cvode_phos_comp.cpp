#include <vector>
#include <zlib.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <string>

#include "SingleCellParameter.hpp"
// #include "CRN_ODE.h"
// #include "Maleckar_ODE.h"
// #include "Grandi_2011_ODE.h"
#include  <iomanip>
#include "EnumSimulationCtrl.hpp"
#include "APInfo.hpp"
#include "stimulus.h"
// #include "lsoda_C.hpp"
// #include "Lsoda_function_wrap.h"


// #include "cvode_function.h"

#include "CVOde_Cell.hpp"
#include "cvode_function.h"

// #define OUT_PHASE

#define S1_Num 100
// #define STATENUM 66
#define OUT_AP
#define IC_FILE


int main(int argc, char const *argv[])
{

	APInfor AP("APD_measure.dat", false, true);
	double stim;
	double v;
	double dv;
	double BCL = 1000.0;



	BCL              = atof(argv[1]);
	int AF           = atoi(argv[2]);
	double ISO       = atof(argv[3]);

	int CaMKII_inhit    = atoi(argv[4]);
	int CaMKII_double = atof(argv[5]);

	// double IKurScale = atof(argv[4]);
	// double CaMKII    = atoi(argv[5]);
	// double ISK_Scale = atof(argv[6]);
	int celltype     = 14;
	double S2         = BCL;
	std::string Mut   = "WT";
	std::string Mode  = "Normal";
	int count         = 0;
	double Total_time = 3000;
	int outfreq       = 10;
	double t          = 0.0;
	double outputTime = 0.0;
	std::ofstream output_file( "HAM_wrap_out.dat");                      // output filename
	std::ofstream output_file_LTCC( "LTCC_current.dat");                      // output filename

	// TypeCell celltype    = LVEPI;
	// celltype = (TypeCell)in_celltype;                            // Type of cell  0 - 5
	double time_step     = 0.1;                                 // Time step

	// HAM_Signalling Cell(AF, ISO, false);
	// Cell.ODE_NUM=40;
	int STATENUM =  184;//Cell.ODE_NUM;
	// Cell.cell_para.IKur_Scale = IKurScale;
	// Cell.cell_para.CaMKII_Inhibition = CaMKII;
	// Cell.cell_para.ISK_Scale = ISK_Scale;



	// Cell.SetAFTypePara(NONE);
	// Cell.SetAFTypePara(AF4);

	// sprintf(filename_out, "results/BCL_%d_Model_%s_Region_%d_FB_%d_%d_Mut_%s_AF_%d_AP.txt", Argin.BCL, Argin.Model_type_char, Argin.region_3D, Argin.FB_type, Argin.FB_number, Argin.mutation_char, Argin.AF_model);
	// std::cout << BCL << " " <<celltype << " " << Mut  << ICFilename << std::endl;

	// Total_time = 50 * BCL;
	// double S2 = 0.5 * BCL;
	Total_time = 620000 + 1000 ; //601000 + 3000; // (S1_Num - 1) * BCL + 2 * BCL;

	double Time_Start = 0.0;

	/*std::ifstream infile("Restart_ICs/ICs.bin", std::ios::binary | std::ios::in);
	if (infile.is_open())
	{
		infile.read(reinterpret_cast<char*>(Cell.y), sizeof(double)*STATENUM);
		infile.close();
		// Cell.V = state[38];
		// Time_Start = Total_time - 3 * BCL - time_step;
	} else {
		std::cout << "ICs file not opened..." << std::endl;
	}*/


	// LSODA lsoda(STATENUM);

	// std::cout << "Time_Start = " << Time_Start << std::endl;
	double tout = 0;

	// for (t = Time_Start; t < Total_time; t = t + time_step)
	t = Time_Start;
	// t = 0;


	std::ofstream CaM_dyad_out( "CaM_dyad_out.dat");                      // output filename
	std::ofstream CaM_sl_out( "CaM_sl_out.dat");                      // output filename
	std::ofstream CaM_cyto_out( "CaM_cyto_out.dat");                      // output filename
	std::ofstream CaMKII_out( "CaMKII_out.dat");                      // output filename
	std::ofstream betaAR_out( "betaAR_out.dat");                      // output filename
	std::ofstream para_out( "para_out.dat");                      // output filename



	// cvode_solver cvode(STATENUM, 1);

	// cvode.set_IC(Cell.y);

	// for (int i = 0; i < STATENUM; i++)
	// std::cout << Ith(cvode.y, i + 1) << std::endl;

	// cvode.initialise_mem(f);
	// cvode.set_user_data(&Cell);


	CVOde_Cell Atrial_wraper(STATENUM, 1, f, true);
	Atrial_wraper.cell.BCL = BCL;
	Atrial_wraper.cell.assign_AF_ISO_CaMKII_para(AF, ISO, CaMKII_double);
	Atrial_wraper.cell.cell_para.CaMKII_Inhibition = CaMKII_inhit;
	// Atrial_wraper.cell.cell_para.No_CaMKII_NaV = true;
	// Atrial_wraper.cell.cell_para.No_CaMKII_K1 = true;


	Atrial_wraper.read_initial_condition("Restart_ICs/ICs_baseline.bin");
	// Atrial_wraper.cell.cell_para.Itof_Scale           = 1.0;//  1.;//atof(argv[2 + 8])
	// Atrial_wraper.cell.cell_para.IKur_Scale           = 0.5;//   0.5;//atof(argv[2 + 9])
	// Atrial_wraper.cell.cell_para.ISK_Scale           =  0.5;//atof(argv[2 + 9])
	// Atrial_wraper.cell.cell_para.ICaL_Scale = 0.9;
	// // Atrial_wraper.cell.cell_para.IKr_Scale = 0.5;


	double scale[22 + 3];

	if (argc == 22 + 3 + 3 + 3) {
		for (int i = 0; i < 22 + 3; ++i)
		{
			scale[i] = atof(argv[i + 3 + 3]);

			// std::cerr << scale[i] << std::endl;
		}
		Atrial_wraper.cell.betaAR_Module.set_scale_parameters(scale);

		Atrial_wraper.cell.betaAR_Module.print_scale_parameters();

		Atrial_wraper.cell.cell_para.cell_CaM_para_dyad.set_scale_parameters(&scale[22]);
		Atrial_wraper.cell.cell_para.cell_CaM_para_sl.set_scale_parameters(&scale[22]);
		Atrial_wraper.cell.cell_para.cell_CaM_para_cyto.set_scale_parameters(&scale[22]);


		std::cerr << "**********************************" << std::endl;

		Atrial_wraper.cell.cell_para.cell_CaM_para_dyad.print_scale_parameters();


	}


	while (outputTime < Total_time) {
		// if (outputTime >= 0.0 && outputTime <= 2.0) stim = -18; else stim = 0.0;
		// if (outputTime > BCL) outputTime = 0.0;

		stim = S1S2(0.0, -12.5, BCL, 10000, S2, outputTime, 5.0);

		if (outputTime > Total_time - 20e3 - 5) {
			Atrial_wraper.cell.ECC_Module.allow_stimulus_app = false;
		}


		/*if (outputTime > Total_time - 30e3 - 5)
			time_step = 0.1;*/
		/*if (outputTime > 100e3) {
			Atrial_wraper.cell.cell_para.cell_betaAR_para.ISO = ISO;
		}
		*/
		/*if (t > 332e3-1) {
			stim = 0;  // stop pacing at the end of the pacing train
		}*/



		// Atrial_wraper.cell.CaM_Module_dyad.print_content();
		// Atrial_wraper.cell.CaM_Module_sl.print_content();
		/*if (t>10 and t<11)
		{
			lsoda.SetIstate( 3);
		}*/
		// Cell.ECC_Module.I_app = -stim;

		outputTime = outputTime + time_step;
		tout = tout + time_step;
		// lsoda.lsoda(lsoda_generic_HAM_Signalling_ODE, Cell.y - 1, &t, tout, &Cell); // lsoda arrays start from 1
		// cvode.solve_single_step(outputTime);
		Atrial_wraper.solve_single_time_step(tout, time_step);
		// CVode(cvode.cvode_mem, outputTime, cvode.y, &t, CV_NORMAL); // 1 time step solution.
		// f(t, cvode.y, cvode.ydot, &Cell);
		// t = outputTime;
		// std::cout << outputTime << std::endl;
		// Cell.V += Cell.dV * time_step;

		Atrial_wraper.cell.ECC_Module.V  = Atrial_wraper.cell.ECC_Module.y[38];


		// has to get values from cvode directly!!!!
		double V = Ith(Atrial_wraper.cvode.y, 38 + 1);
		double Cai = Ith(Atrial_wraper.cvode.y, 37 + 1);
		// std::cout << Cell.ECC_Module.y[38] << std::endl;

		// AP.MeasureAPD90_INa(outputTime, Atrial_wraper.cell.ECC_Module.I_app, BCL, time_step, Atrial_wraper.cell.ECC_Module.y[38], Atrial_wraper.cell.ECC_Module.I_Na, Atrial_wraper.cell.ECC_Module.y[37]);
		AP.MeasureAPD90_INa(outputTime, Atrial_wraper.cell.ECC_Module.I_app, BCL, time_step, V, Atrial_wraper.cell.ECC_Module.I_Na, Cai);

		// only output once
		if (outputTime > Total_time - 50 * BCL - 0.5 * time_step and outputTime < Total_time - 50 * BCL + 0.5 * time_step) {
			if (Mode == "Normal")
			{
				std::ofstream out("Restart_ICs/ICs.bin", std::ios::binary | std::ios::out);
				if (out.is_open())
				{
					// state[38] =  Cell.V;
					out.write(reinterpret_cast<char*>(Atrial_wraper.cell.y), sizeof(double)*STATENUM);
					out.close();
				} else {
					std::cerr << "Restart_ICs/ICs.bin not successfully opened!!!\n" ;
				}

#ifdef IC_FILE


				char ICFilename[200];
				sprintf(ICFilename, "BCL_%.0f_Model_CRN_Region_%d_Mut_%s_AF_%d_ICs.txt", BCL, celltype, Mut.c_str(), AF);
				std::ofstream output_file2(ICFilename);                      // output filename

				for (int i = 0; i < STATENUM; ++i)
				{
					output_file2 << std::setprecision(10) << Atrial_wraper.cell.y[i] << std::endl;
				}
				output_file2.close();


#endif

			}
		}

#ifdef OUT_AP

		if (fabs(Atrial_wraper.cell.ECC_Module.dV) > 1)
			outfreq = int (1 / time_step);
		else
			outfreq = int(4 / time_step);

		// outfreq = 2;

		if (count % outfreq == 0 and outputTime >=  1 * Total_time - 100e3 ) // 0*370e3 - 10)
		{
			// Atrial_wraper.cell.ECC_Module.print_to_file(outputTime, output_file);
			Atrial_wraper.cell.ECC_Module.print_to_file_Vm_Ca(outputTime, output_file);
			// Atrial_wraper.cell.ECC_Module.LTCC.print_to_file(outputTime, output_file_LTCC);

			// Atrial_wraper.cell.CaM_Module_dyad.print_to_file(outputTime, CaM_dyad_out);
			// Atrial_wraper.cell.CaM_Module_sl.print_to_file(outputTime, CaM_sl_out);
			// Atrial_wraper.cell.CaM_Module_cyto.print_to_file(outputTime, CaM_cyto_out);
			// Atrial_wraper.cell.CaMKII_Module.print_to_file(outputTime, CaMKII_out);
			// Atrial_wraper.cell.betaAR_Module.print_to_file(outputTime, betaAR_out);
			// Atrial_wraper.cell.cell_para.print_to_file(outputTime, para_out);




		}
#endif
		count ++;

	}

	output_file.close();



	std::ofstream output_AMP_ratio( "ratio.log");                      // output filename

	output_AMP_ratio << AP.AMP_over_last << std::endl;
	output_AMP_ratio.close();
	// std::cout << S2 << " ";

	AP.ReportLastThree();

	// std::cout << std::setprecision(10) << Atrial_wraper.cell.y[38] << std::endl;
	// std::cout << "done";
	return 0;
}
