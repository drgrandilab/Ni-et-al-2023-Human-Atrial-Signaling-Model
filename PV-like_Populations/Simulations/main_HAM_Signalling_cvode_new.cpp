#include <vector>
#include <zlib.h>
// #include <stdio.h>
// #include <stdlib.h>
#include <math.h>
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
// #define IC_FILE




int main(int argc, char const *argv[])
{

	APInfor AP("APD_measure.dat", false, true);
	double stim;
	double v;
	double dv;
	double BCL = 1000.0;



	BCL              = atof(argv[1]);
	int AF           = 0;
	int Popul_ID = atoi(argv[2]);
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

	std::string output_file_name = "AP.BCL." + std::to_string(BCL) + ".ID." + std::to_string(Popul_ID);

	std::ofstream output_file( output_file_name.c_str(), std::ios::out);                     // output filename
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
	Total_time = 600e3 + 60e3 + 1000 ; //601000 + 3000; // (S1_Num - 1) * BCL + 2 * BCL;

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
	// Atrial_wraper.cell.cell_para.IKr_Scale = 0.;
	// Atrial_wraper.cell.cell_para.IK2p_Scale           = 0.0 ;//

	if (argc == 1 + 5 + 25) {
		Atrial_wraper.cell.cell_para.INa_Scale            =  atof(argv[2 + 4]) ; //     *  1; // ;// 0.991102;// 1;
		Atrial_wraper.cell.cell_para.INaL_Scale           =  atof(argv[2 + 5]) ; //     *  1; // 1.;//0.948192;// 1.2;
		Atrial_wraper.cell.cell_para.INab_Scale           =  atof(argv[2 + 6]) ; //     *  1; // 1.;//0.950898;// 1;//1.077453;
		Atrial_wraper.cell.cell_para.INaK_Scale           =  atof(argv[2 + 7]) ; //     *  1; // 1;//0.932910;// 1;//0.9*1.075329;
		Atrial_wraper.cell.cell_para.Itof_Scale           =  atof(argv[2 + 8]) ; //     *  1; // 1;//1.007522;// 1.041409;
		Atrial_wraper.cell.cell_para.IKur_Scale           =  atof(argv[2 + 9]) ; //     *  1; // 1.07;//1.074665;// 1.028833;
		Atrial_wraper.cell.cell_para.IK2p_Scale           =  atof(argv[2 + 10]); //      * 1; // 1.1;//0.984595;//  1.315084;
		Atrial_wraper.cell.cell_para.IKr_Scale            =  atof(argv[2 + 11]); //      * 1; // 0.95 ; //0.934642;//  0.899033;
		Atrial_wraper.cell.cell_para.IKs_Scale            =  atof(argv[2 + 12]); //      * 1; // 0.95 ; //0.946997;//  1.177648;
		Atrial_wraper.cell.cell_para.IK1_Scale            =  atof(argv[2 + 13]); //      * 1; // 0.95 ; //0.960132;//  0.902799;
		Atrial_wraper.cell.cell_para.IKp_Scale            =  atof(argv[2 + 14]); //      * 1; // 1.05;// 1.045922;//  1.095730;
		Atrial_wraper.cell.cell_para.IKach_Scale          =  atof(argv[2 + 15]); //      * 1; // 0.95 ; //0.891492;//  0.895409;
		Atrial_wraper.cell.cell_para.ISK_Scale            =  atof(argv[2 + 16]); //      * 1; // 1.05;//1.082985;//  1.099031;
		Atrial_wraper.cell.cell_para.ICaL_Scale           =  atof(argv[2 + 17]); //      * 1; // 1.037031;//  1.037185;
		Atrial_wraper.cell.cell_para.ICab_Scale           =  atof(argv[2 + 18]); //      * 1; // 1.073250;//  0.997408;
		Atrial_wraper.cell.cell_para.ICap_Scale           =  atof(argv[2 + 19]); //      * 1; // 1.022257;//  1.020868;
		Atrial_wraper.cell.cell_para.INCX_Scale           =  atof(argv[2 + 20]); //      * 1; // 1.05;//1.15;//1.101249;//  0.972194;
		Atrial_wraper.cell.cell_para.IClCa_Scale          =  atof(argv[2 + 21]); //      * 1; // 1.178030;//  0.912827;
		Atrial_wraper.cell.cell_para.IClb_Scale           =  atof(argv[2 + 22]); //      * 1; // 1.1;//0.956994;//  1.1;
		Atrial_wraper.cell.cell_para.Jrel_Scale           =  atof(argv[2 + 23]); //      * 1; // 1.273543;//  0.980100;
		Atrial_wraper.cell.cell_para.Jserca_Scale         =  atof(argv[2 + 24]); //      * 1; // 1.036176;//  1.013384;
		Atrial_wraper.cell.cell_para.Jleak_Scale          =  atof(argv[2 + 25]); //      * 1;// 1.134138;//  0.899613;
		Atrial_wraper.cell.cell_para.Cleft_Buffer_Scale   =  atof(argv[2 + 26]); //      *  1;//1.111454;//  1.301314;
		Atrial_wraper.cell.cell_para.Cytosol_Buffer_Scale =  atof(argv[2 + 27]); //      *  1; // 1.069467;//  0.940860;
		Atrial_wraper.cell.cell_para.SR_Buffer_Scale      =  atof(argv[2 + 28]); //      *  1; // 1.064127;//  0.958842;
	}

	// Atrial_wraper.cell.cell_para.Jserca_Scale = 1.2;

	double stim_old = 0;

	while (outputTime < Total_time) {
		// if (outputTime >= 0.0 && outputTime <= 2.0) stim = -18; else stim = 0.0;
		// if (outputTime > BCL) outputTime = 0.0;

		stim = S1S2(0.0, -12.5, BCL, 10000, S2, outputTime, 5.0);

		if (outputTime > Total_time - 60e3 - 5) {
			Atrial_wraper.cell.ECC_Module.allow_stimulus_app = false;
		}



		// if (outputTime > Total_time - 30e3 - 5)
		// 	time_step = 0.01;
		/*if (outputTime > 100e3) {
			Atrial_wraper.cell.cell_para.cell_betaAR_para.ISO = ISO;
		}
		*/
		/*if (t > 332e3-1) {
			stim = 0;  // stop pacing at the end of the pacing train
		}*/

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

		if (Atrial_wraper.cell.ECC_Module.allow_stimulus_app)
			AP.MeasureAPD90_INa(outputTime, Atrial_wraper.cell.ECC_Module.I_app, BCL, time_step, V, Atrial_wraper.cell.ECC_Module.I_Na, Cai);

		// only output once
		if (outputTime > Total_time - 60e3 - 0.5 * time_step and outputTime < Total_time - 60e3 + 0.5 * time_step) {
			if (Mode == "Normal")
			{

				std::string output_file_name_IC = "Restart_ICs/ICs.bin." +  std::to_string(Popul_ID);
				std::ofstream out(output_file_name_IC.c_str(), std::ios::binary | std::ios::out);
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

		if (count % outfreq == 0 and outputTime >=  598e3-10) // 0*370e3 - 10)
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



	// std::ofstream output_AMP_ratio( "ratio.log");                      // output filename

	// output_AMP_ratio << AP.AMP_over_last << std::endl;
	// output_AMP_ratio.close();
	// std::cout << S2 << " ";
	std::cout << Popul_ID << " ";
	AP.ReportLastThree();

	// std::cout << std::setprecision(10) << Atrial_wraper.cell.y[38] << std::endl;
	// std::cout << "done";
	return 0;
}
