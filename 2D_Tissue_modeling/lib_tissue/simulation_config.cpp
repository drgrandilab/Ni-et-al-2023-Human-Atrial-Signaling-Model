#include "simulation_config.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

FileSys::FileSys() {
	/*Geometry_file = new char [100];
	Fibre_theta  = new char [100];
	Fibre_phi    = new char [100];
	Pacemap      = new char [100];
	Stim_file    = new char [100];
	SAN_type = new char[100];
	Stim_time_file = new char[100];
	Stim_amp_file = new char[100];
	Fcell_map = new char [100];
	FB_map = new char [100];
	Apicobasal_file = new char [100];*/
}
FileSys::~FileSys() {

	/* delete [] Geometry_file;
	 delete [] Fibre_theta ;
	 delete [] Fibre_phi   ;
	 delete [] Pacemap     ;
	 delete [] Stim_file   ;
	 delete [] SAN_type;
	 delete [] Stim_time_file;
	 delete [] Stim_amp_file;
	 delete [] Fcell_map;
	 delete [] FB_map;
	 delete [] Apicobasal_file;*/
}


Simulation_Config::Simulation_Config() {
	Initilise();
}

void Simulation_Config::Initilise() {

	BCL        = 1000;
	S1_number  = 10;
	S2         = BCL;
	dt         = 0.005;

	Total_time = 2000.0;
	t_start = 0.0;
	Model_type = 2;
	IKur_type  = 0;
	model_out  = 0;
	region     = 3;
	region_3D  = 202;
	AF_model   = 0;
	mutation   = 0;
	FB_type    = 0;
	FB_number  = 0;
	Ggap       = 3.0;
	tau_type   = 1;
	Diff_Scale = 1.0;
	Drug_Scaling = 1.0;


	ISO_con = 0.0;
	CaMKII_db = false;
	CaMKII_inhb = false;
	region_char     = "RAA";
	Model_type_char = "Colman";
	tau_type_char   = "Fast";
	// ICs             = "function_defined";
	Pacemap         = "wholeHeart_091148011-10.pace.gz";
	Stim_type       = "Paced";
	SAN_type        = "Heterogeneous";
	mutation_char   = "None";
	ICs             = "Default"; // Vent_seeman13.theta.gz
	Stim_time_file  = "Geometry/Vent_seeman.time_2.bin";
	Stim_amp_file   = "Geometry/Vent_seeman.volt_2.bin";
	Geometry_file   =  "Geometry/Vent_seeman_Add_RV_Stim_RG.mat.gz";
	Fibre_theta     = "Geometry/Vent_seeman13.theta.gz";
	Fibre_phi       = "Geometry/Vent_seeman13.phi.gz";
	Apicobasal_file = "Geometry/Vent_seeman_Apicobasal_file.bin";
	RVIndex_file    = "Geometry/RV_Gradient_Index.bin";
	IScheamiaIndex_file = "Geometry/Ischeamia_Conditions_None.bin";
	S2_StimLocFile = "Geometry/S2_loc.mat.gz";
	arg["IscheamiaPhase"] = "Phase_1A";
	IKur_Drug_Add = false;
	INa_Drug_Add = false;

	Popl_SF_Number = 0;  // no Popl_scaling factors provided.
	Remodelling_F_Number = 0;
	Modulation_F_Number = 0;

	OneD_OutFile = "OneD_APD_out.dat";
	OneD_OutFile_AppendMode = "std::fstream::app";

	Sim_ID = 0;  // ID of the current simulation;

	PLB_Sim = "Default";
	Output_Folder = "./";

	No_CaMKII_K1   = false;  //"Default";
	No_CaMKII_NaV  = false;  //"Default";
	No_CaMKII_GapG = false;  //"Default";


}


Simulation_Config::Simulation_Config(int argc, char *argv[]) {
	Initilise();
	// Simulation_Config();
	Config_handling(argc, argv);

}
Simulation_Config::~Simulation_Config() {
	/*delete [] Model_type_char;
	delete [] region_char    ;
	delete [] tau_type_char  ;
	delete [] mutation_char  ;
	delete [] ICs            ;
	delete [] Stim_type      ;*/
}


void Simulation_Config::Config_handling(int argc, char *argv[]) {

	int arg_counter;

	arg_counter = 1;
	while (arg_counter < argc) {
		std::string command = argv[arg_counter];
		if (command == "BCL") {
			BCL = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Total_time") {
			Total_time = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Model_type") {
			Model_type_char = argv[arg_counter + 1];
			if (Model_type_char == "Colman") {
				IKur_type = 0;
				model_out = 0;
			}
			else if (Model_type_char == "CNZ") {
				IKur_type = 1;
				model_out = 1;
			}
			else if (Model_type_char == "Colman_v") Model_type = 4;
			else {
				std::cerr << "Invalid Model type\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		} // end Model_type if
		else if (command == "Region") {
			region_char = argv[arg_counter + 1];
			if (region_char == "PM") {
				region = 1;
				region_3D = 12;
			}
			else if (region_char == "CT") {
				region = 2;
				region_3D = 11;
			}
			else if (region_char == "RAA") {
				region = 3;
				region_3D = 202;
			}
			else if (region_char == "AVR") {
				region = 4;
				region_3D = 18;
			}
			else if (region_char == "BB") {
				region = 5;
				region_3D = 15;
			}
			else if (region_char == "LA") {
				region = 6;
				region_3D = 16;
			}
			else if (region_char == "AS") {
				region = 7;
				region_3D = 17;
			}
			else if (region_char == "LAA") {
				region = 8;
				region_3D = 201;
			}
			else if (region_char == "PV") {
				region = 11;
				region_3D = 101;
			}
			else if (region_char == "PV_jones") region = 9;
			else if (region_char == "SAN_C") {
				region = 14;
				region_3D = 10;
			}
			else if (region_char == "SAN_P") {
				region = 15;
				region_3D = 1011;
			}
			arg_counter++;
		}
		else if (command == "AF") {
			AF_model = atoi(argv[arg_counter + 1]);
			if (AF_model > 4) {
				std::cerr << "AF model must be 0, 1, 2, 3 or 4\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "Mutation") {
			mutation_char = argv[arg_counter + 1];

			if (mutation_char == "D322H") mutation = 1;
			else if (mutation_char == "E48G") mutation = 2;
			else if (mutation_char == "A305T") mutation = 3;
			else if (mutation_char == "Y155C") mutation = 4;
			else if (mutation_char == "D469E") mutation = 5;
			else if (mutation_char == "P488S") mutation = 6;
			else if (mutation_char == "None") mutation  = 0;
			else if (mutation_char == "A545P") mutation  = 10;
			else {
				std::cerr << "Invalid mutation\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "FB_type") {
			FB_type = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "FB_number") {
			FB_number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Ggap") {
			Ggap = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Tau_type") {
			tau_type_char = argv[arg_counter + 1];
			if (tau_type_char ==  "Slow") tau_type = 0;
			else if (tau_type_char ==  "Fast") tau_type = 1;
			else {
				std::cerr << "invalid tau type";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}
		else if (command == "BCL") {
			BCL = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "S1") {
			BCL = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "S2") {
			S2 = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "SAN_type") {
			SAN_type = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_type") {
			Stim_type = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "S1_number") {
			S1_number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Diffusion_Scale") {
			Diff_Scale = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "Pacemap") {
			Pacemap = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_Time_File") {
			Stim_time_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Stim_Amp_File") {
			Stim_amp_file = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "Apicobasal_File") {
			Apicobasal_file = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "RVIndex_File") {
			RVIndex_file = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "IScheamiaIndex_File") {
			IScheamiaIndex_file = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "Geometry_File") {
			Geometry_file = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "S2_StimLocFile") {
			S2_StimLocFile = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "Time_Start") {
			t_start = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "IscheamiaPhase") {
			arg["IscheamiaPhase"] = argv[arg_counter + 1];
			arg_counter++;
		} else if (command == "INa_Drug") {
			// arg["Add_INa_Drug"] = argv[arg_counter + 1];
			std::string drug = argv[arg_counter + 1];
			if (drug == "True")
			{
				INa_Drug_Add = true;
				std::cout << "true" << std::endl;
			} else if (drug == "False")
			{
				INa_Drug_Add = false;
			} else {
				std::cerr << drug << " is not a valid INa drug argument. Valid arguments are as follows:\n\t "
				          << "True or False" << std::endl;
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}

		else if (command == "IKur_Drug") {
			std::string drug = argv[arg_counter + 1];
			if (drug == "True")
			{
				IKur_Drug_Add = true;
				std::cout << "true" << std::endl;
			} else if (drug == "False")
			{
				IKur_Drug_Add = false;
			} else {
				std::cerr << drug << " is not a valid IKur drug argument. Valid arguments are as follows:\n\t "
				          << "True or False" << std::endl;
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		}	else if (command == "Drug_Scaling") {
			Drug_Scaling = atof(argv[arg_counter + 1]);
			arg_counter++;
		}
		else if (command == "ICs")  {
			ICs = argv[arg_counter + 1];
			arg_counter++;
		}
		else if (command == "OneD_OutFile")  {
			OneD_OutFile = argv[arg_counter + 1];
			arg_counter++;
		}
		/*else if (command == "OneD_OutFile_AppendMode")  {
			OneD_OutFile_AppendMode = argv[arg_counter + 1];
			arg_counter++;
		}*/
		else if (command == "Sim_ID")  {
			Sim_ID = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}

		else if (command == "Popl_SF_Number")  {
			Popl_SF_Number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}  else if (command == "Popul_scaling_factors")  {
			// Popl_SF_Number = atoi(argv[arg_counter + 1]);

			if (Popl_SF_Number <= 0)
			{
				std::cerr << "Popl_SF_Number = " << Popl_SF_Number << std::endl;
				std::cerr << "Please assign Popl_SF_Number with 'Popl_SF_Number' command first! " << std::endl;
				std::cerr << "Programe exiting ... " << std::endl;
				std::exit(0);
			} else {

				for (int i = 0; i < Popl_SF_Number; ++i)
				{
					Popul_scalingF.push_back(atof(argv[arg_counter + 1]));
					arg_counter++;
				}
			}
		}

		else if (command == "Remodelling_F_Number")  {
			Remodelling_F_Number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}  else if (command == "Remodelling_Factors")  {
			// Popl_SF_Number = atoi(argv[arg_counter + 1]);

			if (Remodelling_F_Number <= 0)
			{
				std::cerr << "Remodelling_F_Number = " << Remodelling_F_Number << std::endl;
				std::cerr << "Please assign Remodelling_F_Number with 'Remodelling_F_Number' command first! " << std::endl;
				std::cerr << "Programe exiting ... " << std::endl;
				std::exit(0);
			} else {

				for (int i = 0; i < Remodelling_F_Number; ++i)
				{
					Remodelling_F.push_back(atof(argv[arg_counter + 1]));
					arg_counter++;
				}
			}
		}

		else if (command == "Modulation_F_Number")  {
			Modulation_F_Number = atoi(argv[arg_counter + 1]);
			arg_counter++;
		}  else if (command == "Modulation_Factors")  {
			// Popl_SF_Number = atoi(argv[arg_counter + 1]);

			if (Modulation_F_Number <= 0)
			{
				std::cerr << "Modulation_F_Number = " << Modulation_F_Number << std::endl;
				std::cerr << "Please assign Modulation_F_Number with 'Modulation_F_Number' command first! " << std::endl;
				std::cerr << "Programe exiting ... " << std::endl;
				std::exit(0);
			} else {

				for (int i = 0; i < Modulation_F_Number; ++i)
				{
					Modulation_F.push_back(atof(argv[arg_counter + 1]));
					arg_counter++;
				}
			}
		} else if (command == "PLB_Sim") {
			PLB_Sim = argv[arg_counter + 1];
			arg_counter++;

			std::vector<std::string> allowed_str {"kmf.1.0.ISO.0.0", "kmf.0.25.ISO.0.0", "kmf.1.0.ISO.1.0", "kmf.0.25.ISO.1.0"};

			// if not found in the list of allowed strings.
			if ( (std::find(allowed_str.begin(), allowed_str.end(), PLB_Sim) == allowed_str.end())) {
				std::cerr << "Input PLB_Sim option error!! >> " <<  PLB_Sim << std::endl;
				std::cerr << "Allowed are: \n" << "kmf.1.0.ISO.0.0" << "\n" << "kmf.0.25.ISO.0.0"
				          << "\n" << "kmf.1.0.ISO.1.0" << "\n" << "kmf.0.25.ISO.1.0" << std::endl;
				std::exit(0);
			}
		} else if (command == "ISO_con") {
			ISO_con = atof(argv[arg_counter + 1]);
			arg_counter++;
		} else if (command == "CaMKII") {
			// ISO_con = atof(argv[arg_counter + 1]);
			std::string temp_str = argv[arg_counter + 1];
			if ( temp_str == "CaMKII_db")
				CaMKII_db = true;
			else if (temp_str == "CaMKII_inhb")
				CaMKII_inhb = true;
			else if (temp_str == "Default") {
				CaMKII_inhb = false;
				CaMKII_db = false;
			}
			else {
				std::cerr << temp_str << " is not a valid CaMKII argument. Valid arguments are as follows:\n" <<
				          "\tCaMKII_db\n" << "\tCaMKII_inhb\n" << "\tDefault\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		} else if (command == "CaMKII_K1_NaV_GapG") {

			// ISO_con = atof(argv[arg_counter + 1]);
			std::string temp_str = argv[arg_counter + 1];
			if ( temp_str == "No_CaMKII_K1")
				No_CaMKII_K1 = true;
			else if (temp_str == "No_CaMKII_NaV")
				No_CaMKII_NaV = true;
			else if (temp_str == "No_CaMKII_GapG") {
				No_CaMKII_GapG = true;
			} else if (temp_str == "Default") {
				No_CaMKII_K1 = false;
				No_CaMKII_NaV = false;
				No_CaMKII_GapG = false;
			}
			else {
				std::cerr << temp_str << " is not a valid CaMKII CaMKII_K1_NaV_GapG. Valid arguments are as follows:\n" <<
				          "\t No_CaMKII_K1\n" << "\tNo_CaMKII_NaV\n" << "\tNo_CaMKII_GapG\n" << "\tDefault\n";
				std::exit(EXIT_FAILURE);
			}
			arg_counter++;
		} else if (command == "Output_Folder") {

			Output_Folder = argv[arg_counter + 1];
			Output_Folder = Output_Folder + "/";

			arg_counter++;


		} else {
			std::cerr << argv[arg_counter] << " is not a valid argument. Valid arguments are as follows:\n\t "
			          << "BCL\n\tTotal_time\n\tModel_type\n\tRegion\n\tAF\n\tMutation\n\tFB_type\n\tFB_number\n\tGgap\n\tTau_type\n"
			          << "please check simulation_config.cpp for more information\n";
			std::exit(EXIT_FAILURE);
		}
		arg_counter ++;
	} // end while

} // end Argument_handling

void Simulation_Config::Report_Config() {



	std::cerr << " S1\t\t = " << BCL << std::endl;
	std::cerr << " S2\t\t = " << S2 << std::endl;
	std::cerr << " Total_time\t = " << Total_time << std::endl;
	std::cerr << " Time_Start\t = " << t_start << std::endl;
	// std::cerr << "dt\t\t = " << dt << std::endl;
	std::cerr << " ISO_con\t = " << ISO_con << std::endl;
	std::cerr << " CaMKII_db\t = " << CaMKII_db << std::endl;
	std::cerr << " CaMKII_inhb\t = " << CaMKII_inhb << std::endl;

	std::cerr << " No_CaMKII_K1\t = " <<  No_CaMKII_K1 << std::endl;
	std::cerr << " No_CaMKII_NaV\t = " << No_CaMKII_NaV << std::endl;
	std::cerr << " No_CaMKII_GapG\t = " << No_CaMKII_GapG << std::endl;
	std::cerr << std::endl;
	std::cerr << std::endl;

	// std::cerr << "Region\t\t = " << region << std::endl;
	// std::cerr << "Region_char\t = " << region_char << std::endl;
	// std::cerr << "AF_model\t = " << AF_model << std::endl;
	// std::cerr << "tau_type_char\t = " << tau_type_char << std::endl;
	// std::cerr << "Diff_Scale\t = " << Diff_Scale << std::endl;
	// if (FB_number != 0)  {
	// 	std::cerr << "Fibroblasts on:\n\t"
	// 	          << "FB_type = " << FB_type
	// 	          << "\n\tFB_number = " << FB_number
	// 	          << "\n\tGgap = " << Ggap << std::endl;
	// }

	// std::cerr << "Pacemap\t\t = " << Pacemap << std::endl
	//           << "Geometry_file\t = " << Geometry_file << std::endl
	//           << "Apicobasal_file\t = " << Apicobasal_file << std::endl
	//           << "Fibre_theta\t = "  << Fibre_theta << std::endl
	//           << "Fibre_phi\t = "  << Fibre_phi << std::endl
	//           <<  "Stim_type\t= " << Stim_type << std::endl
	//           <<  "Stim_Time_File\t= " << Stim_time_file << std::endl
	//           <<  "Stim_Amp_File\t= " << Stim_amp_file << std::endl
	//           <<  "IKur_Drug\t= " << IKur_Drug_Add << std::endl
	//           <<  "INa_Drug\t= " << INa_Drug_Add << std::endl;

	/*	std::cerr << "BCL\t\t = " << BCL << std::endl;
		std::cerr << "Total_time\t = " << Total_time << std::endl;
		std::cerr << "Time_Start\t = " << t_start << std::endl;
		std::cerr << "dt\t\t = " << dt << std::endl;
		std::cerr << "Model_type\t = " << Model_type_char << std::endl;
		std::cerr << "Region\t\t = " << region << std::endl;
		std::cerr << "Region_char\t = " << region_char << std::endl;
		std::cerr << "AF_model\t = " << AF_model << std::endl;
		std::cerr << "tau_type_char\t = " << tau_type_char << std::endl;
		std::cerr << "Diff_Scale\t = " << Diff_Scale << std::endl;
		if (FB_number != 0)  {
			std::cerr << "Fibroblasts on:\n\t"
			          << "FB_type = " << FB_type
			          << "\n\tFB_number = " << FB_number
			          << "\n\tGgap = " << Ggap << std::endl;
		}

		std::cerr << "Pacemap\t\t = " << Pacemap << std::endl
		          << "Geometry_file\t = " << Geometry_file << std::endl
		          << "Apicobasal_file\t = " << Apicobasal_file << std::endl
		          << "Fibre_theta\t = "  << Fibre_theta << std::endl
		          << "Fibre_phi\t = "  << Fibre_phi << std::endl
		          <<  "Stim_type\t= " << Stim_type << std::endl
		          <<  "Stim_Time_File\t= " << Stim_time_file << std::endl
		          <<  "Stim_Amp_File\t= " << Stim_amp_file << std::endl
		          <<  "IKur_Drug\t= " << IKur_Drug_Add << std::endl
		          <<  "INa_Drug\t= " << INa_Drug_Add << std::endl;*/
	// <<  "Stim_Amp_File\t= " << Stim_amp_file << std::endl;
}


void Simulation_Config::Report_All() {
	std::cerr << "BCL\t\t = " << BCL << std::endl;
	std::cerr << "Total_time\t = " << Total_time << std::endl;
	std::cerr << "dt\t\t = " << dt << std::endl;
	std::cerr << "Model_type\t = " << Model_type_char << std::endl;
	std::cerr << "Region\t\t = " << region << std::endl;
	std::cerr << "Region_char\t = " << region_char << std::endl;
	std::cerr << "AF_model\t = " << AF_model << std::endl;
	std::cerr << "tau_type_char\t = " << tau_type_char << std::endl;
	std::cerr << "Diff_Scale\t = " << Diff_Scale << std::endl;
}
