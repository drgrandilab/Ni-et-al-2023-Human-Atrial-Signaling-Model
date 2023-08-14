// Header file for the ICs and RHS functions for the Colman 2013 model

#ifndef ARGUMENTS_H
#define ARGUMENTS_H



#ifdef __INTEL_COMPILER
#define compiler "Intel"
#include <mathimf.h>
#else
#define compiler "GCC"
#include <cmath>
#endif


#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <vector>
#include <algorithm>
class FileSys {
public:
    std::string Pacemap;
    std::string FB_map;
    std::string Fcell_map;
    std::string Geometry_file;
    std::string Fibre_theta;
    std::string Fibre_phi;
    std::string Stim_file;
    std::string SAN_type;
    std::string Stim_time_file;
    std::string Stim_amp_file;
    std::string Apicobasal_file;
    std::string RVIndex_file;
    std::string IScheamiaIndex_file;
    std::string S2_StimLocFile;
    std::string OneD_OutFile;
    FileSys();
    ~FileSys();
};

class Simulation_Config: public FileSys {
public:
    double BCL;
    double Total_time;
    double t_start;
    float dt;
    double ISO_con;
    bool CaMKII_db;
    bool CaMKII_inhb;
    float Ach;
    float Diff_Scale;
    float Drug_Scaling;
    double S2;
    int S1_number;

    int Model_type;
    int model_out;
    int IKur_type;
    int region;
    int region_3D;
    int AF_model;
    int mutation;
    int FB_type;
    int FB_number;
    double Ggap;
    int tau_type;
    bool INa_Drug_Add, IKur_Drug_Add;
    std::string mutation_char;
    std::string region_char;
    std::string tau_type_char;
    std::string ICs;
    std::string Stim_type;
    std::string Model_type_char;
    std::string OneD_OutFile_AppendMode;
    std::string PLB_Sim;
    std::string Output_Folder;


    bool No_CaMKII_K1; //= "Default";
    bool No_CaMKII_NaV; //= "Default";
    bool No_CaMKII_GapG; //= "Default";

    std::map<std::string, std::string> arg;
    std::vector<double> Popul_scalingF;
    int Popl_SF_Number;

    std::vector<double> Remodelling_F;
    int Remodelling_F_Number;

    std::vector<double> Modulation_F;
    int Modulation_F_Number;

    Simulation_Config();
    Simulation_Config(int argc, char *argv[]);  // another constructor
    ~Simulation_Config();
    void Initilise();
    void Config_handling(int argc, char *argv[]);
    void Report_Config();
    void Report_All();

    int Sim_ID;  // ID for identifying/tagging the current simulation 
};




#endif
