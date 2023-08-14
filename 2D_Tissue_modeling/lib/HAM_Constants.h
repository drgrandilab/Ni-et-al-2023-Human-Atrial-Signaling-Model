#ifndef HAM_CONSTANTS_H
#define HAM_CONSTANTS_H





//// Model Parameters
// Constants
const double R = 8314;       // [J/kmol*K]
const double Frdy = 96485;   // [C/mol]
const double Temp = 310;     // [K]
const double FoRT = Frdy / R / Temp;
const double Qpow = (Temp - 310) / 10;



// Geometry
// geom_flag = 0;
// if geom_flag == 0
// Capacitance
double const  Acell = 11e3; // [um^2]
double const  Cmem = Acell * 1e-14; // [F] 110 pF membrane capacitance
//Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10 in ventricles

// Fractional currents in compartments
double const  Fjunc = 0.11;
double const Fsl = 1 - Fjunc;
double const  Fjunc_CaL = 0.9;
double const Fsl_CaL = 1 - Fjunc_CaL;

// Cell dimensions and volume
double const  cellLength = 100;     // cell length [um]113;//100
double const  cellRadius = 10.25;   // cell radius [um]12;//10.25
double const  Vcell = 3.14159265 * cellRadius * cellRadius * cellLength * 1e-15; // [L]
double const  Vmyo = 0.65 * Vcell;
double const  Vsr = 0.035 * Vcell;
double const  Vsl = 0.02 * Vcell;
double const  Vjunc = 0.0539 * 0.01 * Vcell;

// Diffusion rates
double const diffusion_scale = 1.0;
double const J_ca_juncsl = diffusion_scale * 1 / 1.2134e12; // [L/msec]           8.2413e-13
double const J_ca_slmyo = diffusion_scale * 1 / 2.68510e11; // [L/msec]           3.2743e-12
double const J_na_juncsl = 1 / (1.6382e12 / 3 * 100); // [L/msec]   1.8313e-14
double const J_na_slmyo = 1 / (1.8308e10 / 3 * 100); // [L/msec]   1.6386e-12
// else
//     // Capacitance
//     Acell = 11e3; // [um^2]
//     Cmem = Acell*1e-14; // [F] 110 pF membrane capacitance
//     //Cmem = 1.1e-10;   // [F] membrane capacitance 1.3810e-10 in ventricles
//
//     // Fractional currents in compartments
//     Fjunc = 0.11; Fsl = 1-Fjunc;
//     Fjunc_CaL = 0.9; Fsl_CaL = 1-Fjunc_CaL;
//
//     // Cell dimensions and volume
//     cellLength = 100; // cell length [um]
//     cellRadius = 10.25-2.25; // cell radius [um]
//     Vcell = pi*cellRadius^2*cellLength*1e-15; // [L]
//     Vmyo = 0.65*Vcell;
//     Vsr = 0.035*Vcell;
//     Vsl = 0.02*Vcell;
//     Vjunc = 0.0539*0.01*Vcell;
//
//     junctionLength = 15e-3; // junc length [um]
//     junctionRadius = 160e-3; // junc radius [um]
//
//     distJuncSL = 0.5; // dist. junc to SL [um]
//     distSLcyto = 0.45; // dist. SL to cytosol [um] - TO BE INCREASED?
//
//     DcaJuncSL = 1*1.64e-6; // Dca junc to SL [cm^2/sec]
//     DcaSLcyto = 1.22e-6; // Dca SL to cyto [cm^2/sec]
//     DnaJuncSL = 1.09e-5; // Dna junc to SL [cm^2/sec]
//     DnaSLcyto = 1.79e-5; // Dna SL to cyto [cm^2/sec]
//
//     SAsl = Fsl*Acell; // [um^2]
//     Njunc = (Fjunc*Acell)/(pi*junctionRadius^2); // [-]
//     SAjunc = Njunc*pi*2*junctionLength*junctionRadius; // [um^2]
//
//     // Diffusion rates
//     J_ca_juncsl = DcaJuncSL*SAjunc/distJuncSL*1e-10; // [L/msec] [m^2/sec]
//     J_ca_slmyo = DcaSLcyto*SAsl/distSLcyto*1e-10; // [L/msec]
//     J_na_juncsl = DnaJuncSL*SAjunc/distJuncSL*1e-10; // [L/msec]
//     J_na_slmyo = DnaSLcyto*SAsl/distSLcyto*1e-10; // [L/msec]
// end

// Fixed ion concentrations
double const Cli = 15;   // Intracellular Cl  [mM]
double const Clo = 150;  // Extracellular Cl  [mM]
double const Ko = 5.4;   // Extracellular K   [mM]
double const Nao = 140;  // Extracellular Na  [mM]
double const Cao = 1.8;  // Extracellular Ca  [mM]
//Ko = 4;   // Extracellular K   [mM]
//Nao = 147;  // Extracellular Na  [mM]
//Cao = 2;  // Extracellular Ca  [mM]
double const  Mgi = 1;    // Intracellular Mg  [mM]



#endif

