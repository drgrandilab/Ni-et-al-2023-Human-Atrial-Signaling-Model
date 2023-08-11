/*
 * stimulus.h
 * Haibo Ni <qiangzini@gmail.com>
 *
 * generate stimulus according to S1 or S1S2 proctol
 *
 * 
 *	Last Update: Wed 22 Jun 2016 19:12:51 BST
 */

#ifndef STIMULUS_H
#define STIMULUS_H

#include <math.h>
#include <stdio.h>
// #include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
/*  all the times are in ms */
 /* all current in pA/pF */

double S1(double stim_start, double stim, double BCL, double current_time, double stim_duration);

double S1S2(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, double current_time, double stim_duration);
double S1S2_num(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, int S2_num, double current_time, double stim_duration) ;

double S1S2_num_stim(double stim_start_time, double BCL_s1, double stim_current_1, int S1_num, 
						double BCL_s2, double stim_current_2, int S2_num,
						double current_time, double stim_duration) ;
double simple_S1(double stim_start, double stim, double BCL, double current_time, double stim_duration);




/* this class reads in a txt file containning all the stimulation cycle length information and then apply the stimulatin.
useage:

	StimFromInputFile stimulusClass ("stimulus.txt",0);
stim = stimulusClass. ApplyStim(-18, 2.0, 950,t);
 inputs are :: (double stim_current, double stim_duration, double time_start, double current_time) ;
author: Haibo Ni Wed 22 Jun 2016 19:24:25 BST

*/
class StimFromInputFile
{
public:
    StimFromInputFile(const char *, bool report=false);
    ~StimFromInputFile() {};
   double ApplyStim(double stim_current, double stim_duration, double time_start, double current_time) ;
    std::vector<double> CL;
    double BCL;
    int StimNum, TotalNum;
    double StartTime, TimeSinceStim;
    bool stim_ON, stim_ON_prev;
    bool _inprogram_report;
};


/*typedef struct
{
	double ERP;
	double BCL_s1, BCL_s2;
	double ap_amp_ratio;
	double BCL_adjust;
} ERP_measure_struct;
*/

#endif   // END OF STIMULUS.H 

