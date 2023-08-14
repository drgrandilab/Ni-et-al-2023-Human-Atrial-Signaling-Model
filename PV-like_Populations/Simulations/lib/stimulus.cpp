#include <math.h>
#include <stdio.h>

#include "stimulus.h"




/*  implementation of S1  */
double S1(double stim_start_time, double stim_current, double BCL, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) >= 0 ? (current_time - stim_start_time) : -1.0;
    remain =  fmod(time_elapsed, BCL);
    if (remain >= 0 && remain < stim_duration) {
        Istim = stim_current;
    } else {
        Istim = 0.0;
    }

    return Istim;
}

double S1S2(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        }
    }
    else if ((time_elapsed > BCL_s1 * (S1_num - 1) + stim_duration) && (time_elapsed <= BCL_s1 * (S1_num - 1) + BCL_s2  + stim_duration) )
    {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        };
    }
    else
    {
        Istim = 0.0;
    }
    return Istim;
}


double S1S2_num(double stim_start_time, double stim_current, double BCL_s1, int S1_num, double BCL_s2, int S2_num, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        }
    }
    else if ((time_elapsed > BCL_s1 * (S1_num - 1) + stim_duration) && (time_elapsed <= BCL_s1 * (S1_num - 1) + BCL_s2 * S2_num + stim_duration) )
    {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current;
        } else {
            Istim = 0.0;
        };
    }
    else
    {
        Istim = 0.0;
    }
    return Istim;
}

double S1S2_num_stim(double stim_start_time, double BCL_s1, double stim_current_1, int S1_num,
                     double BCL_s2, double stim_current_2, int S2_num,
                     double current_time, double stim_duration) {

    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) > 0 ? (current_time - stim_start_time) : 0;
    if (time_elapsed <= BCL_s1 * (S1_num - 1) + stim_duration)
    {
        remain =  fmod(time_elapsed, BCL_s1);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current_1;
        } else {
            Istim = 0.0;
        }
    }
    else if ((time_elapsed > BCL_s1 * (S1_num - 1) + stim_duration) && (time_elapsed <= BCL_s1 * (S1_num - 1) + BCL_s2 * S2_num + stim_duration) )
    {
        time_elapsed = time_elapsed - BCL_s1 * (S1_num - 1) ;
        remain =  fmod(time_elapsed, BCL_s2);
        if (remain >= 0 && remain <= stim_duration) {
            Istim = stim_current_2;
        } else {
            Istim = 0.0;
        };
    }
    else
    {
        Istim = 0.0;
    }
    return Istim;
}

double simple_S1(double * stim_lapse, double stim_start, double stim, double BCL, double current_time, double stim_duration, double dt) {
    /*double Istim = 0.0;
    if ((* stim_lapse) >= st && (* stim_lapse) <= 5.0) Istim = 12.5; else Istim = 0.0;
    if ((* stim_lapse) > BCL) (* stim_lapse) = 0.0;
    (* stim_lapse) = (* stim_lapse) + dt;*/
    return 0;
}




/*  implementation of S1  */



StimFromInputFile::StimFromInputFile(const char *file, bool report) {
    stim_ON = stim_ON_prev = false;
    StimNum = 0;
    TotalNum = 1;
    // StartTime = TimeSinceStim = 0.0;
    // StartTime = time_start;
    std::ifstream iFile(file);
    CL.push_back(0);
    _inprogram_report = report;
    if (iFile.is_open()) {
        /*        while (true) {
                    double x;
                    iFile >> x;

                    CL.push_back(x + CL[TotalNum - 1]);
                    TotalNum++;
                    std::cerr << x << std::endl;
                    if ( iFile.eof() ) break;
                }*/
        double x;
        while ( iFile >> x) {
            CL.push_back(x + CL[TotalNum - 1]);
            TotalNum++;

            if (_inprogram_report)

                std::cerr << x << std::endl;
            // if ( iFile.eof() ) break;
        }
        if (CL.size() < 2) {
            std::cerr << " no stim infor provided!!" << std::endl;
            std::exit(0);
        }

        if (_inprogram_report)
            for (int i = 0; i <= CL.size()-1; i++) {
                std::cerr << CL[i] << std::endl;
            }
    } else {
        std::cerr << "could not open file " << file << std::endl;
        std::exit(0);
    }
}

double StimFromInputFile::ApplyStim(double stim_current, double stim_duration, double time_start, double current_time) {

    double Istim;
    stim_ON_prev = stim_ON ;
    if (current_time >= time_start + CL[StimNum]  and current_time <= time_start + stim_duration + CL[StimNum] ) {
        Istim = stim_current;
        stim_ON = true;
    } else {
        Istim = 0.0;
        stim_ON = false;
    }

    if (stim_ON == false and stim_ON_prev == true) {
        // stim has just been applied and completed.

        StimNum ++;
        if (StimNum > TotalNum) {
            std::cerr << "StimNum > TotalNum :: " << StimNum << " " << TotalNum << std::endl;
        }
        if (_inprogram_report)
            std::cerr << "Applying stim " << StimNum << " out of " << TotalNum << std::endl;

    }
    return Istim;
}


/*
double StimFromInputFile(double stim_start_time, double stim_current, double BCL, double current_time, double stim_duration) {
    double Istim = 0.0;
    double remain;
    double time_elapsed;

    time_elapsed = (current_time - stim_start_time) >= 0 ? (current_time - stim_start_time) : -1.0;
    remain =  fmod(time_elapsed, BCL);
    if (remain >= 0 && remain < stim_duration) {
        Istim = stim_current;
    } else {
        Istim = 0.0;
    }

    return Istim;
}
*/