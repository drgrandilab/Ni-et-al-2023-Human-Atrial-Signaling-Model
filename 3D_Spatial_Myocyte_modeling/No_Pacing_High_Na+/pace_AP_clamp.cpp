#include "Spatial_Cell.h"
#include "stimulus.h"
#include <math.h>
#include <array>
#include <bits/stdc++.h>
#include <iostream>
using namespace std;

int main(int argc, char const *argv[])
{
        // omp_set_num_threads(14); // core number

        double dt      = 0.01;
        int count      = 0;
        int out_interval = int(1.0 / dt);
        double voltage_rest     = -70;
        double stim_current     = 12.5;
        double stim_duration    = 5;
        double pre_rest         = 4000;        // 5000;
        double bcl              = atof(argv[1]);         // pacing rate – 1HZ
        int iter                = atoi(argv[2]);
        double delay_rest       = 5000;
        double T_total          = pre_rest + bcl * iter + delay_rest;

        int flag_AF             = 0; //atoi(argv[4]); // 0: nSR; 1: cAF // Fri 15 Apr 2022 09:00:06 AM PDT
        int flag_female         = 0; //atoi(argv[5]); // 0: male; 1: female // Fri 15 Apr 2022 09:00:09 AM PDT

        // Spatial_Cell cell(dt, 0, 0, bcl, argv[3]); //, flag_AF, flag_female);


        Spatial_Cell cell(dt, 0, 0, bcl, argv[3]);

        std::string CaMKII_file = std::string(argv[4]);

        cell.CKII_para.read_CaMKII_levels_from_txt(CaMKII_file);
        cell.CKII_para.set_CaMKII_signaling_effects();
        cell.CKII_para.print_to_file(0, "out.dat");
        cell.cc.para->print_to_file(0, "outcc.dat");
        cell.sc.para->print_to_file(0, "outsc.dat");


        // 0: Ramp-step Voltage-clamp
        // 1: Triangular Voltage-clamp (constant duration)
        // 2: -80mV Voltage-clamp
        // 3: Step voltage-clamp
        // 4: AP-clamp
        // 5: Triangular Voltage-clamp (bcl-dependent duration)
        int protocol_flag = 2;

        // Hz   APD90   AP-length
        // 0.5  296     326 
        // 1    267     294
        // 2    227     250
        // 3    210     231
        // 4    186     205
        // 5    164     180
        vector<double> tri_ap_rate_array        = {2000, 1000, 500, 333.333, 250, 200};
        vector<double> tri_ap_duration_array    = {326, 294, 250, 231, 205, 180};
        double tri_ap_duration_const            = 200;

        // Find the BCL-dependent AP length
        if (protocol_flag == 5)
        {
                auto it = find(tri_ap_rate_array.begin(), tri_ap_rate_array.end(), bcl);
  
                // If element was found
                if (it != tri_ap_rate_array.end()) 
                {
                        // calculating the index
                        int index = it - tri_ap_rate_array.begin();
                        tri_ap_duration_const = tri_ap_duration_array[index];
                } else 
                {
                        // If the element is not present in the vector
                        cout << "Error protocol_flag or BCL input" << endl;
                }
        }

        cout << tri_ap_duration_const << endl;

        // int AP_size = int(bcl/dt);
        // double AP_set [AP_size] = { };
        std::vector<double> AP_set;

        int AP_count = 0;

        if (protocol_flag == 4)
        {
                string AP_file          = "./pool_AP/AP_";
                string AP_file_end      = ".txt";
                AP_file.append(argv[1]);
                AP_file.append(AP_file_end);
                // std::cout << AP_file << std::endl;

                std::ifstream inAP(AP_file.c_str());// "./pool_AP/AP_3hz.txt"
                std::string APline;

                double tmp = 0;

                while (std::getline(inAP, APline))
                {
                        std::stringstream ssAP(APline);
                        ssAP >> tmp;
                        AP_set.push_back(tmp);
                        // std::cout << AP_set[AP_count] << std::endl;
                        ++AP_count;
                }
        }

	for (double t = 0; t < T_total; t += dt)
	{
		cell.dt = dt;
		double volt   = -80.0;
                if (protocol_flag==0)
                {
                        if ( t < pre_rest )
                        {
                                volt    = -80.0;
                        }
                        else if ( t < pre_rest + bcl*iter - 1 && ( int ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt ) ) < int( 100 / dt ) ) && ( int(t/dt) < int( T_total/dt - delay_rest /dt) ) )
                        {
                                volt    = -80.0 + 40.0 * ( ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt )) / ( 100 / dt ) );  
                        }
                        else if ( t < pre_rest + bcl*iter - 1 && ( int ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt ) ) < int( 200 / dt ) ) && ( int(t/dt) < int( T_total/dt - delay_rest /dt) ) )
                        {
                                volt    = 10.0;  
                        }
                        else{
                                volt    = -80.0; 
                        }
                } else if (protocol_flag == 1 || protocol_flag == 5)
                {
                        if ( t < pre_rest )
                        {
                                volt    = -80.0;
                        }
                        else if ( t < pre_rest + bcl*iter - 1 && ( int ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt ) ) < int( tri_ap_duration_const / dt ) ) && ( int(t/dt) < int( T_total/dt - delay_rest /dt) ) )
                        {
                                volt    = 10.0 - 90.0 * ( ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt )) / ( tri_ap_duration_const / dt ) );  
                        }
                        else{
                                volt    = -80.0; 
                        }
                } else if (protocol_flag == 2)
                {
                        volt    = -80.0;
                } else if (protocol_flag == 3)
                {
                        if ( t < pre_rest )
                        {
                                volt    = -80.0;
                        }
                        else if ( t < pre_rest + bcl*iter - 1 && ( int ( (int(t/dt) - int(pre_rest / dt) ) % int( bcl / dt ) ) < int( 200 / dt ) ) && ( int(t/dt) < int( T_total/dt - delay_rest /dt) ) )
                        {
                                volt    = -80.0 + 10 * (1 + floor((t - pre_rest) / bcl));
                                std::cout << volt << std::endl;
                        }
                        else{
                                volt    = -80.0; 
                        }
                } else if (protocol_flag == 4)
                {
                        if ( t < pre_rest )
                        {
                                volt    = AP_set[0];
                        } else if ( int(t/dt) < int( T_total/dt - delay_rest /dt))
                        {
                                int AP_set_id = int((t-pre_rest)/dt) % int(bcl/dt);
                                volt    = AP_set[AP_set_id];
                        } else{
                                volt    = AP_set[AP_count-1];
                        }
                } else
                {
                        std::cout << "Error protocol_flag" << std::endl;
                }
                
                cell.pace_Vm_clamp(volt);

		std::string nnn;
		if (count % out_interval == 0) {
			cell.output_Cai(nnn);
		}

                //if (t > pre_rest + bcl*iter - 2000 && count % out_interval == 0)
                //{
		//	cell.output_binary(nnn);
		//}

                if( (t > T_total - delay_rest - 0.5 * dt) && (t < T_total - delay_rest + 0.5 * dt) )
                {
                        cell.sc.get_steady_state();
                        // cell.cc.get_steady_state();  
                }

		count ++;
	}
	return 0;
}
