#include "Spatial_Cell.h"
#include "stimulus.h"

int main(int argc, char const *argv[])
{
        double dt      = 0.01;
        int count      = 0;

        double stim_current     = 12.5;
        double stim_duration    = 5;
        double pre_rest         = 15000;
        double bcl              = atof(argv[1]);
        int iter                = atoi(argv[2]);
        double delay_rest       = 10000;
        double T_total          = pre_rest + bcl * iter + delay_rest;

        Spatial_Cell cell(dt, 0, 0, bcl, argv[3]);

        std::string CaMKII_file = std::string(argv[4]);

        int ISO = 0;
        ISO = atoi(argv[5]);

        cell.CKII_para.read_CaMKII_levels_from_txt(CaMKII_file);
        cell.CKII_para.set_CaMKII_signaling_effects();
        if (ISO)
                cell.CKII_para.set_ISO();
        cell.CKII_para.print_to_file(0, "out.dat");
        cell.cc.para->print_to_file(0, "outcc.dat");
        cell.sc.para->print_to_file(0, "outsc.dat");

        // std::exit(0);
        for (double t = 0; t < T_total; t += dt) {
                cell.dt = dt;
                double stim   = 0;
                if ( t >= pre_rest && t < T_total - delay_rest )
                {
                        stim    = S1(pre_rest, stim_current, bcl, t, stim_duration);
                }

                std::string nnn;
                if (count % 100 == 0)
                {
                        cell.output_Cai(nnn);
                }

                if (count % 200 == 0)
                {
                        cell.output_binary(nnn);
                }

                cell.output_ina_dvdt(nnn, bcl);

                cell.pace(stim);

                if ( (t > T_total - delay_rest - 0.5 * dt) && (t < T_total - delay_rest + 0.5 * dt) )
                {
                        cell.sc.get_steady_state();
                        cell.cc.get_steady_state();
                }

                count ++;
        }

        return 0;
}