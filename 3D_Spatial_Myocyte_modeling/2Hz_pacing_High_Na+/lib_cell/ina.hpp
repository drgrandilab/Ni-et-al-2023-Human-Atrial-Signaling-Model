#ifndef INA_HPP
#define INA_HPP

#include <math.h>

class ina
{
public:
    ina();

    double Update_INa(double dt, double v, double NaV_CKp = 0.12);
    ~ina() {};

    double state1_ina, state2_ina, state3_ina;

    double am, bm;
    double ah, bh;
    double aj, bj;

    double gating;
};


class inaL {

public:
    inaL();
    ~inaL() {};

    double state1_ina;
    double state3_ina;
    double gating;

    double update_INaL(double dt, double v);

};
        // double GNaL = para.INaL_Scale * ( (14.0 / 11.0) / (1 + exp(-(para.NaV_CKp - 0.12) / 0.1)) ) * 0.999763767766896*1.02704368132844 * (1.2 * 0.003)* 9; // [mS/uF] (adjusted) MOD1  // increased GNaL here // 18:29:43, Thu, 07-February-2019, By Haibo
#endif