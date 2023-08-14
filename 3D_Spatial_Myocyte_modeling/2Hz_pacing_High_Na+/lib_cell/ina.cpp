#include "ina.hpp"

ina::ina() {
        state1_ina      = 0.0;
        state2_ina      = 0.0202919327386704;
        state3_ina      = 0.936490564532038;
        gating          = (state1_ina * state1_ina * state1_ina)  * state2_ina * state3_ina;
}

double ina::Update_INa(double dt, double v, double NaV_CKp){

        double Ina_h_shift = (3.25 - 2 * 3.25 * 1 / (1 + exp(-(NaV_CKp - 0.12) / 0.1))); 
        am      = (v == -47.13) * ( 3.2 )
                + (v != -47.13) * ( 0.32 * (v + 47.13) / (1.0 - exp( -0.1 * (v + 47.13) ) ) );
        bm      = 0.08 * exp( -v / 11.0 );

        double vm = v - Ina_h_shift;

        ah      = (vm >= -40.0) * ( 0.0 )
                + (vm < -40.0) * ( 0.135 * exp( -( vm + 80.0 ) / 6.8 ) );
        bh      = (vm >= -40.0) * ( 1.0 / ( 0.13 * ( 1.0 + exp( -(vm + 10.66) / 11.1 ) ) ) )
                + (vm < -40.0) * ((3.56 * exp( 0.079 * vm) + 3.1e5 * exp(0.35 * vm)));

        aj      = (vm >= -40.0) * (0.0) 
                +(vm < -40.0) * (( ( -127140 * exp(0.2444*vm) - 3.474e-5 * exp(-0.04391 * vm)) * (vm + 37.78)) /
                (1.0 + exp( 0.311 * (vm + 79.23) ) ));
        bj      = (vm >= -40.0) * ((0.3 * exp(-2.535e-7*vm)) / (1.0 + exp( -0.1 * (vm + 32.0) )))
                + (vm < -40.0) * ((0.1212 * exp( -0.01052 * vm )) / (1.0 + exp( -0.1378 * (vm + 40.14) )));
        // Rush-Larsen (RL) mothod
        state1_ina = ( am/(am + bm) ) + (state1_ina - ( am/(am + bm) ) ) * exp(-(dt) * (am + bm) );
        state2_ina = ( ah/(ah + bh) ) + (state2_ina - ( ah/(ah + bh) ) ) * exp(-(dt) * (ah + bh) );
        state3_ina = ( aj/(aj + bj) ) + (state3_ina - ( aj/(aj + bj) ) ) * exp(-(dt) * (aj + bj) );

        gating = ( state1_ina * state1_ina * state1_ina )  * state2_ina * state3_ina;

        return gating;

}


inaL::inaL() {
        state1_ina      = 0.0;
        state3_ina      = 0.936490564532038;
        gating          = (state1_ina * state1_ina * state1_ina)  * state3_ina;
}


double inaL::update_INaL(double dt, double v){

        int AF = 0;
        //// Late I_Na
        //GNaL = SA_par(2)*0.0025*AF; // [mS/uF]
        // (current)  //(1.0 + 0.5 / (1 + exp((para.kCKII_Nav_change - 0.3) / -0.05)))

        double am = 0.32 * (v + 47.13) / (1 - exp(-0.1 * (v + 47.13)));
        double bm = 0.08 * exp(-v / 11);
        double inf = 1 / (1 + exp((v + 91) / 6.1));
        double tau = 200;
        state1_ina = ( am/(am + bm) ) + (state1_ina - ( am/(am + bm) ) ) * exp(-(dt) * (am + bm) );
        state3_ina =  ( inf ) + (state3_ina - inf ) * exp(-dt/tau);
        gating          = (state1_ina * state1_ina * state1_ina)  * state3_ina;
        return gating;
}