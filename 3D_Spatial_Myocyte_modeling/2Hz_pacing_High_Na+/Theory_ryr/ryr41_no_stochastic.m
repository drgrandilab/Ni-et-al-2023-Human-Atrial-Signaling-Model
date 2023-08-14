function [Po,K] =  ryr41_no_stochastic(Caj, CaSRj, Ku,Kb,Kcp,KK)


hill_coef = 18; %23 -> 10
rho = 5000 * 1.0 / (1.0 + (KK / CaSRj)^hill_coef);
MM = (sqrt(1 + 8 * rho * 400) - 1) / (4 * rho * 400);

k12 = Ku * Caj^3 / (Kcp^3 + Caj^3) + 0.000001;
k43 = Kb * Caj^3 / (Kcp^3 + Caj^3) + 0.000001;
k14 = MM / 0.5 ;
k21 = 1 / 2;
k23 = k14;
k41 = 1 / 125;
k34 = 1 / 0.3;
k32 = k41 * k12  * k34 / k43 / k21 ;

A = [1, 1, 1, 1;
-k12-k14, k21, 0, k41;
k12, -k21-k23, k32, 0;
0, k23, -k32-k34, k43];

b = [41;0;0;0];
K = [k12 k14 k21 k23 k32 k34 k43 k41];
RyR = A\b;

Po = RyR(2)+RyR(3);







