#include <iostream>
using namespace std;
double d(double x, double h) {//центральная разностная производная
    return ((x + h) - (x - h)) / (2 * h);
}
int main()
{
    //Constants (liq-water g-air)
    double p0_liq0 = 1000.0,
        v_liq = 10e-6,
        c_liq = 4.2,
        lambda_liq = 0.59,
        C_liq = 1500.0,
        T0 = 293.0,
        ksi_g0 = 0.1e-1,
        a0 = 0.1e-2,
        ro0_g = 1.2,
        lambda_g = 2.59 * 10e-2,
        yota = 1.4,
        c_g = 1.003,
        tau = 10e-8,
        h = 0.1e-2,
        ro_g = 1.2928,
        R = 8.31;
    //Formulas
    double
        k_g = lambda_g / (c_g * ro_g),
        Pe_g = 12 * (yota - 1) * (T0 / abs(T_g - T0)) * (a * abs(W) / k_g),
        p_g = ro0_g * R * T_g,
        p_liq = p_o * C_liq * C_liq * (p0_liq - p0_liq0),
        T_g = (p_g / p_0) * (a / a0) * (a / a0) * (a / a0) / T0;
    double Nu_g;
    if (Pe_g >= 100)
        Nu_g = sqrt(Pe_g);
    else
        Nu_g = 10;
    double
        q_g = Nu_g * lambda_g * (T_g - T0) / (2 * a),
        w_a = (p_g - p_liq) / (p0_liq * C_liq * ksi_g0),
        w = w_r + w_a,

}