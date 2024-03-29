﻿#include <stdio.h>;
#include <conio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <windows.h>
#include <limits>
using namespace std;
clock_t clock(void);
void main()
{
    cout << "Liquid - Water, Gas - air" << endl;
    cout << endl;
    //int M1 = 1000; // число точек по z
    int M1 = 500; // число точек по z
    int M2 = 100; // число точек по r
    //описание динамических массивов
    cout << "Fill arrays" << endl;
    long double** T_g = new long double* [M1 + 1];
    for (int i = 0; i < M1 + 1; i++) T_g[i] = new long double[M2 + 1];
    cout << "Array T_g filled" << endl;
    long double** P_liq = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) P_liq[i] = new long double[M2 + 1];
    cout << "Array P_liq filled" << endl;
    long double** Q = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Q[i] = new long double[M2 + 1];
    cout << "Array Q filled" << endl;
    long double** W_A = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W_A[i] = new long double[M2 + 1];
    cout << "Array W_A filled" << endl;
    long double** W = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W[i] = new long double[M2 + 1];
    cout << "Array W filled" << endl;
    long double** Jakobian = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Jakobian[i] = new long double[M2 + 1];
    cout << "Array Jakobian filled" << endl;
    long double** NEW_Vz = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Vz[i] = new long double[M2 + 1];
    cout << "Array New_Vz filled" << endl;
    long double** NEW_Vr = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Vr[i] = new long double[M2 + 1];
    cout << "Array New_Vr filled" << endl;
    long double** NEW_P_liq = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_P_liq[i] = new long double[M2 + 1];
    cout << "Array New_PL filled" << endl;
    long double** NEW_Z = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Z[i] = new long double[M2 + 1];
    cout << "Array New_Z filled" << endl;
    long double** NEW_R = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_R[i] = new long double[M2 + 1];
    cout << "Array New_R filled" << endl;
    long double** Vz = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vz[i] = new long double[M2 + 1];
    cout << "Array Vz filled" << endl;
    long double** Vr = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vr[i] = new long double[M2 + 1];
    cout << "Array Vr filled" << endl;
    long double** Z = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Z[i] = new long double[M2 + 1];
    cout << "Array Z filled" << endl;
    long double** R = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) R[i] = new long double[M2 + 1];
    cout << "Array R filled" << endl;
    long double** Pe = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Pe[i] = new long double[M2 + 1];
    cout << "Array Pe filled" << endl;
    long double** NU = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NU[i] = new long double[M2 + 1];
    cout << "Array Nu filled" << endl;
    long double** A = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) A[i] = new long double[M2 + 1];
    cout << "Array A filled" << endl;
    long double** W_R = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W_R[i] = new long double[M2 + 1];
    cout << "Array w_r filled" << endl;
    long double** P_g = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) P_g[i] = new long double[M2 + 1];
    cout << "Array P_g filled" << endl;
    long double** alpha = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) alpha[i] = new long double[M2 + 1];
    cout << "Array Alfa filled" << endl;
    long double** NEW_A = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_A[i] = new long double[M2 + 1];
    cout << "Array New_A filled" << endl;
    long double** NEW_W_R = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_W_R[i] = new long double[M2 + 1];
    cout << "Array New_W_R filled" << endl;
    long double** NEW_P_g = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_P_g[i] = new long double[M2 + 1];
    cout << "Array New_P_G filled" << endl;
    long double** NEW_alpha = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_alpha[i] = new long double[M2 + 1];
    cout << "Array New_Alfa filled" << endl;
    long double** DJak = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) DJak[i] = new long double[M2 + 1];
    cout << "Array DJak filled" << endl;
    long double** DJ = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) DJ[i] = new long double[M2 + 1];
    cout << "Array DJ filled" << endl;
    long double** RO0 = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) RO0[i] = new long double[M2 + 1];
    cout << "Array Ro0 filled" << endl;
    long double** RO = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) RO[i] = new long double[M2 + 1];
    cout << "Array Ro filled" << endl;
    long double** Cb = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Cb[i] = new long double[M2 + 1];
    cout << "Array Cb filled" << endl;
    long double** Pz = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Pz[i] = new long double[2 * (M2 + 1)];
    cout << "Array Pz filled" << endl;
    long double** Vzz = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vzz[i] = new long double[2 * (M2 + 1)];
    cout << "Array Vzz filled" << endl;
    long double** Vrz = new long double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vrz[i] = new long double[2 * (M2 + 1)];
    cout << "Array Vrz filled" << endl;
    //конец описания динамических массивов
    cout << "Arrays are full" << endl;
    cout << endl;
    long double C_liq = 1500., P_0 = 1.e+5, hZ = 1.e-3, hR = 1.e-3, tau = 1.e-7, T_0 = 293., a0 = 1.e-3, c_g = 1003;
    long double V_liquid = 1e-6, ro_g = 1.2, ro_liq0 = 1000., yota = 1.4, lambda = 2.59 * 1e-2;
    //long double C_liq = 1500., P_0 = 1.e+5, hZ = 1.e-3, hR = 1.e-3, tau = 1.e-7, T_0 = 293., a0 = 1.e-3, c_g = 595;
    //long double V_liquid = 6e-6, ro_g = 5.06, ro_liq0 = 1115., yota = 1.14, lambda = 0.97 * 1e-2;
    long double K_g = lambda / (c_g * ro_g);
    cout << K_g << endl;
    for (int I = 0; I <= M1; I++)
        for (int J = 0; J <= M2; J++)
        {
            Vz[I][J] = 0;
            Vr[I][J] = 0;
            P_liq[I][J] = P_0;
            T_g[I][J] = T_0;
            Jakobian[I][J] = 1.0;
            RO[I][J] = ro_liq0;
            DJak[I][J] = 0.0;
            DJ[I][J] = 0.0;
            NEW_Vz[I][J] = 0;
            NEW_Vr[I][J] = 0;
            NEW_P_liq[I][J] = P_0;
            Z[I][J] = I * hZ;
            R[I][J] = J * hR;
            NEW_Z[I][J] = I * hZ;
            NEW_R[I][J] = J * hR;
            A[I][J] = a0;
            W_R[I][J] = 0.0;
            P_g[I][J] = P_0;
        }
    for (int I = 0; I <= M1; I++)
        for (int J = 0; J <= M2; J++)
        {
            //if (I >= 500 && I <= 600 && J >= 10)
            if (((I-250)*(I-250)+J*J) <= 2500)
                alpha[I][J] = 1e-2;
            else 
                alpha[I][J] = 1e-8;
            RO[I][J] = alpha[I][J] * ro_g + (1 - alpha[I][J]) * ro_liq0;
            RO0[I][J] = alpha[I][J] * ro_g + (1 - alpha[I][J]) * ro_liq0;
            Cb[I][J] = sqrt(yota * P_0 / (alpha[I][J] * ro_liq0));
            if (Cb[I][J] > C_liq) { Cb[I][J] = C_liq; }
        }
    long double T = 0.0;
    int K = 0;
    int numer = 1;
    int numer1 = 1;
    /*for (int J = 0; J <= M2; J++) {
        string Name = "PLR";
        string FileT = to_string(J);
        string TXT = "-1-100-200-300-400-500.txt";
        string FileName = Name + FileT + TXT;
        ofstream File;
        File.open(FileName, ios::ate);
    }*/
    ofstream Alpha;
    Alpha.open("alpha.txt");
    ofstream FileP_liq;
    FileP_liq.open("P.txt");
    while (K <= 10000)
    {
        if ((ceil(K / 100) - numer1) == 0) { cout << "\n" << K; numer1++; }
        long double T = K * tau;
        long double delta_P_liq = 3.e+5; // 0.5МПа
        T = K * tau;
        long double TZR1 = 2.E-4;
        long double TZ = TZR1 / 2.;
        long double ZN11 = 1. / (sqrt(4.0 * log(10.)));
        long double ZN1 = TZ * ZN11;
        if ((T < TZR1 / 2.0))
        {
            for (int J = 0; J <= M2; J++)
            {
                P_liq[0][J] = P_0 + delta_P_liq * exp(-((T - TZ) / (ZN1)) * ((T - TZ) / (ZN1)));
                //P_liq[0][J] = P_0 + delta_P_liq;
            }
        }
        else
        {
            for (int J = 0; J <= M2; J++) { P_liq[0][J] = P_0 + delta_P_liq; }
        }
        for (int I = 1; I <= M1 - 1; I++)
            for (int J = 1; J <= M2 - 1; J++)
            {      
                if (J == 1)
                {
                    Jakobian[I][J] = (Z[I + 1][J] - Z[I][J]) / hZ;
                    NEW_Vr[I][J] = 0;
                    NEW_Vz[I][J] = Vz[I][J] - (tau / (RO[I][J] * Jakobian[I][J])) * ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ));
                }
                else
                    if (J == M2)
                    {
                        Jakobian[I][J] = (R[I][J] / (J * hR)) * ((R[I][J] - R[I][J - 1]) * (Z[I + 1][J] - Z[I][J]) / (hR * hZ) -
                            (R[I + 1][J] - R[I][J]) * (Z[I][J] - Z[I][J - 1]) / (hR * hZ));
                        NEW_Vz[I][J] = Vz[I][J] + (tau / (RO[I][J] * Jakobian[I][J])) * (R[I][J] / ((J)*hR)) *
                            (((P_liq[I][J] - P_liq[I][J - 1]) / (hR)) * ((R[I + 1][J] - R[I][J]) / (hZ)) -
                                ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ)) * ((R[I][J] - R[I][J - 1]) / (hR)));
                        NEW_Vr[I][J] = 0;
                    }
                    else
                        if (J != 1 && J != M2)
                        {
                            Jakobian[I][J] = (R[I][J] / (J * hR)) * ((R[I][J + 1] - R[I][J]) * (Z[I + 1][J] - Z[I][J]) / (hR * hZ) -
                                (R[I + 1][J] - R[I][J]) * (Z[I][J + 1] - Z[I][J]) / (hR * hZ));
                            NEW_Vz[I][J] = Vz[I][J] + (tau / (RO[I][J] * Jakobian[I][J])) * (R[I][J] / ((J)*hR)) *
                                (((P_liq[I][J] - P_liq[I][J - 1]) / (hR)) * ((R[I + 1][J] - R[I][J]) / (hZ)) -
                                    ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ)) * ((R[I][J + 1] - R[I][J]) / (hR)));
                            NEW_Vr[I][J] = Vr[I][J] - (tau / (RO[I][J] * Jakobian[I][J])) * (R[I][J] / ((J)*hR)) *
                                (R[I][J] / ((J)*hR)) * (((P_liq[I][J] - P_liq[I][J - 1]) / (hR)) * ((Z[I + 1][J] - Z[I][J]) / (hZ)) -
                                    ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ)) * ((Z[I][J + 1] - Z[I][J]) / (hR)));
                        }
                NEW_Z[M1][J] = Z[M1][J] + tau * (Vz[M1][J]);
                NEW_R[M1][J] = R[M1][J] + tau * (Vr[M1][J]);
            }
        for (int I = 1; I < M1; I++)
            for (int J = 1; J <= M2; J++)
            {
                T_g[I][J] = (A[I][J] * A[I][J] * A[I][J] * P_g[I][J] * T_0) / (P_0 * a0 * a0 * a0);
                W_A[I][J] = (P_g[I][J] - P_liq[I][J]) / (ro_liq0 * C_liq * pow(alpha[I][J], 1 / 3));
                W[I][J] = W_R[I][J] + W_A[I][J];
                Pe[I][J] = 12. * (yota - 1.) * (T_0 * A[I][J] * fabs(W[I][J])) / (K_g * fabs(T_g[I][J] - (T_0 + 0.0000001)));
                if (Pe[I][J] > 100.) { NU[I][J] = sqrt(Pe[I][J]); }
                else { NU[I][J] = 10.; }
                Q[I][J] = NU[I][J] * lambda * (T_g[I][J] - T_0) / (2. * A[I][J]);
                NEW_A[I][J] = A[I][J] + (tau * W[I][J]);
                NEW_W_R[I][J] = W_R[I][J] + (tau * (((P_g[I][J] - P_liq[I][J]) / ro_liq0) -
                    1.5 * (W_R[I][J] * W_R[I][J]) - 4 * V_liquid * W_R[I][J] / A[I][J]) / A[I][J]);
                NEW_P_g[I][J] = P_g[I][J] - (tau * (3 * yota * P_g[I][J] * W[I][J]
                    + (3 * (yota - 1) * Q[I][J])) / A[I][J]);

                if (isinf(NEW_P_g[I][J])) {
                    cout << "I = " << I << "J = " << J << endl;
                    cout << "Nu = " << NU[I][J] << endl;
                    cout << "T_g = " << T_g[I][J] << endl;
                    cout << "A[I][J] * A[I][J] * A[I][J] * P_g[I][J] * T_0 = " << (A[I][J] * A[I][J] * A[I][J] * P_g[I][J] * T_0) << endl;
                    cout << "A[I][J] * A[I][J] * A[I][J] = " << A[I][J] * A[I][J] * A[I][J] << endl;
                    cout << "Vr = " << Vr[I][J] << endl;
                    cout << "Vz = " << Vz[I][J] << endl;
                    cout << "A = " << A[I][J] << endl;
                    cout << "Wr = " << W_A[I][J] << endl;
                    cout << "Wa = " << W_R[I][J] << endl;
                    cout << "W = " << W[I][J] << endl;
                    cout << "Q = " << Q[I][J] << endl;
                    cout << "P_g = " << P_g[I][J] << endl;
                    cout << "alpha = " << alpha[I][J] << endl;
                    cout << "DJ = " << DJ[I][J] << endl;
                    cout << "DJak = " << DJak[I][J] << endl;
                    cout << "Jakobian = " << Jakobian[I][J] << endl;
                    getchar();
                }

                NEW_Z[M1][J] = Z[M1][J] + tau * (Vz[M1][J]);
                NEW_R[I][M2] = R[I][M2] + tau * (Vr[I][M2]);
                Cb[I][J] = sqrt(yota * P_0 / (alpha[I][J] * ro_liq0));
                if (Cb[I][J] > C_liq) { Cb[I][J] = C_liq; }
                if (J == 1)
                {
                    Jakobian[I][J] = (Z[I + 1][J] - Z[I][J]) / (hZ);
                    DJ[I][J] = (NEW_Vz[I + 1][J] - NEW_Vz[I][J]) / (hZ)+(NEW_Z[I + 1][J] -
                        NEW_Z[I][J]) * (NEW_Vr[I][J + 1] - NEW_Vr[I][J]) / (hZ * hR) -
                        (NEW_Z[I][J + 1] - NEW_Z[I][J]) * (NEW_Vr[I + 1][J] - NEW_Vr[I][J]) / (hZ * hR);
                    DJak[I][J] = (NEW_Vr[I][J] * Jakobian[I][J] / (NEW_R[I][J])) + DJ[I][J];
                }
                else
                    if (J != 1 && J != M2)
                    {
                        DJ[I][J] = ((NEW_Vz[I + 1][J] - NEW_Vz[I][J]) * (NEW_R[I][J + 1] - NEW_R[I][J]) +
                            (NEW_Z[I + 1][J] - NEW_Z[I][J]) * (NEW_Vr[I][J + 1] - NEW_Vr[I][J]) -
                            (NEW_Vz[I][J + 1] - NEW_Vz[I][J]) * (NEW_R[I + 1][J] - NEW_R[I][J]) -
                            (NEW_Z[I][J + 1] - NEW_Z[I][J]) * (NEW_Vr[I + 1][J] - NEW_Vr[I][J])) / (hZ * hR);
                        DJak[I][J] = (NEW_Vr[I][J] * Jakobian[I][J] / (R[I][J])) + (R[I][J] / (J * hR)) * DJ[I][J];
                    }
                NEW_alpha[I][J] = alpha[I][J] + tau * ((3 * alpha[I][J] * W[I][J] / (A[I][J])) - DJak[I][J] * alpha[I][J] / Jakobian[I][J]);
                NEW_P_liq[I][J] = P_liq[I][J] + (tau * C_liq * C_liq * RO[I][J] * (3 * alpha[I][J] * W[I][J] / A[I][J] - DJak[I][J] *
                    (RO0[I][J] / (RO[I][J] * Jakobian[I][J] * Jakobian[I][J]) + alpha[I][J] / Jakobian[I][J])) / (1 - alpha[I][J]));
                // задание условия жесткой стенки на границе z = Lz
                //Vz[M2][J] = 0;
                // задание условия неотражения на границе z = Lz
                NEW_P_liq[M1 - 1][J] = (C_liq * ro_liq0) * (NEW_Vz[M1 - 1][J] - Vz[M1 - 1][J]) + P_liq[M1 - 1][J];
            }
        for (int I = 1; I <= M1; I++) {
            for (int J = 0; J <= M2; J++)
            {
                Vr[I][J] = NEW_Vr[I][J];
                Vz[I][J] = NEW_Vz[I][J];
                P_liq[I][J] = NEW_P_liq[I][J];
                A[I][J] = NEW_A[I][J];
                W_R[I][J] = NEW_W_R[I][J];
                P_g[I][J] = NEW_P_g[I][J];
                alpha[I][J] = NEW_alpha[I][J];
                Z[I][J] = NEW_Z[I][J];
                R[I][J] = NEW_R[I][J];
            }
        }
        for (int J = 0; J <= M2; J++)
        {
            FileP_liq << K << "\t" << P_liq[0][J] << endl;
        }
        for (int I = 0; I <= M1 - 1; I++) { P_liq[I][0] = P_liq[I][1]; }
        for (int I = 0; I <= M1 - 1; I++)
            for (int J = (0); J <= 2 * (M2); J++)
            {
                if (J <= (M2)) { Pz[I][J] = P_liq[I][(M2)-J]; }
                else if (J >= (M2)) { Pz[I][J] = P_liq[I][J - (M2)]; }
            }
        for (int J = 0; J <= M2; J++) {
            string Name = "PLR";
            string FileT = to_string(J);
            //string TXT = "-200-550-800.txt";
            string TXT = "-100-250-400.txt";
            string FileName = Name + FileT + TXT;
            ofstream File;
            File.open(FileName, ios::out | ios::app);
            //File << "\n" << K * tau << "\t" << (Pz[200][J]) << "\t" << (Pz[550][J])  << "\t" << (Pz[800][J]);
            File << "\n" << K * tau << "\t" << (Pz[100][J]) << "\t" << (Pz[250][J]) << "\t" << (Pz[400][J]);
        }
        for (int I = 0; I <= M1 - 1; I++) { Vz[I][0] = Vz[I][1]; }
        for (int I = 0; I <= M1 - 1; I++)
            for (int J = (0); J <= 2 * (M2 - 1); J++)
            {
                if (J <= (M2 - 1)) { Vzz[I][J] = Vz[I][(M2 - 1) - J]; }
                else if (J >= (M2 - 1)) { Vzz[I][J] = Vz[I][J - (M2 - 1)]; }
            }
        for (int I = 0; I <= M1 - 1; I++) { Vr[I][0] = Vr[I][1]; }
        for (int I = 0; I <= M1 - 1; I++)
            for (int J = (0); J <= 2 * (M2 - 1); J++)
            {
                if (J <= (M2 - 1)) { Vrz[I][J] = Vr[I][(M2 - 1) - J]; }
                else if (J >= (M2 - 1)) { Vrz[I][J] = Vr[I][J - (M2 - 1)]; }
            }
        if (K == 0) {
            for (int I = 0; I < M1; I++) {
                for (int J = 0; J < M2; J++)
                    Alpha << "\n" << I << "\t" << J << "\t" << alpha[I][J];
            }
        }
        if (K % 200 == 0)
        {
            string Name = "ZR";
            string FileT = to_string((K / 10000.0));
            string TXT = "msec.txt";
            string FileName = Name + FileT + TXT;
            ofstream file;
            file.open(FileName);
            for (int I = 1; I < M1; I++) {
                for (int J = 1; J <= 2 * (M2 - 1); J++)
                {
                    file << "\n" << I * hZ << "\t" << -((M2)-J)* hR << "\t" << Pz[I][J];
                }
            }
            string Name_ = "PL_Z_";
            string FileT_ = to_string((K / 10000.0));
            string TXT_ = "msec.txt";
            string FileName_ = Name_ + FileT_ + TXT_;
            ofstream file_;
            file_.open(FileName_);
            for (int I = 1; I < M1; I++) {
                file_ << "\n" << I * hZ << "\t" << Pz[I][1] << "\t" << Pz[I][M2 - 1];
            }
            numer++;
        }
        K++;
    }
    cout << "\n";
    cout << "\n";
    cout << clock() / CLOCKS_PER_SEC;
    //удаление динамических массивов
    delete[] T_g;
    delete[] P_liq;
    delete[] Q;
    delete[] W_A;
    delete[] W;
    delete[] Jakobian;
    delete[] NEW_Vz;
    delete[] NEW_Vr;
    delete[] NEW_P_liq;
    delete[] NEW_Z;
    delete[] NEW_R;
    delete[] Vz;
    delete[] Vr;
    delete[] Z;
    delete[] R;
    delete[] Pe;
    delete[] NU;
    delete[] A;
    delete[] W_R;
    delete[] P_g;
    delete[] alpha;
    delete[] NEW_A;
    delete[] NEW_W_R;
    delete[] NEW_P_g;
    delete[] NEW_alpha;
    delete[] DJak;
}
