#include <stdio.h>;
#include <conio.h>
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
using namespace std;
clock_t clock(void);
void main()
{
    cout << "Liquid - Water, Gas - air" << endl;
    cout << endl;
    int M1 = 300; // число точек по z
    int M2 = 30; // число точек по r
    //описание динамических массивов
    cout << "Fill arrays" << endl;
    double** T_g = new double* [M1 + 1];
    for (int i = 0; i < M1 + 1; i++) T_g[i] = new double[M2 + 1];
    cout << "Array T_g filled" << endl;
    double** P_liq = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) P_liq[i] = new double[M2 + 1];
    cout << "Array P_liq filled" << endl;
    double** Q = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Q[i] = new double[M2 + 1];
    cout << "Array Q filled" << endl;
    double** W_A = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W_A[i] = new double[M2 + 1];
    cout << "Array W_A filled" << endl;
    double** W = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W[i] = new double[M2 + 1];
    cout << "Array W filled" << endl;
    double** Jakobian = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Jakobian[i] = new double[M2 + 1];
    cout << "Array Jakobian filled" << endl;
    double** NEW_Vz = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Vz[i] = new double[M2 + 1];
    cout << "Array New_Vz filled" << endl;
    double** NEW_Vr = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Vr[i] = new double[M2 + 1];
    cout << "Array New_Vr filled" << endl;
    double** NEW_P_liq = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_P_liq[i] = new double[M2 + 1];
    cout << "Array New_PL filled" << endl;
    double** NEW_Z = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_Z[i] = new double[M2 + 1];
    cout << "Array New_Z filled" << endl;
    double** NEW_R = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_R[i] = new double[M2 + 1];
    cout << "Array New_R filled" << endl;
    double** Vz = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vz[i] = new double[M2 + 1];
    cout << "Array Vz filled" << endl;
    double** Vr = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vr[i] = new double[M2 + 1];
    cout << "Array Vr filled" << endl;
    double** Z = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Z[i] = new double[M2 + 1];
    cout << "Array Z filled" << endl;
    double** R = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) R[i] = new double[M2 + 1];
    cout << "Array R filled" << endl;
    double** Pe = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Pe[i] = new double[M2 + 1];
    cout << "Array Pe filled" << endl;
    double** NU = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NU[i] = new double[M2 + 1];
    cout << "Array Nu filled" << endl;
    double** A = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) A[i] = new double[M2 + 1];
    cout << "Array A filled" << endl;
    double** W_R = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) W_R[i] = new double[M2 + 1];
    cout << "Array w_r filled" << endl;
    double** P_g = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) P_g[i] = new double[M2 + 1];
    cout << "Array P_g filled" << endl;
    double** alpha = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) alpha[i] = new double[M2 + 1];
    cout << "Array Alfa filled" << endl;
    double** NEW_A = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_A[i] = new double[M2 + 1];
    cout << "Array New_A filled" << endl;
    double** NEW_W_R = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_W_R[i] = new double[M2 + 1];
    cout << "Array New_W_R filled" << endl;
    double** NEW_P_g = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_P_g[i] = new double[M2 + 1];
    cout << "Array New_P_G filled" << endl;
    double** NEW_alpha = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) NEW_alpha[i] = new double[M2 + 1];
    cout << "Array New_Alfa filled" << endl;
    double** DJak = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) DJak[i] = new double[M2 + 1];
    cout << "Array DJak filled" << endl;
    double** DJ = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) DJ[i] = new double[M2 + 1];
    cout << "Array DJ filled" << endl;
    double** RO0 = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) RO0[i] = new double[M2 + 1];
    cout << "Array Ro0 filled" << endl;
    double** RO = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) RO[i] = new double[M2 + 1];
    cout << "Array Ro filled" << endl;
    double** Cb = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Cb[i] = new double[M2 + 1];
    cout << "Array Cb filled" << endl;
    double** Pz = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Pz[i] = new double[2 * (M2 + 1)];
    cout << "Array Pz filled" << endl;
    double** Vzz = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vzz[i] = new double[2 * (M2 + 1)];
    cout << "Array Vzz filled" << endl;
    double** Vrz = new double* [M1 + 1]; //
    for (int i = 0; i < M1 + 1; i++) Vrz[i] = new double[2 * (M2 + 1)];
    cout << "Array Vrz filled" << endl;
    //конец описания динамических массивов
    cout << "Arrays are full" << endl;
    cout << endl;
    double P_0 = 1.e+5, hZ = 1.e-3, hR = 1.e-3, tau = 1.e-7, T_0 = 293., a0 = 1.e-3, c_g = 1003;
    double V_liquid = 1e-6, ro_g = 1.2, ro_liq0 = 1000., yota = 1.4, lambda = 2.59 * 1e-2;
    double K_g = lambda / (c_g * ro_g);
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
            alpha[I][J] = 1.e-2;
            RO[I][J] = alpha[I][J] * ro_g + (1 - alpha[I][J]) * ro_liq0;
            RO0[I][J] = alpha[I][J] * ro_g + (1 - alpha[I][J]) * ro_liq0;
            Cb[I][J] = sqrt(yota * P_0 / (alpha[I][J] * ro_liq0));
        }
    for (int I = 0; I <= M1; I++) { alpha[I][0] = alpha[I][1]; }
    double T = 0.0;
    int K = 0;
    int numer = 1;
    int numer1 = 1;
    ofstream FileDeltaP_liq;
    FileDeltaP_liq.open("PL-1-50-100-150-200-250-300.txt");
    ofstream FileP_liq;
    FileP_liq.open("P.txt");
    while (K <= 50000)
    {
        if ((ceil(K / 100) - numer1) == 0) { cout << "\n" << K; numer1++; }
        double T = K * tau;
        double delta_P_liq = 5.e+5; // 0.5МПа
        T = K * tau;
        double TZR1 = 2.E-4;
        double TZ = TZR1 / 2.;
        double ZN11 = 1. / (sqrt(4.0 * log(10.)));
        double ZN1 = TZ * ZN11;
        if ((T < TZR1))
        {
            for (int J = 0; J <= M2; J++)
            {
                P_liq[0][J] = P_0 + delta_P_liq * exp(-((T - TZ) / (ZN1)) * ((T - TZ) / (ZN1)));
            }
        }
        else
        {
            for (int J = 0; J <= M2; J++) { P_liq[0][J] = P_0 + delta_P_liq; }
        }
        for (int I = 1; I <= M1 - 1; I++)
            for (int J = 1; J <= M2 - 1; J++)
            {      
                Jakobian[I][J] = (R[I][J] / (J * hR)) * ((R[I][J + 1] - R[I][J]) * (Z[I + 1][J] - Z[I][J]) / (hR * hZ) -
                    (R[I + 1][J] - R[I][J]) * (Z[I][J + 1] - Z[I][J]) / (hR * hZ));
                NEW_Vz[I][J] = Vz[I][J] + (tau / (RO[I][J] * Jakobian[I][J])) * (R[I][J] / ((J)*hR)) *
                    (((P_liq[I][J] - P_liq[I][J - 1]) / (hR)) * ((R[I + 1][J] - R[I][J]) / (hZ)) -
                        ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ)) * ((R[I][J + 1] - R[I][J]) / (hR)));
                NEW_Vr[I][J] = Vr[I][J] - (tau / (RO[I][J] * Jakobian[I][J])) * (R[I][J] / ((J)*hR)) *
                    (R[I][J] / ((J)*hR)) * (((P_liq[I][J] - P_liq[I][J - 1]) / (hR)) * ((Z[I + 1][J] - Z[I][J]) / (hZ)) -
                        ((P_liq[I][J] - P_liq[I - 1][J]) / (hZ)) * ((Z[I][J + 1] - Z[I][J]) / (hR)));
                NEW_Z[M1][J] = Z[M1][J] + tau * (Vz[M1][J]);
                NEW_R[M1][J] = R[M1][J] + tau * (Vr[M1][J]);
            }
        for (int I = 1; I <= M1 - 1; I++)
            for (int J = 1; J <= M2 - 1; J++)
            {
                T_g[I][J] = (A[I][J] * A[I][J] * A[I][J] * P_g[I][J] * T_0) / (a0 * a0 * a0 * P_0);
                W_A[I][J] = (P_g[I][J] - P_liq[I][J]) / (ro_liq0 * Cb[I][J] * pow(alpha[I][J], 1 / 3));
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
                NEW_Z[M1][J] = Z[M1][J] + tau * (Vz[M1][J]);
                NEW_R[I][M2] = R[I][M2] + tau * (Vr[I][M2]);
                Cb[I][J] = sqrt(yota * P_0 / (alpha[I][J] * ro_liq0));
                Jakobian[I][J] = (Z[I + 1][J] - Z[I][J]) / (hZ);
                DJ[I][J] = (NEW_Vz[I + 1][J] - NEW_Vz[I][J]) / (hZ)+(NEW_Z[I + 1][J] - NEW_Z[I][J]) * (NEW_Vr[I][J + 1] -
                    NEW_Vr[I][J]) / (hZ * hR) - (NEW_Z[I][J + 1] - NEW_Z[I][J]) * (NEW_Vr[I + 1][J] - NEW_Vr[I][J]) / (hZ * hR);
                DJak[I][J] = (NEW_Vr[I][J] * Jakobian[I][J] / (NEW_R[I][J])) + DJ[I][J];
                NEW_P_liq[I][J] = P_liq[I][J] + (tau * Cb[I][J] * Cb[I][J] * RO[I][J] * (3 * alpha[I][J] * W[I][J] / A[I][J] - DJak[I][J] *
                    (RO0[I][J] / (RO[I][J] * Jakobian[I][J] * Jakobian[I][J]) + alpha[I][J] / Jakobian[I][J])) / (1 - alpha[I][J]));
                // задание условия жесткой стенки на границе z = Lz
                Vz[M2-1][J] = 0;
            }
        for (int I = 1; I <= M1; I++) {
            for (int J = 1; J <= M2; J++)
            {
                Vr[I][J] = NEW_Vr[I][J];
                Vz[I][J] = NEW_Vz[I][J];
                P_liq[I][J] = NEW_P_liq[I][J];
                A[I][J] = NEW_A[I][J];
                W_R[I][J] = NEW_W_R[I][J];
                P_g[I][J] = NEW_P_g[I][J];
                Z[I][J] = NEW_Z[I][J];
                R[I][J] = NEW_R[I][J];
            }
        }
        for (int I = 0; I <= M1 - 1; I++) { P_liq[I][0] = P_liq[I][1]; }
        for (int I = 0; I <= M1 - 1; I++)
            for (int J = (0); J <= 2 * (M2); J++)
            {
                if (J <= (M2)) { Pz[I][J] = P_liq[I][(M2)-J]; }
                else if (J >= (M2)) { Pz[I][J] = P_liq[I][J - (M2)]; }
            }
        FileDeltaP_liq << "\n" << K * tau << "\t" << (Pz[1][15] - P_0) << "\t" << (Pz[50][15] - P_0) << "\t" << (Pz[100][15] - P_0) << "\t" << (Pz[150][15] - P_0) << "\t" << (Pz[200][15] - P_0) << "\t" << (Pz[250][15] - P_0) << "\t" << (Pz[299][15] - P_0);
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
        if (K % 10000 == 0)
        {
            string Name = "ZR";
            string FileT = to_string((int)(K / 10000));
            string TXT = "msec.txt";
            string FileName = Name + FileT + TXT;
            ofstream file;
            file.open(FileName);
            for (int I = 0; I < M1; I++) {
                for (int J = 1; J <= M2; J++)
                {
                    file << "\n" << I * hZ << "\t" << -((M2)-J) * hR << "\t" << Pz[I][J] - P_0;
                    FileP_liq << K << "\t" << P_liq[0][J] << endl;
                }
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
