//
// Created by Faranak on 7/22/23.
//

#include "CaslHHWholeSystemModel.h"

void vdot(double v[], double n[], double u, std::vector<double>& vDotsVec)
{
    double sum_i;
    double alpha[N][N] = {0.0};
    double sum_is[N] = {0.0};

    for (int i = 0; i < N; i++) {
        sum_i = 0;
        for (int j = 0; j < N; j++) {
            if (j != i)
            {
                Generator g(alpha_mu, alpha_sigma, 0, 1);
                double alpha_random = g();
                alpha[i][j] = alpha_random;
            }
            else
            {
                alpha[i][j] = 0;
            }
            sum_i += alpha[i][j] * (v[j] - v[i]);
        }
        sum_is[i] = sum_i;
        Generator g_eta(eta_mu, eta_sigma, 0, 1); // Changed variable name to g_eta to differentiate from the previous generator
        double eta = g_eta(); // Added () to call the operator function and generate the random value
        vDotsVec.push_back(fv(v[i], n[i]) + eta + (1.0 / N) * sum_is[i] + u);
    }
}

void ndot(double v[], double n[], std::vector<double>& nDotsVec)
{
    for (int i = 0; i < N; i++)
    {
        nDotsVec.push_back(fn(v[i], n[i]));
    }
}

