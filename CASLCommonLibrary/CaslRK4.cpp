//
// Created by Faranak on 7/22/23.
//

#include "CaslRK4.h"
#include <vector>
#include <functional>
#include <iostream>

void rk4(const std::function<void(double, double, double, std::vector<double>&)>& f,
         std::vector<double>& y0, double t0, double tf, int NTimeS, double tyData[][3])
{
    double h = static_cast<double>((tf - t0) / (NTimeS - 1));
    std::vector<double> k1(2), k2(2), k3(2), k4(2);
    std::vector<double> yTemp = y0;
    std::vector<double> y_0sVec;
    std::vector<double> y_1sVec;
    std::vector<double> tsVec;

    for (int i = 0; i < NTimeS; i++)
    {
        // Calculate the next time step
        double next_t = t0 + h;

        // If the next time step is greater than tf, adjust h to match the remaining time
        if (next_t > tf) {
            h = tf - t0;
            next_t = tf;
        }

        // Update t0 to the calculated next_t
        t0 = next_t;

        y_0sVec.push_back(y0[0]);
        y_1sVec.push_back(y0[1]);
        tsVec.push_back(t0);

        f(t0, y0[0], y0[1], k1);
        f(t0 + h / 2.0, y0[0] + k1[0] * h / 2.0, y0[1] + k1[1] * h / 2.0, k2);
        f(t0 + h / 2.0, y0[0] + k2[0] * h / 2.0, y0[1] + k2[1] * h / 2.0, k3);
        f(t0 + h, y0[0] + k3[0] * h, y0[1] + k3[1] * h, k4);

        yTemp[0] = y0[0] + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * h / 6.0;
        yTemp[1] = y0[1] + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * h / 6.0;

        y0 = yTemp;
//        std::cout << y0[0] << std::endl;
    }

    int numRows = tsVec.size();
    for (int row = 0; row < numRows; row++)
    {
        tyData[row][0] = tsVec[row];
        tyData[row][1] = y_0sVec[row];
        tyData[row][2] = y_1sVec[row];
    }
}
