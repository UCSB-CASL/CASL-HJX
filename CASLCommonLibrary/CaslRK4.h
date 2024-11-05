//
// Created by Faranak on 7/22/23.
//

#ifndef CASLUNIFORM_CASLRK4_H
#define CASLUNIFORM_CASLRK4_H

#include <functional>
#include <vector>

void rk4(const std::function<void(double, double, double, std::vector<double>&)>& f,
         std::vector<double>& y0, double t0, double tf, int NTimeS, double tyData[][3]);

#endif //CASLUNIFORM_CASLRK4_H

