//
// Created by faranak on 8/3/23.
//

/*
 * H:
     * v = phi!: The Cost-to-go function
     * if abs(vx) <= 2KuMax:
     * H = grad(v).F(z) - (1/4K^2)vx^2
     * else:
     * H = grad(v).F(z) + uMax^2 - |vx|uMax/K
 * maxAbsH1:
     * vx = max(abs(phi_x_min), abs(phi_x_max))
     * alpha_x = abs(fx) + max(uMax/K, 0.5/Ks^2 vx)
 * maxAbsH2:
     * alpha_y = abs(fy)
 */

#ifndef CASL_HAMILTONIAN_HH_MODEL_CPP
#define CASL_HAMILTONIAN_HH_MODEL_CPP

#include <cmath>
#include <vector>

#include "CaslHamiltonianHHModel.h"

CaslHamiltonianHHModel::CaslHamiltonianHHModel(CaslGrid2D& grid, CaslArray2D<double> & fx, CaslArray2D<double> & fy):
        _fx(fx), _fy(fy), CaslHamiltonian2D(grid){}

double CaslHamiltonianHHModel::H(const double phi_x, const double phi_y, const int  i, const int j, const double t) {

    if (fabs(phi_x) <= 2 * Ks * uMax) {
        return -_fx(i, j) * phi_x - _fy(i, j) * phi_y
        + ((.25 / pow(Ks, 2))) * pow(phi_x, 2);
    }
    else {
        return -_fx(i, j) * phi_x - _fy(i, j) * phi_y
        - pow(uMax, 2) + fabs(phi_x) * uMax / Ks;
    }
}

double CaslHamiltonianHHModel::maxAbsH1(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    double phi_x_1 = std::max(fabs(phi_x_min), fabs(phi_x_max));
    return fabs(_fx(i,j)) + std::max(uMax/Ks, (0.5 / pow(Ks, 2)) * phi_x_1);
}

double CaslHamiltonianHHModel::maxAbsH2(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    return fabs(_fy(i,j));
}

#endif // CASL_HAMILTONIAN_HH_MODEL_CPP