//
// Created by faranak on 8/3/23.
//

/*
 * H:
 *   v = phi : cost-to-go function
 *
 *   if abs(vx) <= 2*K*uMax:
 *       H = grad(v)·F(z) - (1/(4K^2)) vx^2
 *   else:
 *       H = grad(v)·F(z) + uMax^2 - |vx|uMax/K
 *
 * maxAbsH1:
 *   vx = max(abs(phi_x_min), abs(phi_x_max))
 *   alpha_x = abs(fx) + max(uMax/K, 0.5/K^2 * vx)
 *
 * maxAbsH2:
 *   alpha_y = abs(fy)
 */

#ifndef CASL_HAMILTONIAN_HH_MODEL_CPP
#define CASL_HAMILTONIAN_HH_MODEL_CPP

#include <cmath>
#include <vector>
#include <algorithm>

#include "CaslHamiltonianHHModel.h"

CaslHamiltonianHHModel::CaslHamiltonianHHModel(CaslGrid2D& grid,
                                               CaslArray2D<double>& fx,
                                               CaslArray2D<double>& fy)
    : _fx(fx), _fy(fy), CaslHamiltonian2D(grid) {}

double CaslHamiltonianHHModel::H(const double phi_x,
                                 const double phi_y,
                                 const int i,
                                 const int j,
                                 const double t) {
    (void)t;

    const double fx = _fx(i, j);
    const double fy = _fy(i, j);

    if (std::fabs(phi_x) <= 2.0 * Ks * uMax) {
        return -fx * phi_x - fy * phi_y
               + (0.25 / (static_cast<double>(Ks) * static_cast<double>(Ks))) * phi_x * phi_x;
    } else {
        return -fx * phi_x - fy * phi_y
               - uMax * uMax
               + std::fabs(phi_x) * uMax / Ks;
    }
}

double CaslHamiltonianHHModel::maxAbsH1(const double phi_x_min,
                                        const double phi_x_max,
                                        const double phi_y_min,
                                        const double phi_y_max,
                                        const int i,
                                        const int j,
                                        const double t) {
    (void)phi_y_min;
    (void)phi_y_max;
    (void)t;

    const double phi_x_1 = std::max(std::fabs(phi_x_min), std::fabs(phi_x_max));
    return std::fabs(_fx(i, j))
           + std::max(uMax / Ks, (0.5 / (static_cast<double>(Ks) * static_cast<double>(Ks))) * phi_x_1);
}

double CaslHamiltonianHHModel::maxAbsH2(const double phi_x_min,
                                        const double phi_x_max,
                                        const double phi_y_min,
                                        const double phi_y_max,
                                        const int i,
                                        const int j,
                                        const double t) {
    (void)phi_x_min;
    (void)phi_x_max;
    (void)phi_y_min;
    (void)phi_y_max;
    (void)t;

    return std::fabs(_fy(i, j));
}

#endif // CASL_HAMILTONIAN_HH_MODEL_CPP