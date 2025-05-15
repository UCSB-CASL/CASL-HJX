// // CaslHamiltonianLQR2D.cpp
// // Created by Faranak Rajabi on 12/24/23
//
// #ifndef CASL_HAMILTONIAN_LQR2D_CPP
// #define CASL_HAMILTONIAN_LQR2D_CPP
//
// #include "CaslHamiltonianLQR2D.h"
// #include <cmath>
//
// // Constructor: Initialize 2D Hamiltonian with system dynamics fx1, fx2
// CaslHamiltonianLQR2D::CaslHamiltonianLQR2D(CaslGrid2D& grid,
//                                            CaslArray2D<double>& fx1,
//                                            CaslArray2D<double>& fx2)
//     : _fx1(fx1), _fx2(fx2), CaslHamiltonian2D(grid) {}
//
// // Hamiltonian function H(phi_x, phi_y, i, j, t)
// // For the 2D LQR system: min_u { 0.5 * (x1^2 + x2^2) + ∇φ · f(x,u) }
// // where f = [x2; u] and u* = -φ_y
//
// double CaslHamiltonianLQR2D::H(double phi_x, double phi_y, int i, int j, double t) {
//     double x1 = _grid.x(i);
//     double x2 = _grid.y(j);
//
//     double H = 0.5 * (x1 * x1 + x2 * x2)         // State cost
//              + x2 * phi_x                        // phi_x * dx/dt = phi_x * x2
//              - 0.5 * phi_y * phi_y;              // min_u { phi_y * u + 0.5 * u^2 } = -0.5 * phi_y^2
//
//     return -1.0 * H; // Return negative of Hamiltonian per HJ equation convention
// }
//
// // Maximum magnitude of ∂H/∂φ_x = fx1 = x2
// double CaslHamiltonianLQR2D::maxAbsH1(double phi_x_min, double phi_x_max,
//                                       double phi_y_min, double phi_y_max,
//                                       int i, int j, double t) {
//     return std::abs(_fx1(i, j));  // = |x2|
// }
//
// // Maximum magnitude of ∂H/∂φ_y = -φ_y
// double CaslHamiltonianLQR2D::maxAbsH2(double phi_x_min, double phi_x_max,
//                                       double phi_y_min, double phi_y_max,
//                                       int i, int j, double t) {
//     return std::max(std::abs(phi_y_min), std::abs(phi_y_max));
// }
//
// #endif // CASL_HAMILTONIAN_LQR2D_CPP
//
// Created by Faranak Rajabi on 12/24/23.
//

// CaslHamiltonianLQR4D.cpp

#ifndef CASL_HAMILTONIAN_LQR2D_CPP
#define CASL_HAMILTONIAN_LQR2D_CPP

#include "CaslHamiltonianLQR2D.h"
#include <cmath>

CaslHamiltonianLQR2D::CaslHamiltonianLQR2D(CaslGrid2D &grid, CaslArray2D<double> &fx1, CaslArray2D<double> &fx2) :
        _fx1(fx1), _fx2(fx2),
        CaslHamiltonian2D(grid) {}

double CaslHamiltonianLQR2D::H(double phi_x, double phi_y, int i, int j, double t) {
    // min H on u = (0.5) * x' * Q * x- Jx3^2/2 + Jx1 * (x3) - Jx4^2/2 + Jx2 * (x4)
    // fx1 = x2, fx2 = 0.
    /*
     * Hamiltonian H:
        min H on u: (0.5) * (x1^2 + x2^2) + x2 * Jx1 + - Jx2^2 / 2
        X2*conj(Jx1) - (Jx2*conj(Jx2))/2 + (X1*conj(X1))/2 + (X2*conj(X2))/2

        dHdJx:
          X2
        -Jx2

        J* Analytical:
        X1*((K1_1*conj(X1))/2 + (K2_1*conj(X2))/2) + X2*((K1_2*conj(X1))/2 + (K2_2*conj(X2))/2)

        u* Analytical:
        - K2_1*X1 - K2_2*X2

        Kdot:
        [      -K1_2*K2_1,       K1_1 - K1_2*K2_2]
        [K1_1 - K2_1*K2_2, - K2_2^2 + K1_2 + K2_1]
     */
    double Hamiltonian = - (pow(phi_y, 2) / 2.0) + phi_x * _fx1(i, j) +
                         (1.0 / 2.0) * (pow(_grid.x(i), 2) + pow(_grid.y(j), 2));

    return (-1.0) * Hamiltonian;
}

double CaslHamiltonianLQR2D::maxAbsH1(double phi_x_min, double phi_x_max,
                                      double phi_y_min, double phi_y_max,
                                      int i, int j, double t) {
    // dHdJx:
    //  X2 (fx1)
    //-Jx2
    return fabs(_fx1(i, j));
}

double CaslHamiltonianLQR2D::maxAbsH2(double phi_x_min, double phi_x_max,
                                      double phi_y_min, double phi_y_max,
                                      int i, int j, double t) {
    double phi_y_1 = std::max(fabs(phi_y_min), fabs(phi_y_max));
    return phi_y_1;
}

#endif //CASL_HAMILTONIAN_LQR4D_CPP
