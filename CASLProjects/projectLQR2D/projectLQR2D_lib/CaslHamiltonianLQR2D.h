//
// Created by Faranak Rajabi on 12/24/23.
//

/*
 * CaslHamiltonianLQR2D:
 *
 * Description:
 *   This class represents the Hamiltonian for a 2D Linear Quadratic Regulator (LQR) system.
 *   It extends the CaslHamiltonian2D class and provides specific methods for a 2D LQR system,
 *   including the calculation of the Hamiltonian, the maximum absolute value of the Hamiltonian
 *   gradients, and the differential equation for the evolution of the control matrix K.
 *
 * Members:
 *   - _fx1: CaslArray2D<double> - Array representing the first component of the system dynamics.
 *   - _fx2: CaslArray2D<double> - Array representing the second component of the system dynamics.
 *   - _grid: CaslGrid2D - Grid representing the discretized state space.
 *
 * Methods:
 *   - H(double phi_x, double phi_y, int i, int j, double t): Calculates the Hamiltonian for the 2D LQR system.
 *   - maxAbsH1(double phi_x_min, double phi_x_max, double phi_y_min, double phi_y_max, int i, int j, double t):
 *       Calculates the maximum absolute value of the first component of the Hamiltonian gradient.
 *   - maxAbsH2(double phi_x_min, double phi_x_max, double phi_y_min, double phi_y_max, int i, int j, double t):
 *       Calculates the maximum absolute value of the second component of the Hamiltonian gradient.
 *   - Other inherited methods from CaslHamiltonian2D.
 *
 * Usage:
 *   // Example usage of CaslHamiltonianLQR2D
 *   CaslGrid2D grid(initialize grid parameters);
 *   CaslArray2D<double> fx1(initialize fx1 values);
 *   CaslArray2D<double> fx2(initialize fx2 values);
 *   CaslHamiltonianLQR2D hamiltonian(grid, fx1, fx2);
 *   double hamiltonianValue = hamiltonian.H(provide arguments);
 *   // Perform other operations as needed.
 *
 * Author:
 *   Faranak Rajabi
 *   Created on: 12/24/23
 */

#ifndef CASL_HMAILTONIAN_LQR2D_H
#define CASL_HMAILTONIAN_LQR2D_H

#include "CaslHamiltonian2D.h"
#include "CaslArray2D.h"
#include <iostream>

class CaslHamiltonianLQR2D : public CaslHamiltonian2D {
public:
    CaslArray2D<double> &_fx1;
    CaslArray2D<double> &_fx2;

    explicit CaslHamiltonianLQR2D(CaslGrid2D &grid, CaslArray2D<double> &fx1, CaslArray2D<double> &fx2);

    double H(double phi_x, double phi_y, int i, int j, double t) override;

    double maxAbsH1(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;

    double maxAbsH2(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;
};

#include "CaslHamiltonianLQR2D.cpp"

#endif //CASL_HMAILTONIAN_LQR2D_H
