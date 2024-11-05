//
// Created by Frederic Gibou on 1/5/23.
//

#ifndef CASL_HAMILTONIAN2D_CPP
#define CASL_HAMILTONIAN2D_CPP

#include "CaslHamiltonian2D.h"

CaslHamiltonian2D::CaslHamiltonian2D(CaslGrid2D &grid) : _grid(grid) {}

void CaslHamiltonian2D::undefinedHamiltonianErrorMessage() {
    cout << "CASL ERROR in Hamiltonian2D - this hamiltonian is not defined." <<
         "\n Exiting." << endl;
    exit(1);
}

double CaslHamiltonian2D::H(const double phi_x, const double phi_y, const int  i, const int j, const double t) {
    undefinedHamiltonianErrorMessage();
    return 0;
}

double CaslHamiltonian2D::maxAbsH1(const double phi_x_min, const double phi_x_max,
                                   const double phi_y_min, const double phi_y_max,
                                   const int i, const int j, const double t) {
    undefinedHamiltonianErrorMessage();
    return 0;
}

double CaslHamiltonian2D::maxAbsH2(const double phi_x_min, const double phi_x_max,
                                   const double phi_y_min, const double phi_y_max,
                                   const int i, const int j, const double t) {
    undefinedHamiltonianErrorMessage();
    return 0;
}

#endif // CASL_HAMILTONIAN2D_CPP