//
// Created by Frederic Gibou on 1/6/23.
//

#include "CaslHamiltonianReinitialization2D.h"
#include <cmath>

CaslHamiltonianReinitialization2D::CaslHamiltonianReinitialization2D(CaslGrid2D& grid, CaslArray2D<double> & signPhi0) :
    _signPhi0(signPhi0), CaslHamiltonian2D(grid){}

double CaslHamiltonianReinitialization2D::H(const double phi_x, const double phi_y, const int  i, const int j, const double t) {
    return _signPhi0(i,j) * ( sqrt(phi_x*phi_x + phi_y*phi_y) - 1 );
}

double CaslHamiltonianReinitialization2D::maxAbsH1(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    if(phi_y_min*phi_y_max <= 0) return 1; // i.e. take abs_phi_y=0

    double abs_phi_x = std::max(fabs(phi_x_min), fabs(phi_x_max) );
    double abs_phi_y = std::min(fabs(phi_y_min), fabs(phi_y_max) );
    return abs_phi_x / sqrt(abs_phi_x * abs_phi_x + abs_phi_y * abs_phi_y);
}

double CaslHamiltonianReinitialization2D::maxAbsH2(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    if(phi_x_min*phi_x_max <= 0) return 1; // i.e. take abs_phi_x=0

    double abs_phi_x = std::min(fabs(phi_x_min), fabs(phi_x_max) );
    double abs_phi_y = std::max(fabs(phi_y_min), fabs(phi_y_max) );
    return abs_phi_y / sqrt(abs_phi_x * abs_phi_x + abs_phi_y * abs_phi_y);
}