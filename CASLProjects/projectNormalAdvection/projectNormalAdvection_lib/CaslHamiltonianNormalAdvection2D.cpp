//
// Created by Frederic Gibou on 1/6/23.
//

#include "CaslHamiltonianNormalAdvection2D.h"
#include <cmath>

CaslHamiltonianNormalAdvection2D::CaslHamiltonianNormalAdvection2D(CaslGrid2D& grid, CaslArray2D<double> & vn) :
        _vn(vn), CaslHamiltonian2D(grid){}

double CaslHamiltonianNormalAdvection2D::H(const double phi_x, const double phi_y, const int  i, const int j, const double t) {
    return _vn(i,j) * sqrt( phi_x*phi_x + phi_y*phi_y );
}

double CaslHamiltonianNormalAdvection2D::maxAbsH1(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    if(phi_y_min*phi_y_max <= 0) return fabs(_vn(i,j)); // i.e. take abs_phi_y=0

    double abs_phi_x = std::max(fabs(phi_x_min), fabs(phi_x_max) );
    double abs_phi_y = std::min(fabs(phi_y_min), fabs(phi_y_max) );
    return fabs(_vn(i,j)) * abs_phi_x / sqrt(abs_phi_x * abs_phi_x + abs_phi_y * abs_phi_y);
}

double CaslHamiltonianNormalAdvection2D::maxAbsH2(const double phi_x_min, const double phi_x_max,
                                                  const double phi_y_min, const double phi_y_max,
                                                  const int i, const int j, const double t) {
    if(phi_x_min*phi_x_max <= 0) return fabs(_vn(i,j)); // i.e. take abs_phi_x=0

    double abs_phi_x = std::min(fabs(phi_x_min), fabs(phi_x_max) );
    double abs_phi_y = std::max(fabs(phi_y_min), fabs(phi_y_max) );
    return fabs(_vn(i,j)) * abs_phi_y / sqrt(abs_phi_x * abs_phi_x + abs_phi_y * abs_phi_y);
}