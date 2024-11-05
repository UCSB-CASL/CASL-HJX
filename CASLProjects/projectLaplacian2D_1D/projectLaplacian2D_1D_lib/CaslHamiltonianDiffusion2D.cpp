//
// Created by Faranak Rajabi on 2/18/24.
//

#include "CaslHamiltonianDiffusion2D.h"

CaslHamiltonianDiffusion2D::CaslHamiltonianDiffusion2D(CaslGrid2D& grid, CaslArray2D<double> & c1, CaslArray2D<double> & c2) :
        _c1(c1), _c2(c2), CaslHamiltonian2D(grid) {}

double CaslHamiltonianDiffusion2D::H(const double phi_x, const double phi_y, const int i, const int j, const double t) {
    return  -_c1(i,j) * _grid.x(i) * phi_x - _c2(i,j) * _grid.y(j) * phi_y;
}

double CaslHamiltonianDiffusion2D::maxAbsH1(const double phi_x_min, const double phi_x_max, const double phi_y_min,
                                            const double phi_y_max, const int i, const int j, const double t) {
    return fabs(_c1(i, j));
}

double CaslHamiltonianDiffusion2D::maxAbsH2(const double phi_x_min, const double phi_x_max, const double phi_y_min,
                                            const double phi_y_max, const int i, const int j, const double t) {
    return fabs(_c2(i, j));
}