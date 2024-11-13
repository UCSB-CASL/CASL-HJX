//
// Created by Faranak on 8/3/23.
//

#ifndef HJPRACTICE_CASLHAMILTONIANHHMODEL_H
#define HJPRACTICE_CASLHAMILTONIANHHMODEL_H

#include "CaslHamiltonian2D.h"
#include "CaslArray2D.h"
#include "CaslHHWholeSystemModel.h"
#include <iostream>

class CaslHamiltonianHHModel : public CaslHamiltonian2D {
public:
    // Need to upgrade
    CaslArray2D<double> &  _fx; // (1/K)fv
    CaslArray2D<double> &  _fy; // fn
    const double uMax = Ib;
    const double Ks = N;

    explicit CaslHamiltonianHHModel(CaslGrid2D& grid, CaslArray2D<double> & fx, CaslArray2D<double> & fy);

    double H(double phi_x, double phi_y, int  i, int j, double t) override;

    double maxAbsH1(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;

    double maxAbsH2(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;
};

#endif //HJPRACTICE_CASLHAMILTONIANHHMODEL_H
