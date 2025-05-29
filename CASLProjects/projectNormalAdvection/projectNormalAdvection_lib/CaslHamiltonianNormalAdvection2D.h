//
// Created by Frederic Gibou on 1/6/23.
//

#ifndef CASL_HAMILTONIAN_NORMAL_ADVECTION_2D_H
#define CASL_HAMILTONIAN_NORMAL_ADVECTION_2D_H

#include "../../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../../CASLCommonLibrary/CaslArray2D.h"

class CaslHamiltonianNormalAdvection2D : public CaslHamiltonian2D {
public:
    CaslArray2D<double> & _vn;  // normal velocity

    explicit CaslHamiltonianNormalAdvection2D(CaslGrid2D& grid, CaslArray2D<double> & vn);

    double H(double phi_x, double phi_y, int  i, int j, double t) override;

    double maxAbsH1(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;

    double maxAbsH2(double phi_x_min, double phi_x_max,
                    double phi_y_min, double phi_y_max,
                    int i, int j, double t) override;
};


#endif //CASL_HAMILTONIAN_NORMAL_ADVECTION_2D_H
