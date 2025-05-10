//
// Created by Faranak Rajabi on 2/18/24.
//

#ifndef CASL_HAMILTONIAN_BURGERS_2D_H
#define CASL_HAMILTONIAN_BURGERS_2D_H

#include "../../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../../CASLCommonLibrary/CaslArray2D.h"

class CaslHamiltonianBurgers2D : public CaslHamiltonian2D {
public:
    CaslArray2D<double> _c1;
    CaslArray2D<double> _c2;
    const double _constant = 0.0; // Modify if needed

    explicit CaslHamiltonianBurgers2D(CaslGrid2D& grid, CaslArray2D<double> & c1, CaslArray2D<double> & c2);

    double H(const double phi_x, const double phi_y, const int  i, const int j, const double t) override;

    double maxAbsH1(const double phi_x_min, const double phi_x_max,
                    const double phi_y_min, const double phi_y_max,
                    const int i, const int j, const double t) override;

    double maxAbsH2(const double phi_x_min, const double phi_x_max,
                    const double phi_y_min, const double phi_y_max,
                    const int i, const int j, const double t) override;
};


#endif //CASL_HAMILTONIAN_BURGERS_2D_H
