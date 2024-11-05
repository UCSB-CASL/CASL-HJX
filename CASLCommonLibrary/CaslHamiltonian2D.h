//
// Created by Frederic Gibou on 1/5/23.
//

#ifndef CASL_HAMILTONIAN_2D_H
#define CASL_HAMILTONIAN_2D_H

#include <iostream>
using namespace std;

#include "CaslGrid2D.h"

class CaslHamiltonian2D {
protected:
    CaslGrid2D& _grid;

    explicit CaslHamiltonian2D(CaslGrid2D &grid);

private:
    static void undefinedHamiltonianErrorMessage();

public:
    virtual double H(double phi_x, double phi_y, int  i, int j, double t);

    virtual double maxAbsH1(double phi_x_min, double phi_x_max,
                            double phi_y_min, double phi_y_max,
                            int i, int j, double t);

    virtual double maxAbsH2(double phi_x_min, double phi_x_max,
                            double phi_y_min, double phi_y_max,
                            int i, int j, double t);
};

#endif // CASL_HAMILTONIAN_2D_H
