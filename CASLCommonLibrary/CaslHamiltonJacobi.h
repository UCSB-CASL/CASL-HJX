//
// Created by Frederic Gibou on 1/1/23.
//

#ifndef CASL_HAMILTON_JACOBI_H
#define CASL_HAMILTON_JACOBI_H

#include "CaslOptions.h"

class CaslHamiltonJacobi {
protected:
    CaslOptionNumericalFirstDerivative _firstDerivativeScheme;

    // inline so compiled once. Otherwise, put all the implementation in the .h and do not write the .cpp:
    inline CaslHamiltonJacobi();
    inline explicit CaslHamiltonJacobi(CaslOptionNumericalFirstDerivative firstDerivativeScheme);
};

#include "CaslHamiltonJacobi.cpp"

#endif // CASL_HAMILTON_JACOBI_H
