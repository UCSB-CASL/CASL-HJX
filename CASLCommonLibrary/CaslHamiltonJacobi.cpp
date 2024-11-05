//
// Created by Frederic Gibou on 1/1/23.
//

#ifndef CASL_HAMILTON_JACOBI_CPP
#define CASL_HAMILTON_JACOBI_CPP

#include "CaslHamiltonJacobi.h"

CaslHamiltonJacobi::CaslHamiltonJacobi()
{
    _firstDerivativeScheme = WENO5;
}

CaslHamiltonJacobi::CaslHamiltonJacobi(CaslOptionNumericalFirstDerivative firstDerivativeScheme)
{
    _firstDerivativeScheme = firstDerivativeScheme;
}

#endif // CASL_HAMILTON_JACOBI_CPP