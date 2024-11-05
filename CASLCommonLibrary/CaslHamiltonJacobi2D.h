//
// Created by Frederic Gibou on 1/5/23.
//

#ifndef CASL_HAMILTON_JACOBI2D_H
#define CASL_HAMILTON_JACOBI2D_H

#include <fstream>

#include "CaslGrid2D.h"
#include "CaslArray2D.h"
#include "CaslHamiltonJacobi.h"
#include "CaslHamiltonian2D.h"

template<class T>  class CaslHamiltonJacobi2D : public CaslHamiltonJacobi {
private:
    CaslGrid2D                         & _grid;                    // grid
    T                                  & _dt;                      // time step
    T                                  & _time;                    // time
    CaslOptionNumericalFirstDerivative   _firstDerivativeScheme;   // Upwind or ENO2 or ENO3 or WENO5
    CaslHamiltonian2D                  & _hamiltonian;             // Hamiltonian
public:
    CaslHamiltonJacobi2D(CaslGrid2D & grid, CaslHamiltonian2D & hamiltonian, T & dt, T & time, CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5);

    virtual ~CaslHamiltonJacobi2D();

    CaslArray2D<T> computeNumericalHamiltonian(const CaslArray2D<double>& un);
    // Before this function is called, one should have defined
    // the values of un in the paddings:
    void eulerStep(const CaslArray2D<double>& un, CaslArray2D<double>& unp1);
    // Before this function is called, one should have defined
    // the values of un in the paddings:
    double findCFL(const CaslArray2D<double>& un, const double time);
    // Before this function is called, one should have defined
    // the values of un in the paddings:
    void computeOptimalControl(const CaslArray2D<double>& un, int K, double uMax, CaslArray2D<double>& uStar);

private:
    void computeDxMinusDxPlus      (const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus);
    void computeDxMinusDxPlusUpwind(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus);
    void computeDxMinusDxPlusENO2  (const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus);
    void computeDxMinusDxPlusENO3  (const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus);
    void computeDxMinusDxPlusWENO5 (const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus);

    void computeDyMinusDyPlus      (const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus);
    void computeDyMinusDyPlusUpwind(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus);
    void computeDyMinusDyPlusENO2  (const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus);
    void computeDyMinusDyPlusENO3  (const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus);
    void computeDyMinusDyPlusWENO5 (const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus);

};

#include "CaslHamiltonJacobi2D.cpp"

#endif // CASL_HAMILTON_JACOBI2D_H