//
// Created by Faranak Rajabi on 2/8/24.
//

/**
 * @brief A template class for solving second-order partial differential equations in 2D space.
 *
 * This class provides functionality to solve PDEs with both second-order and first-order derivatives in 2D space,
 * using various numerical schemes.
 * It supports different boundary conditions and numerical methods for solving the underlying partial differential equations.
 *
 * @details Solution steps for Hamiltonian2D-diffusion1D type problems:
 * 1. Solve the Hamiltonian term for the n+1 time step, take it as a temporal mid-point for solution in time (we call it u*)
 * 2. Use the solution from step one (u*) for the Laplacian term solution (as part of the RHS).
 *
 * Equation:
 *    u_t + H(x, y, u_x, u_y) - D * u_xx = F (I)
 *    or
 *    u_t + H(x, y, u_x, u_y) - D * u_yy = F (II)
 *
 * Example:
 *    Using LLF for Hamiltonian term, and BTCS for Laplacian term, for (I):
 *    [-alpha * u(i+1, j) + (1 + 2 * alpha) u(i, j) -alpha * u(i-1, j)] at (time_step_n+1) =
 *    [-dt * (u*(i, j) + F(i, j)) + u(i, j)] at (time_step_n)
 *
 *    Matrix form:
 *    [A * U] at (time_step_n+1) = [U - dt * H + F] at (time_step_n)
 *
 *  Where:
 *    alpha = D * dt / dx^2
 *    u*(i, j) = u_n(i, j) - _dt * numericalHamiltonian;
 *
 * @tparam T The type parameter for the template.
 */

#ifndef CASL_SECOND_ORDER_DERIVATIVE_H
#define CASL_SECOND_ORDER_DERIVATIVE_H

#include <functional>
#include <fstream>

#include "CaslArray2D.h"
#include "CaslGrid2D.h"
#include "CaslHamiltonian2D.h"
#include "CaslHamiltonJacobi2D.h"

template <class T> class CaslSecondOrderDerivative2D {
private:
    CaslGrid2D                                  & _grid;                      // grid
    T                                           & _dt;                        // time step
    T                                           & _time;                      // time
    T                                           & _diffusionCoefficient;      // Diffusion coefficient
    CaslHamiltonJacobi2D<T>                     & _HJSolver;                  //
    std::function<double(CaslGrid2D&, int, int, double, int)> & _inhomogeneousTermFunction;
    CaslOptionNumericalSecondDerivative           _secondDerivativeScheme;    // Upwind or ENO2 or ENO3 or WENO5
    CaslOptionSecondOrderTermDirection            _secondOrderTermDirection;  // heat equation dimension(x, y, xy)
    CaslOptionPaddingWith                         _boundaryCondition;         // boundary condition for the second-order derivative term

public:
    /**
      * @brief Constructor for CaslSecondOrderDerivative2D.
      *
      * @param grid                     The grid for the computation.
      * @param dt                       The time step.
      * @param time                     The current time.
      * @param diffusionCoefficient     The diffusion coefficient.
      * @param inhomogeneousTerm        The inhomogeneous term.
      * @param HJSolver                 The Hamilton-Jacobi solver.
      * @param secondOrderTermDirection The direction for the second-order term (x, y, xy).
      * @param secondDerivativeScheme   The numerical scheme for second-order derivatives.
      * @param boundaryCondition        The boundary condition for the second-order derivative term.
     */
      CaslSecondOrderDerivative2D(CaslGrid2D & grid, T & dt, T & time, T & diffusionCoefficient,
            CaslOptionSecondOrderTermDirection secondOrderTermDirection, CaslOptionNumericalSecondDerivative secondDerivativeScheme, CaslOptionPaddingWith boundaryCondition,
            CaslHamiltonJacobi2D<T> & HJSolver, std::function<double(CaslGrid2D&, int, int, double, int)> inhomogeneousTermFunction);

     ~CaslSecondOrderDerivative2D();
    double findCFLHamiltonianDiffusion           (const CaslArray2D<double>& un, double time);               // dt based on the CFL Condition for advection-diffusion tyoe equatiosn
    void forwardTimeCentralSpacing               (const CaslArray2D<double>& un, CaslArray2D<double>& unp1); // FTCS
    void forwardTimeCentralSpacingWithMixedDerivs(const CaslArray2D<double>& un, CaslArray2D<double>& unp1); // FTCS with mixed derivatives
    void backwardTimeCentralSpacing              (CaslArray2D<double> &un, CaslArray2D<double> &unp1);       // BTCS
    void crankNicolson                           (const CaslArray2D<double> &un, CaslArray2D<double> &unp1); // Crank Nicolson (average of FTCS and BTCS)

private:
    void heatEquationDimension               (int & dimension);

    void computeBoundaryPolyCoeffs           (CaslArray2D<double>& coeffsVec, int polynomialDegree, double x);
    void computeLagrangeCoefficients         (CaslArray2D<double>& coeffsVec, double x, int start, int end, int i);

    void computeTriDiagonalMatrixOneDCase    (CaslArray2D<double> &triDiagonalMatrix_1D);
    void triDiagonalMatrixSolver             (CaslArray2D<double>& A, CaslArray2D<double>& RHS, CaslArray2D<double>& X);
    void invertMatrix                        (CaslArray2D<T>& A, CaslArray2D<T>& inverseA);
    void invertMatrixAndMultiply             (CaslArray2D<T>& A, CaslArray2D<T>& RHS, CaslArray2D<T>& solution);

    void computeDxx                          (const CaslArray2D<double>& un, CaslArray2D<double>& Dxx);
    void computeDxxForwardTimeCentralSpacing (const CaslArray2D<double>& un, CaslArray2D<double>& Dxx);
    void computeDxxBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dxx);
    void computeDxxCrankNicolson             (const CaslArray2D<double>& un, CaslArray2D<double>& Dxx);

    void computeDyy                          (const CaslArray2D<double>& un, CaslArray2D<double>& Dyy);
    void computeDyyForwardTimeCentralSpacing (const CaslArray2D<double>& un, CaslArray2D<double>& Dyy);
    void computeDyyBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dyy);
    void computeDyyCrankNicolson             (const CaslArray2D<double>& un, CaslArray2D<double>& Dyy);

    void computeDxy                          (const CaslArray2D<double>& un, CaslArray2D<double>& Dxy);
    void computeDxyForwardTimeCentralSpacing (const CaslArray2D<double>& un, CaslArray2D<double>& Dxy);
    void computeDxyBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dxy);
    void computeDxyCrankNicolson             (const CaslArray2D<double>& un, CaslArray2D<double>& Dxy);

    void computeDyx                          (const CaslArray2D<double>& un, CaslArray2D<double>& Dyx);
    void computeDyxForwardTimeCentralSpacing (const CaslArray2D<double>& un, CaslArray2D<double>& Dyx);
    void computeDyxBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dyx);
    void computeDyxCrankNicolson             (const CaslArray2D<double>& un, CaslArray2D<double>& Dyx);

};

#include "CaslSecondOrderDerivative2D_1D.cpp"

#endif //CASL_SECOND_ORDER_DERIVATIVE_H
