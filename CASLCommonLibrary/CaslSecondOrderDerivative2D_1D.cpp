//
// Created by Faranak Rajabi on 2/8/24.
//

#ifndef CASL_SECOND_ORDER_DERIVATIVE_CPP
#define CASL_SECOND_ORDER_DERIVATIVE_CPP

#include "CaslSecondOrderDerivative2D_1D.h"

template <class T>
CaslSecondOrderDerivative2D<T>::CaslSecondOrderDerivative2D(
        CaslGrid2D &grid,
        T &dt,
        T &time,
        T &diffusionCoefficient,
        CaslOptionSecondOrderTermDirection secondOrderTermDirection,
        CaslOptionNumericalSecondDerivative secondDerivativeScheme,
        CaslOptionPaddingWith boundaryCondition,
        CaslHamiltonJacobi2D<T> &HJSolver,
        std::function<double(CaslGrid2D&, int, int, double, int)> inhomogeneousTermFunction
) :
        _grid(grid),
        _dt(dt),
        _time(time),
        _diffusionCoefficient(diffusionCoefficient),
        _inhomogeneousTermFunction(inhomogeneousTermFunction),
        _HJSolver(HJSolver),
        _secondOrderTermDirection(secondOrderTermDirection),
        _secondDerivativeScheme(secondDerivativeScheme),
        _boundaryCondition(boundaryCondition)
{}

template<class T>
CaslSecondOrderDerivative2D<T>::~CaslSecondOrderDerivative2D() = default;

template<class T>
double CaslSecondOrderDerivative2D<T>::findCFLHamiltonianDiffusion(const CaslArray2D<double> &un, const double _time) {
    // Equation     : ∂u/∂t + H(u, ∇u) - D ∇·(∇u) = F
    // CFL Condition: Δt {(|H1| / Δx) + (|H2| / Δy) + (2b / Δx²) + (2b / Δy²)} < CFL Number
    T dx = _grid.dx(), dy = _grid.dy();
    double CFLBasedOnDiffusion, CFLBasedOnHamiltonian;
    CFLBasedOnHamiltonian = _HJSolver.findCFL(un, _time); // (|H1| / Δx) + (|H2| / Δy)

    if (_secondDerivativeScheme == ForwardTimeCentralSpacing) {
        CFLBasedOnDiffusion = (2.0 * _diffusionCoefficient / (dx * dx)) + (2.0 * _diffusionCoefficient / (dy * dy));
    }
    else { // Since other schemes are unconditionally stable
        // For the pure diffusion case
        if (CFLBasedOnHamiltonian == 0) {
            CFLBasedOnDiffusion = 0.5 * std::min(_grid.dx(), _grid.dy());
        }
        else {
            CFLBasedOnDiffusion = 0.0;
        }
    }
    return (1.0 / (CFLBasedOnDiffusion + CFLBasedOnHamiltonian)); // One can multiply this out by the desired CFL number
}

template<class T>
void CaslSecondOrderDerivative2D<T>::forwardTimeCentralSpacing(const CaslArray2D<double> &un, CaslArray2D<double> &unp1) {
    int nX = _grid.nX(), nY = _grid.nY();

    CaslArray2D<double> Dxx(nX, nY), Dyy(nX, nY);

    computeDxx(un, Dxx);
    computeDyy(un, Dyy);

    // u_np1 = u_n + alpha * (uxx + uyy + uxy + uyx)
    // This is one iteration of the central differencing for the heat/ diffusion equation
    // This representation's only for generality of the solver

    //  Update the solution:
    double secondOrderTerm;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            secondOrderTerm = Dxx(i, j) + Dyy(i, j);
            unp1(i, j) = un(i, j) + _dt * _diffusionCoefficient * secondOrderTerm;
        }
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::forwardTimeCentralSpacingWithMixedDerivs(const CaslArray2D<double> &un, CaslArray2D<double> &unp1) {
    int nX = _grid.nX(), nY = _grid.nY();

    CaslArray2D<double> Dxx(nX, nY), Dyy(nX, nY), Dxy(nX, nY), Dyx(nX, nY);

    computeDxx(un, Dxx);
    computeDyy(un, Dyy);
    computeDxy(un, Dxy);
    computeDyx(un, Dyx);

    // u_np1 = u_n + alpha * (uxx + uyy + uxy + uyx)
    // This equation has no physical meaning when we have mixed-derivatives
    // This representation's only for generality of the solver

    //  Euler step: Update the solution:
    double secondOrderTerm;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            secondOrderTerm = Dxx(i, j) + Dyy(i, j) + Dxy(i, j) + Dyx(i, j);
            unp1(i, j) = un(i, j) + _dt * _diffusionCoefficient * secondOrderTerm;
        }
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::backwardTimeCentralSpacing(CaslArray2D<double>& un, CaslArray2D<double>& unp1, T tolerance, int iterations) {

    int nX = _grid.nX(), nY = _grid.nY();

    // Check on direction of diffusion term
    int dimension;
    heatEquationDimension(dimension);

    switch (dimension) {
        case 1: {  // Diffusion term: u_xx
            CaslArray2D<double> Dxx(un.nX(), un.nX()), RHS(un.nX(), un.nY()), RHS_copy(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);
                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;
                }
            }

            // Compute the coefficient matrix (Dxx)
            computeDxx(un, Dxx);

            // For the following BCs, the coefficient matrix will be tri-diagonal which is solvable using triDiagonal solvers
            if (_boundaryCondition == withLinearExtrapolation || _boundaryCondition == withConstantExtrapolation) {
                // symmetric tri-diagonal matrix on LHS * unp1
                triDiagonalMatrixSolver(Dxx, RHS, unp1);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // non-tri-diagonal matrix
                if (_diffusionCoefficient != 0.0) {
                    for (int j = 1; j<= RHS.nY(); j++) {
                        RHS(1, j) = 0.0;
                        RHS(nX, j) = 0.0;
                    }
                }

                invertMatrixAndMultiply(Dxx, RHS, unp1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation)
            {
                // non-tri-diagonal matrix
                invertMatrixAndMultiply(Dxx, RHS, unp1);
            }

            break;
        }

        case 2: {  // Diffusion term: u_yy
            CaslArray2D<double> Dyy(nY, nY), RHS(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);

                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;
                }
            }

            // Compute the coefficient matrix (Dyy)
            computeDyy(un, Dyy);

            // For the following BCs, the coefficient matrix will be tri-diagonal which is solvable using triDiagonal solvers
            if (_boundaryCondition == withLinearExtrapolation || _boundaryCondition == withConstantExtrapolation) {
                // symmetric tri-diagonal matrix on LHS * unp1
                triDiagonalMatrixSolver(Dyy, RHS, unp1);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // non-tri-diagonal matrix
                if (_diffusionCoefficient != 0.0) {
                    for (int j = 1; j<= RHS.nY(); j++) {
                        RHS(1, j) = 0.0;
                        RHS(nX, j) = 0.0;
                    }
                }
                invertMatrixAndMultiply(Dyy, RHS, unp1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation) {
                // non-tri-diagonal matrix
                invertMatrixAndMultiply(Dyy, RHS, unp1);
            }

            break;
        }

        case 12: { // 2D Diffusion term
            CaslArray2D<double> RHS(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);

                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;
                }
            }

            double dx    = _grid.dx();
            double dy    = _grid.dy();
            double alpha_x = _diffusionCoefficient * ( _dt / (dx * dx));
            double alpha_y = _diffusionCoefficient * ( _dt / (dy * dy));
            double beta  = 1.0 + 2.0 * (alpha_x + alpha_y);

            if (_boundaryCondition == withConstantExtrapolation) {
                DPMatrix2D<double> LaplacianSolver(nX, nY);
                LaplacianSolver.constantBoundary(beta, -alpha_x, -alpha_y);

                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);
            }

            if (_boundaryCondition == withLinearExtrapolation){
                DPMatrix2D<double> LaplacianSolver(nX, nY);
                LaplacianSolver.linearBoundary(beta, -alpha_x, -alpha_y);

                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // TODO: implement ILU with less sparse structure (condition number too high)
//                DPMatrixPeriodic2D<double> LaplacianSolver(nX, nY);
//                LaplacianSolver.periodicBoundary(beta, -alpha_x, -alpha_y);
//
//                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);

                std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid heat equation dimension for withPeriodicCondition! EXITING." << std::endl;
                exit(1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation) {
                // TODO: less sparse structure for non-adjacent cells
                std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid heat equation dimension for withQuadraticExtrapolation! EXITING." << std::endl;
                exit(1);
            }

            break;
        }

        default:
            std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid dimension for the heat equation! EXITING." << std::endl;
            exit(1);
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::crankNicolson(const CaslArray2D<double>& un, CaslArray2D<double>& unp1, T tolerance, int iterations){
    int nX = _grid.nX(), nY = _grid.nY();

    // Check on direction of diffusion term
    int dimension;
    heatEquationDimension(dimension);

    switch (dimension) {
        case 1: {  // Diffusion term: u_xx
            CaslArray2D<double> DxxExplicit(nX, nX), Dxx(nX, nX), RHS(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);


            // Halve diffusion coefficient temporarily for crank-nicolson
            _diffusionCoefficient = 0.5 * _diffusionCoefficient;

            // Compute explicit u_xx term
            computeDxxForwardTimeCentralSpacing(un, DxxExplicit);

            // Compute the coefficient matrix (Dxx)
            computeDxx(un, Dxx);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);
                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;

                    // Add explicit u_xx term
                    RHS(i, j) += _diffusionCoefficient * _dt * DxxExplicit(i, j);
                }
            }

            // Return diffusion coefficient to normal
            _diffusionCoefficient += _diffusionCoefficient;

            // For the following BCs, the coefficient matrix will be tri-diagonal which is solvable using triDiagonal solvers
            if (_boundaryCondition == withLinearExtrapolation || _boundaryCondition == withConstantExtrapolation) {
                // symmetric tri-diagonal matrix on LHS * unp1
                triDiagonalMatrixSolver(Dxx, RHS, unp1);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // non-tri-diagonal matrix
                if (_diffusionCoefficient != 0.0) {
                    for (int j = 1; j<= RHS.nY(); j++) {
                        RHS(1, j) = 0.0;
                        RHS(nX, j) = 0.0;
                    }
                }

                invertMatrixAndMultiply(Dxx, RHS, unp1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation)
            {
                // non-tri-diagonal matrix
                invertMatrixAndMultiply(Dxx, RHS, unp1);
            }

            break;
        }

        case 2: {  // Diffusion term: u_yy
            CaslArray2D<double> DyyExplicit(nY, nY), Dyy(nY, nY), RHS(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);

            // Halve diffusion coefficient temporarily for crank-nicolson
            _diffusionCoefficient = 0.5 * _diffusionCoefficient;

            // Compute explicit u_xx term
            computeDyyForwardTimeCentralSpacing(un, DyyExplicit);

            // Compute the coefficient matrix (Dyy)
            computeDyy(un, Dyy);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);

                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;

                    // Add explicit u_yy term
                    RHS(i, j) += _dt * _diffusionCoefficient * DyyExplicit(i, j);
                }
            }

            // Revert diffusion coefficient
            _diffusionCoefficient += _diffusionCoefficient;

            // For the following BCs, the coefficient matrix will be tri-diagonal which is solvable using triDiagonal solvers
            if (_boundaryCondition == withLinearExtrapolation || _boundaryCondition == withConstantExtrapolation) {
                // symmetric tri-diagonal matrix on LHS * unp1
                triDiagonalMatrixSolver(Dyy, RHS, unp1);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // non-tri-diagonal matrix
                if (_diffusionCoefficient != 0.0) {
                    for (int j = 1; j<= RHS.nY(); j++) {
                        RHS(1, j) = 0.0;
                        RHS(nX, j) = 0.0;
                    }
                }
                invertMatrixAndMultiply(Dyy, RHS, unp1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation) {
                // non-tri-diagonal matrix
                invertMatrixAndMultiply(Dyy, RHS, unp1);
            }

            break;
        }

        case 12: { // 2D Diffusion term
            CaslArray2D<double> RHS(nX, nY), Dxx(nX, nY), Dyy(nX, nY);
            un.fillPaddingPoints(_boundaryCondition);
            CaslArray2D<T> numericalHamiltonian = _HJSolver.computeNumericalHamiltonian(un);

            computeDxxForwardTimeCentralSpacing(un, Dxx);
            computeDyyForwardTimeCentralSpacing(un, Dyy);

            // Compute the right-hand side (RHS)
            for (int i = 1; i <= nX; i++) {
                for (int j = 1; j <= nY; j++) {
                    // For test cases that there is inhomogeneous BC at x = L(i = nX) || x = 0(i = 1):
                    // Change the whole RHS(i, j) with the function formula
                    RHS(i, j) = un(i, j);

                    // Add numerical hamiltonian
                    // Compute the numerical Hamiltonian at time step: n, solving explicitly
                    RHS(i, j) -= _dt * numericalHamiltonian(i, j);

                    // Compute the heat source term at time step: n
                    auto sourceTerm_ij = _inhomogeneousTermFunction(_grid, i, j, _time, dimension);
                    RHS(i, j) += _dt * sourceTerm_ij;

                    // Add explicit term
                    RHS(i, j) += _dt * 0.5 * _diffusionCoefficient * (Dxx(i, j) + Dyy(i, j));
                }
            }

            double dx    = _grid.dx();
            double dy    = _grid.dy();
            double alpha_x = 0.5 * _diffusionCoefficient * ( _dt / (dx * dx));
            double alpha_y = 0.5 * _diffusionCoefficient * ( _dt / (dy * dy));
            double beta  = 1.0 + 2.0 * (alpha_x + alpha_y);

            if (_boundaryCondition == withConstantExtrapolation) {
                DPMatrix<double> LaplacianSolver(nX, nY);
                LaplacianSolver.constantBoundary(beta, -alpha_x, -alpha_y);

                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);
            }

            if (_boundaryCondition == withLinearExtrapolation){
                DPMatrix<double> LaplacianSolver(nX, nY);
                LaplacianSolver.linearBoundary(beta, -alpha_x, -alpha_y);

                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);
            }

            if (_boundaryCondition == withPeriodicCondition) {
                // TODO: implement ILU with less sparse structure (condition number too high)
//                DPMatrixPeriodic2D<double> LaplacianSolver(nX, nY);
//                LaplacianSolver.periodicBoundary(beta, -alpha_x, -alpha_y);
//
//                LaplacianSolver.conjGrad(unp1, RHS, 1e-9, 1000);

                std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid heat equation dimension for withPeriodicCondition! EXITING." << std::endl;
                exit(1);
            }

            if (_boundaryCondition == withQuadraticExtrapolation) {
                // TODO: less sparse structure for non-adjacent cells
                std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid heat equation dimension for withQuadraticExtrapolation! EXITING." << std::endl;
                exit(1);
            }

            break;
        }

        default:
            std::cerr << "In CaslSecondOrderDerivative2D::backwardTimeCentralSpacing, Invalid dimension for the heat equation! EXITING." << std::endl;
            exit(1);
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxx(const CaslArray2D<double> &un, CaslArray2D<double> &Dxx) {
    if (_secondDerivativeScheme == ForwardTimeCentralSpacing)  {computeDxxForwardTimeCentralSpacing (un, Dxx); return;}
    if (_secondDerivativeScheme == CrankNicolson)              {computeDxxCrankNicolson             (un, Dxx); return;}
    if (_secondDerivativeScheme == BackwardTimeCentralSpacing) {computeDxxBackwardTimeCentralSpacing(un, Dxx); return;}

    std::cout << "In CaslSecondOrderDerivative2D::ComputeDxx, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxxForwardTimeCentralSpacing(const CaslArray2D<double> &un, CaslArray2D<double>& Dxx) {
    T dx = _grid.dx(); int nX = _grid.nX(), nY = _grid.nY();
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            Dxx(i, j) = ( un(i,j) - 2 * un(i, j) + un(i-1,j) ) / (pow(dx, 2));
        }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxxBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dxx) {
    // creating the tri-diagonal matrix for the heat term for implicit Euler
    int nX = _grid.nX();

    // In this case while the Dxx is basically the central differencing in the next time step, but for convenience,
    // We consider the Dxx to be the coefficient matrix directly
    CaslArray2D<double> triDiagonalMatrixOneD_x(un.nX(), un.nX());
    computeTriDiagonalMatrixOneDCase(triDiagonalMatrixOneD_x);
    Dxx = triDiagonalMatrixOneD_x;
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxxCrankNicolson(const CaslArray2D<double> &un, CaslArray2D<double> &Dxx) {

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyy(const CaslArray2D<double> &un, CaslArray2D<double> &Dyy) {
    if (_secondDerivativeScheme == ForwardTimeCentralSpacing) { computeDyyForwardTimeCentralSpacing(un, Dyy); return;}
    if (_secondDerivativeScheme == CrankNicolson  )     { computeDyyCrankNicolson(un, Dyy); return;}
    if (_secondDerivativeScheme == BackwardTimeCentralSpacing  )     { computeDyyBackwardTimeCentralSpacing(un, Dyy); return;}

    std::cout << "In CaslSecondOrderDerivative2D::ComputeDyy, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyyForwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dyy){
    T dy = _grid.dy(); int nX = _grid.nX(), nY = _grid.nY();
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            Dyy(i, j) = ( un(i,j) - 2 * un(i, j) + un(i-1,j) ) / (pow(dy, 2));
        }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyyBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dyy){
    // creating the tri-diagonal matrix for the heat term for implicit Euler
    int nY = _grid.nY();

    // In this case while the Dxx is basically the central differencing in the next time step, but for convenience,
    // We consider the Dxx to be the coefficient matrix directly
    CaslArray2D<double> triDiagonalMatrixOneD_y(nY, nY);
    computeTriDiagonalMatrixOneDCase(triDiagonalMatrixOneD_y);
    Dyy = triDiagonalMatrixOneD_y;
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyyCrankNicolson(const CaslArray2D<double>& un, CaslArray2D<double>& Dyy){

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxy(const CaslArray2D<double>& un, CaslArray2D<double>& Dxy){
    if (_secondDerivativeScheme == ForwardTimeCentralSpacing) { computeDxyForwardTimeCentralSpacing(un, Dxy); return;}
    if (_secondDerivativeScheme == CrankNicolson  )     { computeDxyCrankNicolson(un, Dxy); return;}
    if (_secondDerivativeScheme == BackwardTimeCentralSpacing  )     { computeDxyBackwardTimeCentralSpacing(un, Dxy); return;}

    std::cout << "In CaslSecondOrderDerivative2D::ComputeDxy, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxyForwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dxy){
    T dx = _grid.dx(), dy = _grid.dy();
    int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            Dxy(i, j) =   (un(i + 1, j + 1) - 2 * un(i, j) + un(i - 1, j - 1)) / (4.0 * dx * dy)
                        + (un(i + 1, j - 1) - 2 * un(i, j) + un(i - 1, j + 1)) / (4.0 * dx * dy);
        }
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxyBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dxy){

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDxyCrankNicolson(const CaslArray2D<double>& un, CaslArray2D<double>& Dxy){

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyx(const CaslArray2D<double>& un, CaslArray2D<double>& Dyx){
    if (_secondDerivativeScheme == ForwardTimeCentralSpacing) { computeDyxForwardTimeCentralSpacing(un, Dyx); return;}
    if (_secondDerivativeScheme == CrankNicolson  )     { computeDyxCrankNicolson(un, Dyx); return;}
    if (_secondDerivativeScheme == BackwardTimeCentralSpacing  )     { computeDyxBackwardTimeCentralSpacing(un, Dyx); return;}

    std::cout << "In CaslSecondOrderDerivative2D::ComputeDyx, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyxForwardTimeCentralSpacing(const CaslArray2D<double> &un,
                                                                         CaslArray2D<double> &Dyx) {
    T dx = _grid.dx(), dy = _grid.dy();
    int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            Dyx(i, j) =   (un(i + 1, j + 1) - 2 * un(i, j) + un(i - 1, j - 1)) / (4.0 * dx * dy)
                        + (un(i - 1, j + 1) - 2 * un(i, j) + un(i + 1, j - 1)) / (4.0 * dx * dy);
        }
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyxBackwardTimeCentralSpacing(const CaslArray2D<double>& un, CaslArray2D<double>& Dyx) {

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeDyxCrankNicolson(const CaslArray2D<double>& un, CaslArray2D<double>& Dyx) {

}

/* Helper private member functions */
template<class T>
void CaslSecondOrderDerivative2D<T>::heatEquationDimension(int &dimension) {
    if (_secondOrderTermDirection == X) {
        dimension = 1;
    }
    else if (_secondOrderTermDirection == Y) {
        dimension = 2;
    }
    else if (_secondOrderTermDirection == XY) {
        dimension = 12;
    }
    else {
        std::cerr << "Error: Invalid value for _secondOrderTermDirection. Expected X, Y, or XY." << std::endl;
    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeTriDiagonalMatrixOneDCase(CaslArray2D<double>& triDiagonalMatrix_1D) {
    if (triDiagonalMatrix_1D.nX() != triDiagonalMatrix_1D.nY()) {
        std::cerr << "In CaslSecondOrderDerivative2D::computeTriDiagonalMatrixOneDCase, the input matrix is NOT a squared matrix! EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    int dimension;
    heatEquationDimension(dimension);

    if (dimension == 1) {
        if (triDiagonalMatrix_1D.nX() != nX) {
            std::cerr << "In CaslSecondOrderDerivative2D::computeTriDiagonalMatrixOneDCase, the tri-diagonal matrix for the second derivative in x-direction size doesn't match with grid size! EXITING." << std::endl;
            exit(1);
        }

        double dx    = _grid.dx();
        double alpha = _diffusionCoefficient * (_dt / (dx * dx));
        double beta  = 1.0 + 2.0 * alpha;

        for (int i = 1; i <= nX; ++i) {
            triDiagonalMatrix_1D(i, i) = beta;
        }

        for (int i = 1; i < nX; ++i) {
            triDiagonalMatrix_1D(i       , i + 1) = -alpha;
            triDiagonalMatrix_1D(i + 1, i    ) = -alpha;
        }

        // So far the above case was true for the case that we have dirichlet boundary condition as a general case
        // Change the matrix based on the other boundary condition cases
        if (_diffusionCoefficient != 0) {
            if (_boundaryCondition == withLinearExtrapolation) {
                triDiagonalMatrix_1D(1, 1      ) = 1.0;
                triDiagonalMatrix_1D(1, 2      ) = 0.0;
                triDiagonalMatrix_1D(nX, nX    ) = 1.0;
                triDiagonalMatrix_1D(nX, nX - 1) = 0.0;
            }

            else if (_boundaryCondition == withConstantExtrapolation) {
                triDiagonalMatrix_1D(1, 1  ) = 1.0 + alpha;
                triDiagonalMatrix_1D(nX, nX) = 1.0 + alpha;
            }

            else if (_boundaryCondition == withPeriodicCondition) {
                triDiagonalMatrix_1D(1, 1 ) =  1.0;
                triDiagonalMatrix_1D(1, 2 ) =  0.0;
                triDiagonalMatrix_1D(1, nX) = -1.0;

                triDiagonalMatrix_1D(nX, 1     ) =  0.0;
                triDiagonalMatrix_1D(nX, 2     ) = -0.5;
                triDiagonalMatrix_1D(nX, nX    ) =  1.0;
                triDiagonalMatrix_1D(nX, nX - 1) = -0.5;
            }

            else if (_boundaryCondition == withQuadraticExtrapolation) {
                // Refilling matrix for first row:
                triDiagonalMatrix_1D(1, 1) = 1.0 - alpha;
                triDiagonalMatrix_1D(1, 2) = 2.0 * alpha;
                triDiagonalMatrix_1D(1, 3) =     - alpha;

                // Refilling matrix for last row:
                triDiagonalMatrix_1D(nX, nX - 2) = -     alpha;
                triDiagonalMatrix_1D(nX, nX - 1) = 2.0 * alpha;
                triDiagonalMatrix_1D(nX, nX    ) = 1.0 - alpha;
            }
        }
    }

    else if (dimension == 2) {
        if (triDiagonalMatrix_1D.nX() != nY) {
            std::cerr << "Error: the tri-diagonal matrix for the second derivative in y-direction size doesn't match with grid size!" << std::endl;
            exit(1);
        }

        double dy = _grid.dy();
        double alpha = _diffusionCoefficient * (_dt / (dy * dy));
        double beta = 1.0 + 2.0 * alpha;

        for (int j = 1; j <= nY; ++j) {
            triDiagonalMatrix_1D(j, j) = beta;
        }

        for (int j = 1; j < nY; ++j) {
            triDiagonalMatrix_1D(j, j + 1) = -alpha;
            triDiagonalMatrix_1D(j + 1, j) = -alpha;
        }

        if (_boundaryCondition == withLinearExtrapolation) {
            triDiagonalMatrix_1D(1, 1) = 1.0;
            triDiagonalMatrix_1D(1, 2) = 0.0;
            triDiagonalMatrix_1D(nY, nY) = 1.0;
            triDiagonalMatrix_1D(nY, nY - 1) = 0.0;
        }

        if (_boundaryCondition == withPeriodicCondition) {
            triDiagonalMatrix_1D(1, 1 ) =  1.0;
            triDiagonalMatrix_1D(1, 2 ) =  0.0;
            triDiagonalMatrix_1D(1, nY) = -1.0;

            triDiagonalMatrix_1D(nY, 1     ) =  0.0;
            triDiagonalMatrix_1D(nY, 2     ) = -0.5;
            triDiagonalMatrix_1D(nY, nY    ) =  1.0;
            triDiagonalMatrix_1D(nY, nY - 1) = -0.5;
        }

        if (_boundaryCondition == withConstantExtrapolation) {
            triDiagonalMatrix_1D(1, 1) = 1.0 + alpha;
            triDiagonalMatrix_1D(nY, nY) = 1.0 + alpha;
        }

        if (_boundaryCondition == withQuadraticExtrapolation) {
            // Refilling matrix for first row:
            triDiagonalMatrix_1D(1, 1) = 1.0 - alpha;
            triDiagonalMatrix_1D(1, 2) = 2.0 * alpha;
            triDiagonalMatrix_1D(1, 3) =     - alpha;

            // Refilling matrix for last row:
            triDiagonalMatrix_1D(nY, nY - 2) = -     alpha;
            triDiagonalMatrix_1D(nY, nY - 1) = 2.0 * alpha;
            triDiagonalMatrix_1D(nY, nY    ) = 1.0 - alpha;
        }

    }
}

template<class T>
void CaslSecondOrderDerivative2D<T>::triDiagonalMatrixSolver(CaslArray2D<double>& A, CaslArray2D<double>& RHS, CaslArray2D<double>& X) {
    // Check if matrices have compatible sizes
    if (A.nX() != A.nY() || RHS.nX() != A.nX() || X.nX() != A.nX()) {
        std::cerr << "IN CaslCaslSecondOrderDerivative2D:: Incompatible matrix sizes! EXITING." << std::endl;
        exit(1); // or throw an exception
    }

    // Determines the size
    int nX = RHS.nX();
    int nY = RHS.nY();

    // Forward elimination
    for (int i = 2; i <= nX; ++i) {
        double w = A(i, i - 1) / A(i - 1, i - 1);
        A(i, i) = A(i, i) - w * A(i - 1, i);

        for (int j = 1; j <= nY; ++j) {
            RHS(i, j) = RHS(i, j) - w * RHS(i - 1, j);
        }
    }

    // Backward substitution
    for (int j = 1; j <= nY; ++j) {
        X(nX, j) = RHS(nX, j) / A(nX, nX);

        for (int i = nX - 1; i >= 1; --i) {
            X(i, j) = (RHS(i, j) - A(i, i + 1) * X(i + 1, j)) / A(i, i);
        }
    }
}


template<class T>
void CaslSecondOrderDerivative2D<T>::invertMatrix(CaslArray2D<T>& A, CaslArray2D<T>& inverseA) {

    int n = A.nX();

    CaslArray2D<T> I(n, n);
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=n; j++) {
            I(i,j) = (i==j) ? 1 : 0;
        }
    }

    // Perform Gaussian elimination
    for(int i=1; i<n; i++) {
        int maxRow = i;
        T maxVal = abs(A(i,i));
        for(int k=i+1; k<=n; k++) {
            if(abs(A(k,i)) > maxVal) {
                maxRow = k;
                maxVal = abs(A(k,i));
            }
        }
        // Swap rows
        if(maxRow != i) {
            for(int k=1; k<=n; k++) {
                swap(A(maxRow,k), A(i,k));
                swap(I(maxRow,k), I(i,k));
            }
        }
        // Elimination
        for(int k=i+1; k<=n; k++) {
            T c = -A(k,i)/A(i,i);
            for(int j=1; j<=n; j++) {
                A(k,j) += c * A(i,j);
                I(k,j) += c * I(i,j);
            }
        }
    }

    // Back substitution
    for(int i=n; i>=1; i--) {
        for(int k=i-1; k>=1; k--) {
            T c = -A(k,i)/A(i,i);
            for(int j=1; j<=n; j++) {
                A(k,j) += c * A(i,j);
                I(k,j) += c * I(i,j);
            }
        }
        T c = 1/A(i,i);
        for(int j=1; j<=n; j++) {
            A(i,j) *= c;
            I(i,j) *= c;
        }
    }

    inverseA = I;

}


template<class T>
void CaslSecondOrderDerivative2D<T>::invertMatrixAndMultiply(CaslArray2D<T>& A, CaslArray2D<T>& RHS, CaslArray2D<T>& solution) {

    int n = A.nX();

    // Create identity matrix
    CaslArray2D<T> I(n, n);
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=n; j++) {
            I(i,j) = (i==j) ? 1 : 0;
        }
    }

    // Gaussian elimination
    for(int i=1; i<n; i++) {
        int maxRow = i;
        T maxVal = abs(A(i,i));
        for(int k=i+1; k<=n; k++) {
            if(abs(A(k,i)) > maxVal) {
                maxRow = k;
                maxVal = abs(A(k,i));
            }
        }
        // Swap rows
        if(maxRow != i) {
            for(int k=1; k<=n; k++) {
                swap(A(maxRow,k), A(i,k));
                swap(I(maxRow,k), I(i,k));
            }
        }

        // Elimination
        for(int k=i+1; k<=n; k++) {
            T c = -A(k,i)/A(i,i);
            for(int j=1; j<=n; j++) {
                A(k,j) += c * A(i,j);
                I(k,j) += c * I(i,j);
            }
        }

    }

    // Back substitution
    for(int i=n; i>=1; i--) {
        for(int k=i-1; k>=1; k--) {
            T c = -A(k,i)/A(i,i);
            for(int j=1; j<=n; j++) {
                A(k,j) += c * A(i,j);
                I(k,j) += c * I(i,j);
            }
        }
        T c = 1/A(i,i);
        for(int j=1; j<=n; j++) {
            A(i,j) *= c;
            I(i,j) *= c;
        }
    }

    // Multiply inverse matrix by RHS
    for(int i=1; i<=n; i++) {
        for(int j=1; j<=RHS.nY(); j++) {
            T sum = 0;
            for(int k=1; k<=n; k++) {
                sum += I(i,k) * RHS(k,j);
            }
            solution(i,j) = sum;
        }
    }

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeBoundaryPolyCoeffs(CaslArray2D<double>& coeffsVec, int polynomialDegree, double x) {
    int numberOfPoints = polynomialDegree + 1;
    int dimension;
    heatEquationDimension(dimension);

    if (coeffsVec.nX() != numberOfPoints || coeffsVec.nY() != 1) {
        std::cerr << "Error: The coefficient vector has an incorrect size. Expected size = polynomial degree + 1. Exiting." << std::endl;
        exit(1);
    }

    int nX = _grid.nX();
    int nY = _grid.nY();

    if (dimension == 1) {
        if (numberOfPoints > nX) {
            std::cerr << "Error: Polynomial degree is out of bounds. Expected polynomial degree <= grid.nX() - 1. Exiting." << std::endl;
            exit(1);
        }

        if (x < 1) {  // There is a ghost point to the left
            computeLagrangeCoefficients(coeffsVec, polynomialDegree, x, 1, numberOfPoints);
        }
        else if (x > nX) {  // There is a ghost point to the right
            computeLagrangeCoefficients(coeffsVec, polynomialDegree, x, nX - numberOfPoints + 1, nX);
        }
    }

    else if (dimension == 2) {
        if (numberOfPoints > nY) {
            std::cerr << "Error: Polynomial degree is out of bounds. Expected polynomial degree <= grid.nY() - 1. Exiting." << std::endl;
            exit(1);
        }

        if (x < 1) {  // There is a ghost point to the left
            computeLagrangeCoefficients(coeffsVec, polynomialDegree, x, 1, numberOfPoints);
        }
        else if (x > nY) {  // There is a ghost point to the right
            computeLagrangeCoefficients(coeffsVec, polynomialDegree, x, nY - numberOfPoints + 1, nY);
        }
    }

    // Need generalization
//    else if (_secondOrderTermDirection == XY) {
//
//    }

}

template<class T>
void CaslSecondOrderDerivative2D<T>::computeLagrangeCoefficients(CaslArray2D<double> &coeffsVec, double x, int start,
                                                                 int end, int i) {
    double term;

    for (int i = start; i <= end; i++) {
        term = 1;
        for (int j = start; j <= end; j++) {
            if (i != j) {
                term *= (x - j) / (i - j);
            }
        }
        // The order of this vector's elements is increasing based on the node index on the grid
        // example, i = 0: coeffsVec(1, 1) = 3 (for u_1), coeffsVec(2, 1) = -3 (for u_2), coeffsVec(3, 1) = 1 (for u_1)
        // example, i = nX + 1: coeffsVec(1, 1) = 1 (for u_nX-2), coeffsVec(2, 1) = -3 (for u_nX-1), coeffsVec(3, 1) = 3 (for u_nX)
        coeffsVec(i - start + 1, 1) = term;
    }
}

#endif // CASL_SECOND_ORDER_DERIVATIVE_CPP