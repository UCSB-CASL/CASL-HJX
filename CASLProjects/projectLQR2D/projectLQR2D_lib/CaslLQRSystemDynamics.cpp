//
// Created by Faranak Rajabi on 12/24/23.
//

#ifndef CASL_LQR_SYSTEM_DYNAMICS_CPP
#define CASL_LQR_SYSTEM_DYNAMICS_CPP

#include "CaslLQRSystemDynamics.h"
#include <cmath>

// Constructor
CaslLQRSystemDynamics::CaslLQRSystemDynamics() {
}

// Destructor
CaslLQRSystemDynamics::~CaslLQRSystemDynamics() {
}

// Helper function to compute time-dependent Riccati solution
void CaslLQRSystemDynamics::solveRiccati(double t, double& P11, double& P12, double& P22) const {
    double T = 5.0;
    double tau = T - t;  // Time-to-go

    // For LQR with A = [0 1; 0 0], B = [0; 1], Q = I, R = 1
    // The analytical solution is:
    P11 = tau + (1.0/3.0)*pow(tau,3);
    P12 = (1.0/2.0)*pow(tau,2);
    P22 = tau;
}

// Member functions
void CaslLQRSystemDynamics::LQRDynamics(double x1, double x2, double &fx1, double &fx2) {
    /*
     * xdot = [0 1] * x + [0]   *
     * xddot = [0 0] * xdot + [1] * u
     * */
    fx1 = x2;
    fx2 = 0;
}

// Exact solution computation with time dependence
double CaslLQRSystemDynamics::exactSolution(double x1, double x2, double t) const {
    double P11, P12, P22;
    solveRiccati(t, P11, P12, P22);

    // Value function V(x,t) = (1/2)xᵀP(t)x
    return 0.5 * (P11 * x1 * x1 + 2.0 * P12 * x1 * x2 + P22 * x2 * x2);
}

// Get optimal control
double CaslLQRSystemDynamics::getOptimalControl(double x1, double x2, double t) const {
    double P11, P12, P22;
    solveRiccati(t, P11, P12, P22);

    // u*(x,t) = -R⁻¹BᵀP(t)x = -(P12*x1 + P22*x2)
    return -(P12 * x1 + P22 * x2);
}

// Add error computation
double CaslLQRSystemDynamics::computeMaxError(const CaslGrid2D& grid,
                                              const CaslArray2D<double>& numericalSolution,
                                              double t) {
    double maxError = 0.0;

    // Loop over grid points
    for(int i = 1; i <= grid.nX(); i++) {
        for(int j = 1; j <= grid.nY(); j++) {
            double x1 = grid.x(i);
            double x2 = grid.y(j);

            // Compute exact solution at this point
            double exactVal = exactSolution(x1, x2, t);

            // Compute error
            double error = fabs(exactVal - numericalSolution(i,j));

            // Update max error
            maxError = std::max(maxError, error);
        }
    }

    return maxError;
}

// Add L2 error computation
double CaslLQRSystemDynamics::computeL2Error(const CaslGrid2D& grid,
                                             const CaslArray2D<double>& numericalSolution,
                                             double t) {
    double l2Error = 0.0;

    // Loop over grid points
    for(int i = 1; i <= grid.nX(); i++) {
        for(int j = 1; j <= grid.nY(); j++) {
            double x1 = grid.x(i);
            double x2 = grid.y(j);

            // Compute exact solution at this point
            double exactVal = exactSolution(x1, x2, t);

            // Compute squared error
            double error = pow(exactVal - numericalSolution(i,j), 2);

            // Add to L2 error
            l2Error += error * grid.dx() * grid.dy(); // Include grid cell area
        }
    }

    return sqrt(l2Error); // Return square root for L2 norm
}

// Verify HJB equation (optional)
double CaslLQRSystemDynamics::verifyHJB(double x1, double x2, double t) const {
    double P11, P12, P22;
    solveRiccati(t, P11, P12, P22);

    double T = 5.0;
    double tau = T - t;

    // Compute -∂V/∂t
    double dP11dt = -1.0 - tau*tau;
    double dP12dt = -tau;
    double dP22dt = -1.0;
    double Vt = -0.5 * (dP11dt*x1*x1 + 2.0*dP12dt*x1*x2 + dP22dt*x2*x2);

    // Compute optimal control
    double u = -(P12 * x1 + P22 * x2);

    // Compute Hamiltonian
    double H = 0.5*(x1*x1 + x2*x2) + x2*(P11*x1 + P12*x2) - 0.5*u*u;

    // Should be approximately zero if solution is correct
    return Vt + H;
}

#endif //CASL_LQR_SYSTEM_DYNAMICS_CPP