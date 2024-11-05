//
// Created by Faranak Rajabi at 01/17/24
//

/**
 * Solving the Advection-Diffusion types equations
 *
 * Equations:
 *  u_t = D * u_xx
 *
 * */

// u(x, 0) = sin(pi * x)
// 0 <= x <= 1.0
// BC: Dirichlet
// u(x, t) = sin(pi * (x - ct)) * exp(-pi^2 * D * t)

// u_t + H = 0;
// ui_n+1 - ui_n / dt = -H --> ui_n+1 = ui_n - H * dt

// ui_n+1 - ui_n / dt + H_hat(Hamiltonian i, n) - D / dx^2 (ui+1 - 2ui + ui-1)_n+1 = F;
// [-alpha * ui+1 + (1+2alpha) ui - alpha * u_i-1]_n+1 = ui - dt*H_hat + dt*F
// Ax = b.

// alpha = D * dt / dx^2,

// Todo: Test the example 4(advection-diffusion) with periodic BC after the deadline
// Todo: Test the Fokker-Plank equation with linear extrapolation BC
// Todo: Test the linear SDE example by Ian-Mitchell
// Todo: Test the Fokker-Plank equation with a known solution p = exp(-t) * (cos(x) + sin(y))

// Test case 5: Linear SDE example by Ian Mitchel
// dphi/dt + a.dphi/dx + (1/2)trace(b^2D2_x*phi) = 0
// a = 1, b = 1, tf = 0.25, dim = 1
// phi = E[g(x)]
// where x(t) solves: dx = ax * dt + b dB;
// we assume that E[x(T)] = x, and Var(x) = E[x62] - E[x]^2.
// a : c1 (initial velocity), 0.5 * b^2 = D (Diffusion Coefficient)

#include <cmath>
#include <iostream>
#include <fstream>

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "projectLaplacian2D_1D_lib/CaslInitialProfiles2D.h"
#include "projectLaplacian2D_1D_lib/CaslHamiltonianDiffusion2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"

using namespace std;

double heatSourceFunction(CaslGrid2D & grid, int i, int j, double time, int dimension);
CaslArray2D<double> hamiltonianReal(CaslGrid2D & grid, CaslArray2D<double> & c1, double tFinal, double time);
void   maxErrorOnGrid    (CaslArray2D<double> u_approximated, CaslArray2D<double> u_exact, double & max_error, double & max_relative_error, int & i_max_error, int & j_max_error);

int main() {
    string testFolderName = "linearSDE"; // NEVER USE THIS: "heat1D_test", advection_diffusion_cos1D

    // Steps:
    // 1. Grid
    // 2. Padding
    // 3. Time set up
    // 4. Initial function definition
    // 5. CaslSecondOrder2D constructor
    //

    // Grid set up:
    int nX = 400, nY = 2;

    // For example 1 & 2:
     double xMin = -2.0, xMax = 2.0,
            yMin = 0.0, yMax = xMax;

    // For example 3:
    // double xMin = 0.0  , xMax = 2.0,
    //        yMin = 0.0  , yMax = 2.0;

    // For example 5: (linear SDE)
    //    % Create the grid.
    //    if(nargin < 5)
    //        g.dim = 1;
    //    else
    //        g.dim = dim;
    //    end
    //    g.min = -2;
    //    g.dx = 1 / 50;
    //    g.max = +2;
    //    %g.bdry = @addGhostExtrapolate;
    //    g.bdry = @addGhostExtrapolate2;
    //    g = processGrid(g);

//    double xMin = -2.0  , xMax = 2.0,
//           yMin = -2.0  , yMax = 2.0;

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    // Time set up:
    double initialT = 0.0, finalT = 0.25, current_time = initialT;
    double dt;
    // int num_time_steps = 2 * std::max(grid.nX(), grid.nY());
    // For accuracy test:
     dt = 0.5 * pow(std::min(grid.dx(), grid.dy()), 1);
    // dt = (finalT - initialT) / num_time_steps;

    // Phi arrays set up:
    const int nPadsX = 3;
    const int nPadsY = 3;

    CaslArray2D<double> phi_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // For example 5: (linear SDE) we have two functions: E[x] and Var[x]
    CaslArray2D<double> phi2_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi2_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi2_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi2_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi2_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // Initial function:
    // For example 1 & 2:
    // double L = 1;
    // For example 3:
    // double L = 2;
    // For example 1, 2, 3 & 5:
     double b_squared = 1.0;
     double diffusionC = 0.5 * b_squared;
//    double diffusionC = 1.0 / pow(M_PI, 2);
    // double diffusionC = 1;
    double L = 1.0;
    int dimension = 1.0; // X-direction

    // Inputs for the hamiltonian constructor, and initial functions if needed
    // For example 1, 2 & 3: Hamiltonian's all zeros, change in the following loop if needed
    // CaslArray2D<double> fx(nX, nY), fy(nX, nY);
    // For example 4: linear advection
    // For example 5: linear SDE a : c1
    CaslArray2D<double> c1(nX, nY); CaslArray2D<double> c2(nX, nY);
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            c1(i, j) = 1.0;
            c2(i, j) = 0.0;
        }
    }

    CaslInitialProfiles2D<double> initialProfile(grid, dt, current_time, b_squared, L, c1);
    // For example 1:
    // initialProfile.heatSinIP1D(phi_n, dimension);

    // For example 2:
    // initialProfile.heatSinePolyIP1D(phi_n, 1, current_time);

    // For example 3:
    // initialProfile.heatNonSmoothLineIP1D(phi_n, 1);

    // For example 4:
    // initialProfile.advectionDiffusionCosIP1D(phi_n, 1);

    // For example 5:
     initialProfile.linearSDE(phi_n, phi2_n, dimension, current_time);

    // initialProfile.heatCosIP1D(phi_n, 1);

    // initialProfile.advectionDiffusionCosIP1D(phi_n, 1);

    // initialProfile.advectionSquaredIP1D(phi_n, 1);

    // initialProfile.advectionSinIP1D(phi_n, dimension);

    // initialProfile.advectionDiffusionSinIP1D(phi_n, dimension);

    // initialProfile.advectionDiffusionExpIP1D(phi_n, dimension);

    // initialProfile.advectionDiffusionCosIP1D(phi_n, dimension);

    // initialProfile.linearSDEForward(phi_n, dimension);

    // Export phi_0 to MATLAB
    CaslCppToMATLAB2D cppToMatlab;
    cppToMatlab.exportDataToMatlab(grid, phi_n, "phi_0");

    // Inhomogeneous term to the heat equation
    // Note: double heatSourceFunction(CaslGrid2D & grid, int i, int j, double time, int dimension), pass it to the CaslSecondOrderDerivative2D instance;
    // For example 1: No heat source, all zeros.
    // For example 2: u(x, 0) = x.
    std::ofstream hHatRealOfstrm("H_hat_real.dat");
    CaslArray2D<double> hamiltonainRealArray = hamiltonianReal(grid, c1, finalT, current_time);

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            hHatRealOfstrm << hamiltonainRealArray(i, j) << " ";
        }
        hHatRealOfstrm << endl;
    }

    // BC: Dirichlet --> linear extrapolation
    CaslOptionPaddingWith boundaryCondition = withQuadraticExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
    CaslOptionSecondOrderTermDirection secondOrderTermDirection = X;
    CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

    CaslHamiltonianDiffusion2D hamiltonianForDiffusion(grid, c1, c2);
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianForDiffusion, dt, current_time, firstDerivativeScheme);
    CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, current_time, diffusionC, secondOrderTermDirection, secondDerivativeScheme, boundaryCondition, HJSolver, heatSourceFunction);

    // Parameters required for the simulation
    int i = 0;
    int time_step_number = 0;
    double CFL_number = 0.1;

    string folderName        = "../__Output/" + testFolderName + "/";
    string maxErrorsFileName = folderName     + "/maxErrors_"  + to_string(nX) + "_" + to_string(nY) + ".dat";
    ofstream maxErrors(maxErrorsFileName);
    maxErrors << "time" << " " << "phi_exact" << " " << "phi_approx" << " " << "Error" << std::endl;

    CaslArray2D<double> phi_n_exact  (phi_n.nX(), phi_n.nY());
    CaslArray2D<double> phi_n_exact_2(phi_n.nX(), phi_n.nY());

    // Backward in time simulation:
    while (current_time < finalT) {
        // Compute dt using findCFL function:
         // dt = (CFL_number * LaplacianSolver.findCFLHamiltonianDiffusion(phi_n, current_time));
         // dt = HJSolver.findCFL(phi_n, current_time);
         if (current_time + dt > finalT) {
             dt = finalT - current_time;
             current_time = finalT;
         }

        // if (stepCounter * exportThreshold > exportCounter)
        if (time_step_number % 5 == 0 || current_time == finalT) {
                cout << "phi_" + std::to_string(i) + ", Export to Matlab at currentTime = " << current_time << endl;
                std::string phi_fileName = folderName + "phi_" + std::to_string(i) + ".dat";
                cppToMatlab.exportDataToMatlab(grid, phi_n, phi_fileName);
                std::string phi2_fileName = folderName + "phi2_" + std::to_string(i) + ".dat";
                cppToMatlab.exportDataToMatlab(grid, phi2_n, phi2_fileName);

                // Check on the error:
                // Phi_exact for example 1:
                // initialProfile.heatSinIP1DExact(phi_n_exact, 1, current_time);
                // Phi_exact for example 2:
                // initialProfile.heatSinePolyIP1DNonHomogeneousExact(phi_n_exact, 1, current_time);
                // Phi_exact for example 3:
                // initialProfile.heatSinePolyIP1DNonHomogeneousExact(phi_n_exact, 1, current_time);
                // Phi_exact for example 4;
                // Phi_exact_1 and phi_exact_2 for example 5:

                // initialProfile.linearSDEForwardExact(phi_n_exact, dimension,current_time);
                // initialProfile.heatCosIP1DExact(phi_n_exact, dimension, current_time);
                // initialProfile.advectionDiffusionCosIP1DExact(phi_n_exact, 1, current_time);
                // initialProfile.advectionSinIP1DExact(phi_n_exact, 1, current_time);
                // initialProfile.advectionDiffusionSinIP1DExact(phi_n_exact, dimension, current_time);
                // initialProfile.advectionDiffusionExpIP1DExact(phi_n_exact, dimension, current_time);
                // initialProfile.advectionDiffusionCosIP1DExact(phi_n_exact, dimension, current_time);
                initialProfile.linearSDEExact(phi_n_exact, phi_n_exact_2, dimension, current_time);
                std::string phiExactFileName = folderName + "phi_exact" + to_string(i) + ".dat";
                cppToMatlab.exportDataToMatlab(grid, phi_n_exact, phiExactFileName);
                std::string phi2ExactFileName = folderName + "phi2_exact" + to_string(i) + ".dat";
                cppToMatlab.exportDataToMatlab(grid, phi_n_exact_2, phi2ExactFileName);
                double maxError1, maxError2, maxRelError;
                // Note: i_maxError, j_maxError computed at the location of max error on grid
                int i_maxError, j_maxError;
                maxErrorOnGrid(phi_n, phi_n_exact, maxError1, maxRelError, i_maxError, j_maxError);

                cout << "Max Error 1: " << maxError1 << endl;
                cout << "Max Relative Error: " << maxRelError << endl;

                maxErrors << current_time << " " << phi_n_exact(i_maxError, j_maxError) << " " << phi_n(i_maxError, j_maxError) <<
                          " " << maxError1 << endl;

                int i_maxError2, j_maxError2;
                double maxRel2;
                maxErrorOnGrid(phi2_n, phi_n_exact_2, maxError2, maxRel2, i_maxError2, j_maxError2);
                cout << "Max Error 2: " << maxError2 << endl;
                maxErrors << current_time << " " << phi_n_exact_2(i_maxError, j_maxError) << " " << phi2_n(i_maxError, j_maxError) <<
                          " " << maxError2 << endl;
                i++;
        }

        // TVD-RK3:
//        phi_n.fillPaddingPoints(boundaryCondition);
//        HJSolver.eulerStep(phi_n, phi_np1);

//        phi_np1.fillPaddingPoints(boundaryCondition);
//        HJSolver.eulerStep(phi_np1, phi_np2);
//
//        phi_np12 = .75 * phi_n + .25 * phi_np2;
//
//        phi_np12.fillPaddingPoints(boundaryCondition);
//        HJSolver.eulerStep(phi_np12, phi_np32);
//
//        phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;
//
        LaplacianSolver.backwardTimeCentralSpacing(phi_n, phi_np1);
        LaplacianSolver.backwardTimeCentralSpacing(phi2_n, phi2_np1);

        LaplacianSolver.backwardTimeCentralSpacing(phi_np1, phi_np2);
        LaplacianSolver.backwardTimeCentralSpacing(phi2_np1, phi2_np2);

        phi_np12 = .75 * phi_n + .25 * phi_np2;
        phi2_np12 = .75 * phi2_n + .25 * phi2_np2;

        LaplacianSolver.backwardTimeCentralSpacing(phi_np12, phi_np32);
        LaplacianSolver.backwardTimeCentralSpacing(phi2_np12, phi2_np32);

        phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;
        phi2_np1 = 1. / 3 * phi2_n + 2. / 3 * phi2_np32;

        // Update for next time step:
        current_time += dt;
        phi_n = phi_np1;
        phi2_n = phi2_np1;
        time_step_number++;
    }

    return EXIT_SUCCESS;
}

void maxErrorOnGrid(CaslArray2D<double> u_approximated, CaslArray2D<double> u_exact, double & max_error, double & max_relative_error, int & i_max_error, int & j_max_error) {
    double error_ij;
    double relative_error_ij;
    max_error = -INFINITY;
    max_relative_error = -INFINITY;
    for (int i = 1; i <= u_exact.nX(); i++) {
        for (int j = 1; j <= u_exact.nY(); j++) {
            error_ij = fabs(u_exact(i, j) - u_approximated(i, j));
            relative_error_ij = fabs((u_exact(i, j) - u_approximated(i, j)) / u_exact(i, j)) * 100.0;
            if (error_ij >= max_error) {
                max_error = error_ij;
                i_max_error = i;
                j_max_error = j;
            }
            if(relative_error_ij >= max_relative_error) {
                max_relative_error = relative_error_ij;
            }
        }
    }
}

double heatSourceFunction(CaslGrid2D & grid, int i, int j, double time, int dimension) {
    // For example 1, 3, 4 & 5: (heat1D_x), homogenous
     return 0.0;

    // For example 2: (heat1D_x_nonhomogeneous)
    // double x = (dimension == 1) ? grid.x(i) : grid.y(j);
    // return x;
}

CaslArray2D<double> hamiltonianReal(CaslGrid2D & grid, CaslArray2D<double> & c1, double tFinal, double time) {
    auto nX = grid.nX(), nY = grid.nY();
    CaslArray2D<double> hamiltonian_Real(nX, nY);
    double hamiltonian_t;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
             hamiltonian_t = exp(c1(i, j) * (tFinal - time));
             hamiltonian_Real(i, j) = hamiltonian_t;
        }
    }
    return hamiltonian_Real;
}