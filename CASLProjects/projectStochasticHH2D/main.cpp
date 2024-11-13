//
// Created by Faranak Rajabi on 02/26/2024.
//

/**
 * Description:
 * - Hamilton-Jacobi-Bellman (HJB) Equation with stochastic forcing term:
 *   ∂u/∂t + H(t, ∇u) + (D)∇²u = 0
 *
 * - Modified HJB Equation for our stochastic HH model:
 *   dv/dt + H(z, t, v, ∇v) + (D)d²v/dx² = 0
 *
 * Numerical Solution Procedure:
 * 1. Calculate ∇v using WENO method
 * 2. Compute H_hat using LLF method with quadratic extrapolation BC
 * 3. Discretize the diffusion term using backward time central differencing with quadratic extrapolation BC
 * 4. Integrate v(z, t) using TVD RK method
 *
 * Note: The PDE is solved backward in time starting with a terminal payoff at the final time, otherwise there will be an
 * unstable diffusion term.
 *
 * Numerical Hamiltonian (The same as the deterministic case):
 * Given by the formula:
 * (vⁿ⁺¹ - vⁿ) / Δt + H_hatⁿ(v_x_minus, v_x_plus, v(y,z))
 *
 * LLF solution to H_hat is the same as the deterministic case:
 * H( (v_x_minus + v_x_plus) / 2, (v_y_minus + v_y_plus) / 2 )
 *   - α_x( (v_x_plus + v_x_minus) / 2 )
 *   - α_y( (v_y_plus + v_y_minus) / 2 )
 * where:
 *  α_x = max|H₁|
 *  α_y = max|H₂|
 *
 * The only difference is the additional diffusion term discretized using BTCS:
 * (D/dx²)*(v(i+1) - 2*v(i) + v(i-1)) @ time step n + 1
 *
 * Definitions:
 *       z == (x, y) = ((1/K)v, n) implies ż = F(z) + Bu
 *       B = [(1/K), 0].T
 *       F(z) = [fx(z), fy(z)].T
 *       fx = (1/K)fv(Kx, y)
 *       fy = fn(Kx, y)
 *       H₁ = ∂H/∂vx = (-1/2K²) * vx + fx
 *       H₂ = ∂H/∂vy = fy
 */

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <array>

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"
#include "projectStochasticHH2D_lib/CaslHHNeuronModel.h"
#include "projectStochasticHH2D_lib/CaslHamiltonianHHModel.h"

using namespace std;
namespace fs = std::filesystem;

double heatSourceFunction                (CaslGrid2D & grid, int i, int j, double time, int dimension);
double calculateIC                       (double x, double y, double gama, double xTarget, double yTarget, double sigma     );
double minCalculator                     (const CaslGrid2D& grid, const CaslArray2D<double>& un                             );
double maxCalculator                     (const CaslGrid2D& grid, const CaslArray2D<double>& un                             );
void   exportToMatlab                    (const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName);
void   save2DTRUE2D                      (const CaslGrid2D& grid, const CaslArray2D<double>& phi                            );
void   findMinMax                        (CaslArray2D<double> f, double & minVal, double & maxVal, int & iMin, int & jMin   );
void   ensureDirectoryExists             (const std::string& directoryPath);
CaslArray2D<double> computeOptimalControl(const CaslArray2D<double> & phi_x, CaslGrid2D grid, double uMax, double K);
CaslArray2D<double> computePhi_x         (CaslArray2D<double> & phi_n, CaslGrid2D grid, CaslOptionPaddingWith BC, int derivativeDirection);

int main() {

    auto start_time = chrono::high_resolution_clock::now();

    // Project Structure:
    // Set Directories
    //    CASLHJB2D/
    //    └── CASLProjects/
    //      └── projectStochastic2D/
    //          └── __Output/
    //              ├── phi/
    //              ├── uStar/
    string projectPath = fs::current_path().string() + "/";  // Gets current working directory
    string outputFolder = projectPath + "__Output/";
    string phiFolder = outputFolder + "phi/";
    string uStarFolder = outputFolder + "uStar/";

    ensureDirectoryExists(outputFolder);
    ensureDirectoryExists(phiFolder);
    ensureDirectoryExists(uStarFolder);

    // Pre-defined Constants:
    const double K      = 100.0;  // scale factor
    const double uMax   = 10.0 ;  // maximum current
    const double vSpike = 44.8 ;  // mV 44.2
    const double nSpike = 0.459;  //
    const double vTh    = -20.0;  // Threshold mean voltage that triggers input
    const double tTh    = 0.0  ;  // Threshold time that activates input
    double       D      = 1.0;  // Mean of white noise(stochastic term coefficient)
    double diffusionC = D / (K * K);

    // Grid Set up:
    const double vMin = -100.0, vMax = 100.0;
    const double nMin = 0.0   , nMax = 1.0  ;
    auto xMin = (1.0 / K) * vMin;
    auto xMax = (1.0 / K) * vMax;
    auto yMin = nMin;
    auto yMax = nMax;

    const int nX = 320, nY = 320;
    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    const int nPadsX = 3; // 3 pads on each end since we will use at most ENO3 or WENO5.
    const int nPadsY = 3; // 3 pads on each end since we will use at most ENO3 or WENO5.

    CaslArray2D<double> phi_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_1   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // Time Set Up:
    double tInitial = 0, tFinal = 7.0, currentTime = tInitial;
    double dt = 0.0;
    double deltaT = 2.221035105269855e-04; // In mohammad's code
    double deltaT2 = 5.670321693340182e-05;

    // Periodic Orbit Set:
    std::vector<double> vSetInit, nSetInit;
    std::vector<double> InitC = {nSpike, vSpike};

    int DtFactor = 400;
    int tPlotEnd = DtFactor + 1;
    int numTimeSteps = 2 * tPlotEnd;
    double tvnData[numTimeSteps][3];

    std::vector<double> vDataVec, nDataVec;

    std::ofstream tvnFile("output_tvn_PO.txt");
    for (int i = 0; i < numTimeSteps; i++) {
        tvnFile << tvnData[i][0] << " "
                << tvnData[i][1] << " "
                << tvnData[i][2] << std::endl;
        vDataVec.push_back(tvnData[i][1]);
        nDataVec.push_back(tvnData[i][2]);
    }
    tvnFile.close();

    // Initial Function @ Phase-less Set = Target Point
    const double v_pl = -59.6;  // phase-less target point
    const double n_pl = 0.403;  // phase-less target point
    auto xTarget = (1.0 / K) * v_pl;
    auto yTarget = n_pl;

    const double gama  = 1000.0;
    const double sigma = 0.001 ;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            phi_n(i, j) = calculateIC(grid.x(i), grid.y(j), gama, xTarget, yTarget, sigma);
        }
    }

    // Defining F(z) on the whole grid, from HH-Neuorn Model class:
    CaslHHNeuronModel HH2DModel;
    CaslArray2D<double> fxVec(nX, nY), fyVec(nX, nY);
    for (int i = 1; i <= grid.nX(); i++)
        for (int j = 1; j <= grid.nY(); j++) {
            fxVec(i, j) = (1.0 / K) * HH2DModel.calculateVdotSN(K * grid.x(i), grid.y(j));
            fyVec(i, j) =             HH2DModel.calculateNdotSN(K * grid.x(i), grid.y(j));

        }

    // Chose which Hamiltonian and Laplacian to consider:
    CaslOptionPaddingWith               boundaryCondition        = withLinearExtrapolation;
    CaslOptionNumericalFirstDerivative  firstDerivativeScheme    = WENO5                  ;
    CaslOptionSecondOrderTermDirection  secondOrderTermDirection = X                         ;
    CaslOptionNumericalSecondDerivative secondDerivativeScheme   = BackwardTimeCentralSpacing;

    CaslHamiltonianHHModel              hamiltonian    (grid, fxVec, fyVec);
    CaslHamiltonJacobi2D       <double> HJSolver       (grid, hamiltonian, dt, currentTime, firstDerivativeScheme);
    CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, currentTime, diffusionC, secondOrderTermDirection, secondDerivativeScheme, boundaryCondition, HJSolver, heatSourceFunction);

    // Export data to MATLAB:
    CaslCppToMATLAB2D CppToMatlab;

    // Backward in time simulation:
    const double CFLNumber         = 0.9;
    int          time_step_number  = 0;
    int          output_number     = 0;
    double minVal, maxVal;
    int    iMin  , jMin;

    CaslArray2D<double> phi_x(nX, nY);
    CaslArray2D<double> uStar(nX, nY);
    // Compute the u_star based on phi_0 == u_star0
    // uStar = computePhi_x(phi_n, grid, boundaryCondition, 1, uMax);


    while (currentTime < tFinal) {
        dt = (CFLNumber * LaplacianSolver.findCFLHamiltonianDiffusion(phi_n, currentTime));
//        if (currentTime == tInitial) {
//            dt = deltaT;
//        }
//        else if (currentTime + dt <= tFinal) {
//            dt = deltaT2;
//        }
        if (currentTime + dt > tFinal) {
            currentTime = tFinal;
            dt = currentTime - tFinal;
        }

        // Computing the uStar
        phi_x = computePhi_x(phi_n, grid, boundaryCondition, 1);
        uStar = computeOptimalControl(phi_x, grid, uMax, K);

        if (time_step_number % 200 == 0 || currentTime == tFinal) {
            // Compute the minimum value point
            findMinMax(phi_n, minVal, maxVal, iMin, jMin);

            // Exporting the Cost-to-go data(phi) and print out on console
            cout << "dt: " << dt << endl;
            cout << "phi_" + std::to_string(output_number) + ", Export to Matlab at currentTime = " << currentTime
                 << endl;
            cout << "phi(" + std::to_string(iMin) + ", " + std::to_string(jMin) + "): " << minVal << endl;

            std::string phi_fileName = phiFolder + "phi_" + std::to_string(output_number) + ".dat";
            CppToMatlab.exportDataToMatlab(grid, phi_n, phi_fileName);

            // Exporting the optimal control(u_star) and print out on console
//            cout << "uStar_" + std::to_string(output_number) + ", Export to Matlab at currentTime = " << currentTime
//                 << endl;
            cout << "uStar(" + std::to_string(iMin) + ", " + std::to_string(jMin) + "): " << uStar(iMin, jMin) << endl;

            std::string phi2_fileName = uStarFolder + "uStar_" + std::to_string(output_number) + ".dat";
            // uStar needs modification
            CppToMatlab.exportDataToMatlab(grid, uStar, phi2_fileName);
            output_number++;
        }

        // TVD-RK3
        LaplacianSolver.backwardTimeCentralSpacing(phi_n, phi_np1);

        LaplacianSolver.backwardTimeCentralSpacing(phi_np1, phi_np2);
        phi_np12 = .75 * phi_n + .25 * phi_np2;

        LaplacianSolver.backwardTimeCentralSpacing(phi_np12, phi_np32);
        phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;

        // Update for next time step:
        currentTime += dt;
        phi_n = phi_np1;
        time_step_number++;
    }

    auto end_time = std::chrono::high_resolution_clock::now();

    // Calculate the duration between start and end
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    // Output the elapsed time in microseconds
    std::cout << "Time for whole simulation: " << duration.count() << " microseconds" << std::endl;

    return EXIT_SUCCESS;
}


double calculateIC(double x, double y, double gama, double xTarget, double yTarget, double sigma) {
    return gama * (1 - exp(-(pow((x - xTarget), 2) / sigma) -
                            (pow((y - yTarget), 2) / sigma)));
}

double minCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un) {
    int nX, nY;
    nX = grid.nX();
    nY = grid.nY();

    double minVal = 100000;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            if (un(i, j) <= minVal) {
                minVal = un(i, j);
            }
        }
    }
    return minVal;
}

double maxCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un) {
    int nX, nY;
    nX = grid.nX();
    nY = grid.nY();

    double maxVal = -100000;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            if (un(i, j) >= maxVal) {
                maxVal = un(i, j);
            }
        }
    }
    return maxVal;
}

void findMinMax(CaslArray2D<double> f, double & minVal, double & maxVal, int & iMin, int & jMin) {
    int nX = f.nX(), nY = f.nY();
    minVal =  std::numeric_limits<double>::infinity();
    maxVal = -std::numeric_limits<double>::infinity();

    for (int i = 1; i <= nX; i++) for (int j = 1; j <= nY; j++) {
            if (f(i, j) <= minVal) {
                minVal = f(i, j);
            }
            if (f(i, j) >= maxVal) {
                maxVal = f(i, j);
            }
        }

    bool foundMin = false;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            if (f(i, j) == minVal) {
                iMin = i;
                jMin = j;
                foundMin = true;
                break;
            }
        }
        if (foundMin) break;
    }

}

double heatSourceFunction(CaslGrid2D & grid, int i, int j, double time, int dimension) {
    return 0.0;
}

CaslArray2D<double> computePhi_x(CaslArray2D<double> & phi_n, CaslGrid2D grid, CaslOptionPaddingWith BC, int derivativeDirection) {
    int nX = grid.nX(), nY = grid.nY();
    double dx = grid.dx(), dy = grid.dy();

//    phi_n.extrapolate with linear

    // phi_x = phi(i+1) - phi(i-1) / 2dx
    CaslArray2D<double> phi_n_x(nX, nY);
    if (BC == withLinearExtrapolation) {
        phi_n.fillPaddingPoints(withLinearExtrapolation);
    }
    else if (BC == withQuadraticExtrapolation) {
        phi_n.fillPaddingPoints(withQuadraticExtrapolation);
    }

    if (derivativeDirection == 1) {
        for (int i = 1; i <= nX; i++) {
            for (int j = 1; j <= nY; j++) {
                    phi_n_x(i, j) = (phi_n(i + 1, j) - phi_n(i - 1, j)) / 2.0 / dx;
                // Done with x direction
            }
        }
    }
    return phi_n_x;
}

CaslArray2D<double> computeOptimalControl(const CaslArray2D<double> & phi_x, CaslGrid2D grid, double uMax, double K) {
    int nX = grid.nX(), nY = grid.nY();
    CaslArray2D<double> uStar(nX, nY);

    std::ofstream phi_x_ofstream("phi_x");
    std::ofstream uStar_ofstream("uSTar");

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {

            if (fabs(phi_x(i, j)) <= (2.0 * K * uMax)) {
                uStar(i, j) = (-0.5 * phi_x(i, j)) / K;
            }
            if (fabs(phi_x(i, j)) > (2.0 * K * uMax)) {
                if (phi_x(i, j) > 0) {
                    uStar(i, j) = -uMax;
                }
                if (phi_x(i, j) <= 0) {
                    uStar(i, j) = uMax;
                }
            }
            phi_x_ofstream << phi_x(i, j) << " ";
            uStar_ofstream << uStar(i, j) << " ";
        }
        uStar_ofstream << std::endl;
        phi_x_ofstream << std::endl;
    }
    return uStar;
}

void ensureDirectoryExists(const std::string& directoryPath) {
    try {
        if (fs::exists(directoryPath)) {
            std::cout << "Directory already exists: " << directoryPath << std::endl;
            return;
        }

        fs::create_directories(directoryPath);
        std::cout << "Directory created: " << directoryPath << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error with directory: " << directoryPath << ". " << e.what() << std::endl;
    }
}