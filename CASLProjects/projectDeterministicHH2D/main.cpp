//
// Created by Faranak Rajabi on 7/30/23.
//

/**
 * Description:
 *   - Hamilton-Jacobi (HJ) Equation:
 *       dphi/dt + H(t, grad(phi)) = 0
 *   - Modified Hamilton-Jacobi (HH) Equation (specific to our case):
 *       dv/dt + H(z, t, v, grad(v)) = 0
 *
 * Numerical Solution Procedure:
 *   1. Calculate grad(v) using ENO method
 *   2. Compute H_hat using LLF method
 *   3. Integrate v(z, t) to find v(z, t) using TVD RK method
 *
 * Numerical Hamiltonian:
 *   Given by the formula:
 *       (phi^{n+1} - phi^{n}) / delta_t + H_hat^{n}(phi_x_minus, phi_x_plus, phi_(y, z))
 *   LLF solution to H_hat is given by:
 *       H( (phi_x_minus + phi_x_plus) / 2, (phi_y_minus + phi_y_plus) / 2 )
 *       - alpha_x( (phi_x_plus + phi_x_minus) / 2 )
 *       - alpha_y( (phi_y_plus + phi_y_minus) / 2 )
 *   where:
 *       alpha_x = max |H1| over Ix, Iy
 *       alpha_y = max |H2| over Ix, Iy
 *       H1 represents the partial derivative dH/dphi_x
 *       H2 represents the partial derivative dH/dphi_y
 *
 * Numerical Hamiltonian (specific to our case):
 *   For |vx| <= 2Ku_max,
 *       H = grad(v).T * F(z) - (1/4K^2) * vx^2
 *   For |vx| > 2Ku_max,
 *       H = grad(v).T * F(z) + u_max^2 - |vx|u_max/K
 *   Definitions:
 *       z == (x, y) = ((1/K)v, n) implies zdot = F(z) + Bu
 *       B = [(1/K), 0].T
 *       F(z) = [fx(z), fy(z)].T
 *       fx = (1/K)fv(Kx, y)
 *       fy = fn(Kx, y)
 *       H1 = dH/dvx = (-1/2K^2) * vx + fx
 *       H2 = dH/dvy = fy
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <array>
#include <filesystem>
#include <chrono>

#include "CaslGrid2D.h"
#include "CaslArray2D.h"
#include "CaslHamiltonJacobi2D.h"
#include "projectDeterministicHH2D_lib/CaslHHWholeSystemModel.h"
#include "projectDeterministicHH2D_lib/CaslHamiltonianHHModel.h"

using namespace std;
namespace fs = std::filesystem;

double  calculateIC          (double x, double y, double gama, double xTarg, double yTarg, double sigma);
double  minCalculator        (const CaslGrid2D& grid, const CaslArray2D<double>& un);
double  maxCalculator        (const CaslGrid2D& grid, const CaslArray2D<double>& un);
void    save2DTRUE2D         (const CaslGrid2D& grid, const CaslArray2D<double> & phi);
void    exportToMatlab       (const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName );
void    ensureDirectoryExists(const std::string& directoryPath);

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
    const int scaleFac = Ks;        // scale factor
    const double uMax = Ib;         // maximum current
    const double vSpike = 44.8;     // mV 44.2
    const double nSpike = 0.459;    // 0.465

    // Grid Set up:
    const double vMin = -100, vMax = 100;
    const double xMin = (1.0 / scaleFac) * vMin, xMax = (1.0 / scaleFac) * vMax;
    const double yMin = 0, yMax = 1;
    const int nX = 321, nY = 321;

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    const int nPadsX = 3; // 3 pads on each end since we will use at most ENO3 or WENO5.
    const int nPadsY = 3; // 3 pads on each end since we will use at most ENO3 or WENO5.

    CaslArray2D<double> phi_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_1   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // Steps:
    /*
     * 1. using vs, ns for periodic orbit set
     * 2. passing that to penalty as the initial condition
     * 3. Backward integration TVD RK4, Tend:-dt:0
     */

    // Time Set Up:
    double tInitial = 0, tFinal = 7, currentTime = tInitial;

    // Periodic Orbit Set:
    std::vector<double> vSetInit, nSetInit;
    double vInit = vSpike;
    double nInit = nSpike;
    std::vector<double> InitC = {vInit, nInit};

    // Periodic orbit simulation parameters:
    int DtFactor = 400;
    int tPlotEnd = DtFactor + 1;
    int numTimeSteps = 2 * tPlotEnd;
    double tvnData[numTimeSteps][3];

    // Using Ts(Spiking time) to get the whole periodic orbit:
//    CaslRK42DSolver solver;
//    solver.rk4(hhs, InitC, tInitial, Ts, numTimeSteps, tvnData);

    std::vector<double> vDataVec, nDataVec;

    // Export to Matlab for Plotting:
    std::ofstream tvnFile("output_tvn.txt");
    for (int i = 0; i < numTimeSteps; i++) {
        tvnFile << tvnData[i][0] << " "
                << tvnData[i][1] << " "
                << tvnData[i][2] << std::endl;
        vDataVec.push_back(tvnData[i][1]);
        nDataVec.push_back(tvnData[i][2]);
    }
    tvnFile.close();

    // Initial Function:
    // Phase-less Set = Target Point
    const double v_pl = -59.6;  // phase-less target point
    const double n_pl = 0.403;  // phase-less target point

    // Global Parameters in FinalState Class:
    double gama = 1000;
    double xTarg = (1.0 / scaleFac) * v_pl;
    // double xTarg = (1.0 / scaleFac) * vSpike;
    double yTarg = n_pl;
    // double yTarg = nSpike;
    double sigma = 0.001;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            phi_n(i, j) = calculateIC(grid.x(i), grid.y(j), gama, xTarg, yTarg, sigma);
        }
    }

    double minPhi_0 = minCalculator(grid, phi_n);
    double maxPhi_0 = maxCalculator(grid, phi_n);
    cout << "min: " << minPhi_0 << ", max: " << maxPhi_0 << endl;

    int iMin, jMin;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            if (phi_n(i, j) == minPhi_0) {
                cout << "(i, j)_min: " << i << ", " << j << endl;
                iMin = i;
                jMin = j;
            }
        }
    }

    std::ofstream penaltyFile("phi0_PenaltyFOnGrid_data.txt");
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            penaltyFile << phi_n(i, j);
            if (j != nY) {
                penaltyFile << "\t";  // separate the data with tabs for .txt format
            }
        }
        penaltyFile << "\n";  // start a new line for each i value
    }
    penaltyFile.close();

    // Defining F(z) on the whole grid:
    CaslArray2D<double> fxVec(grid.nX(), grid.nY()), fyVec(grid.nX(), grid.nY());
    for (int i = 1; i <= grid.nX(); i++)
        for (int j = 1; j <= grid.nY(); j++) {
            fxVec(i, j) = (1.0 / scaleFac) * fv(scaleFac * grid.x(i), grid.y(j));
            fyVec(i, j) = fn(scaleFac * grid.x(i), grid.y(j));
        }

    // Chose which Hamiltonian to consider:
    CaslHamiltonianHHModel hamiltonian(grid, fxVec, fyVec);
    CaslOptionPaddingWith boundaryConditions = withLinearExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

    // HJSolver Initializing Parameters:
    double CFLNum = 0.9;
    double dt; // double dt = 0.000549828;
    CaslArray2D<double> uStar(nX, nY);
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, currentTime, firstDerivativeScheme);

    // Define exporting to Matlab parameters:
    int exportCounter = 0;              // To count how many time steps exported to MATLAB
    int i = 0;
    int stepCounter = 1;                // To count total steps
    std::vector<CaslArray2D<double>> uStarData;

    vector<double> phiMinData;
    phiMinData.push_back(phi_n(iMin, jMin));

    // Backward in time simulation:
    while (currentTime < tFinal) {
        // Compute dt using findCFL function:
        dt = (CFLNum * HJSolver.findCFL(phi_n, currentTime));

        // Finding the optimal control:
        HJSolver.computeOptimalControl(phi_n, scaleFac, uMax, uStar);

        // if (stepCounter * exportThreshold > exportCounter)
        if (currentTime + dt <= tFinal) {
            if (exportCounter % 120 == 0) {
                cout << "Export to Matlab at currentTime = " << currentTime << endl;
                std::string phi_fileName = phiFolder + "phi_" + std::to_string(i) + ".dat";
                std::string uStar_fileName = uStarFolder + "uStar_" + std::to_string(i) + ".dat";

                exportToMatlab(grid, phi_n, phi_fileName);
                exportToMatlab(grid, uStar, uStar_fileName);

                // Saving uStar data for forward integration
                uStarData.push_back(uStar);
                cout << "dt: " << dt<< endl;
                cout << "phi_" + std::to_string(i) << ":" << phi_n(iMin, jMin) << std::endl;
                phiMinData.push_back(phi_n(iMin, jMin));

                cout << "uStar: " << uStar(iMin, jMin) << std::endl;
                i++;
            }
        }

        if (currentTime + dt > tFinal) {
            dt = tFinal - currentTime;
            currentTime = tFinal;
            cout << "Export to Matlab at currentTime = " << currentTime << endl;
            std::string phi_fileName = phiFolder + "phi_" + std::to_string(i) + ".dat";
            std::string uStar_fileName = uStarFolder + "uStar_" + std::to_string(i) + ".dat";

            exportToMatlab(grid, phi_n, phi_fileName);
            exportToMatlab(grid, uStar, uStar_fileName);

            // Saving uStar data for forward integration
            uStarData.push_back(uStar);
            cout << "dt: " << dt<< endl;
            cout << "phi_" + std::to_string(i) << ":" << phi_n(iMin, jMin) << std::endl;
            phiMinData.push_back(phi_n(iMin, jMin));

            cout << "uStar: " << uStar(iMin, jMin) << std::endl;
            i++;
        }

        // TVD-RK3:
        phi_n.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_n, phi_np1);

        phi_np1.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_np1, phi_np2);

        phi_np12 = .75 * phi_n + .25 * phi_np2;

        phi_np12.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_np12, phi_np32);

        phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;

        // Update for next time step:
        currentTime += dt;
        phi_n = phi_np1;
        stepCounter++;
        exportCounter++;
    }

    // ofstream phiOnMinFile("HH2DphiOfMin_Upwind.dat");
    // ofstream phiOnMinFile("HH2DphiOfMin_WENO5.dat");
    // ofstream phiOnMinFile("HH2DphiOfMin_ENO3.dat");
    // for (int inx = 1; inx <= size(phiMinData); inx++) {
    //    phiOnMinFile << phiMinData[inx] << "\n";
    // }
    // phiOnMinFile.close();
    // save2DTRUE2D(grid, phi_n);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto sim_tim = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

    std::cout << "Time for the whole simulation: "  << setprecision(4) << sim_tim.count() << " minutes" << std::endl;
}

void save2DTRUE2D(const CaslGrid2D & grid, const CaslArray2D<double> & phi) {
    ofstream out;
    out.open("2DTRUE2D");
    if (out.fail()) {
        cout << "Error opening 2DTRUE2D. ABORT." << endl; exit(1);
    }
    for (int i = 1; i <= grid.nX(); i++) {
        for (int j = 1; j <= grid.nY(); j++) {
            out << phi(i,j) << " ";
        }
        out << endl;
    }
}

double calculateIC(double x, double y, double gama, double xTarg, double yTarg, double sigma) {
    return gama * (1 - exp(-(pow((x - xTarg), 2) / sigma) -
                           (pow((y - yTarg), 2) / sigma)));
}

void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName) {
    std::ofstream ofStream; ofStream.open(fileName);
    if (ofStream.is_open()) {
        for (int i = 1; i <= un.nX(); ++i) {
            for (int j = 1; j <= un.nY(); ++j) ofStream << std::setprecision(16) << un(i,j) << " ";
            ofStream << std::endl;
        }
        ofStream.close(); return;
    }
    // If file could not be open:
    std::cout << "File " << fileName << " could not be opened. EXITING." << std::endl; exit(1);
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
