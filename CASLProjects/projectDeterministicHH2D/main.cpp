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
 *   1. Calculate grad(v) using ENO/WENO method
 *   2. Compute H_hat using LLF method
 *   3. Integrate in time using TVD RK3
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <array>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <limits>
#include <vector>
#include <string>
#include <cstdlib>

#include "CaslGrid2D.h"
#include "CaslArray2D.h"
#include "CaslHamiltonJacobi2D.h"
#include "projectDeterministicHH2D_lib/CaslHHWholeSystemModel.h"
#include "projectDeterministicHH2D_lib/CaslHamiltonianHHModel.h"

using namespace std;
namespace fs = std::filesystem;

double calculateIC(double x, double y, double gama, double xTarg, double yTarg, double sigma);
double minCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un);
double maxCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un);
void save2DTRUE2D(const CaslGrid2D& grid, const CaslArray2D<double>& phi);
void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName);
void ensureDirectoryExists(const std::string& directoryPath);

int main() {
    auto start_time = chrono::high_resolution_clock::now();

    string projectPath  = fs::current_path().string() + "/";
    string outputFolder = projectPath + "__Output/";
    string phiFolder    = outputFolder + "phi/";
    string uStarFolder  = outputFolder + "uStar/";

    ensureDirectoryExists(outputFolder);
    ensureDirectoryExists(phiFolder);
    ensureDirectoryExists(uStarFolder);

    // Pre-defined Constants:
    const int scaleFac   = Ks;
    const double uMax    = Ib;
    const double vSpike  = 44.8;
    const double nSpike  = 0.459;

    // Grid setup:
    const double vMin = -100.0, vMax = 100.0;
    const double xMin = (1.0 / scaleFac) * vMin;
    const double xMax = (1.0 / scaleFac) * vMax;
    const double yMin = 0.0, yMax = 1.0;
    const int nX = 80, nY = 80;

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    const int nPadsX = 3;
    const int nPadsY = 3;

    CaslArray2D<double> phi_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_1   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // Time setup:
    double tInitial = 0.0, tFinal = 7.0, currentTime = tInitial;

    // Periodic orbit set-up
    std::vector<double> vSetInit, nSetInit;
    double vInit = vSpike;
    double nInit = nSpike;
    std::vector<double> InitC = {vInit, nInit};

    // Periodic orbit storage
    int DtFactor = 400;
    int tPlotEnd = DtFactor + 1;
    int numTimeSteps = 2 * tPlotEnd;

    // IMPORTANT:
    // This array is zero-initialized so it is not garbage if solver call stays disabled.
    double tvnData[2 * 401][3] = {{0.0}};
    std::vector<double> vDataVec, nDataVec;

    // If you want real orbit data, uncomment and provide Ts + solver:
    // CaslRK42DSolver solver;
    // solver.rk4(hhs, InitC, tInitial, Ts, numTimeSteps, tvnData);

    // Only export periodic orbit data if it was actually filled.
    // Right now this just writes zeros unless you enable the solver above.
    {
        std::ofstream tvnFile("output_tvn.txt");
        for (int k = 0; k < numTimeSteps; k++) {
            tvnFile << tvnData[k][0] << " "
                    << tvnData[k][1] << " "
                    << tvnData[k][2] << std::endl;
            vDataVec.push_back(tvnData[k][1]);
            nDataVec.push_back(tvnData[k][2]);
        }
    }

    // Initial function: phase-less target point
    const double v_pl = -59.6;
    const double n_pl = 0.403;

    double gama  = 1000.0;
    double xTarg = (1.0 / scaleFac) * v_pl;
    double yTarg = n_pl;
    double sigma = 0.001;

    if (sigma <= 0.0) {
        std::cerr << "ERROR: sigma must be positive.\n";
        return 1;
    }

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            phi_n(i, j) = calculateIC(grid.x(i), grid.y(j), gama, xTarg, yTarg, sigma);
            if (!std::isfinite(phi_n(i, j))) {
                std::cerr << "Nonfinite initial condition at (i,j)=(" << i << "," << j << ")"
                          << " x=" << grid.x(i) << " y=" << grid.y(j)
                          << " phi=" << phi_n(i, j) << "\n";
                return 1;
            }
        }
    }

    double minPhi_0 = minCalculator(grid, phi_n);
    double maxPhi_0 = maxCalculator(grid, phi_n);
    cout << "min: " << minPhi_0 << ", max: " << maxPhi_0 << endl;

    int iMin = 1, jMin = 1;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            if (phi_n(i, j) == minPhi_0) {
                cout << "(i, j)_min: " << i << ", " << j << endl;
                iMin = i;
                jMin = j;
            }
        }
    }

    {
        std::ofstream penaltyFile("phi0_PenaltyFOnGrid_data.txt");
        for (int i = 1; i <= nX; i++) {
            for (int j = 1; j <= nY; j++) {
                penaltyFile << phi_n(i, j);
                if (j != nY) penaltyFile << "\t";
            }
            penaltyFile << "\n";
        }
    }

    // Define F(z) on whole grid
    CaslArray2D<double> fxVec(grid.nX(), grid.nY()), fyVec(grid.nX(), grid.nY());

    for (int i = 1; i <= grid.nX(); i++) {
        for (int j = 1; j <= grid.nY(); j++) {
            const double v = scaleFac * grid.x(i);
            const double n = grid.y(j);

            fxVec(i, j) = (1.0 / scaleFac) * fv(v, n);
            fyVec(i, j) = fn(v, n);

            if (!std::isfinite(fxVec(i, j)) || !std::isfinite(fyVec(i, j))) {
                std::cerr << "Bad dynamics at (i,j)=(" << i << "," << j << ")\n"
                          << "  v=" << v << " n=" << n << "\n"
                          << "  fx=" << fxVec(i, j) << " fy=" << fyVec(i, j) << "\n";
                return 1;
            }
        }
    }

    CaslHamiltonianHHModel hamiltonian(grid, fxVec, fyVec);
    CaslOptionPaddingWith boundaryConditions = withLinearExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

    double CFLNum = 0.9;
    double dt = 0.0;
    CaslArray2D<double> uStar(nX, nY);
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, currentTime, firstDerivativeScheme);

    int exportCounter = 0;
    int iExport = 0;
    int stepCounter = 1;
    std::vector<CaslArray2D<double>> uStarData;

    vector<double> phiMinData;
    phiMinData.push_back(phi_n(iMin, jMin));

    while (currentTime < tFinal) {
        const double cflVal = HJSolver.findCFL(phi_n, currentTime);
        if (!std::isfinite(cflVal) || cflVal <= 0.0) {
            std::cerr << "Invalid CFL value at step " << stepCounter
                      << " currentTime=" << currentTime
                      << " cflVal=" << cflVal << "\n";
            return 1;
        }

        dt = CFLNum * (1.0 / cflVal);

        if (!std::isfinite(dt) || dt <= 0.0) {
            std::cerr << "Invalid dt at step " << stepCounter
                      << " currentTime=" << currentTime
                      << " dt=" << dt << "\n";
            return 1;
        }

        HJSolver.computeOptimalControl(phi_n, scaleFac, uMax, uStar);

        if (currentTime + dt <= tFinal) {
            if (exportCounter % 120 == 0) {
                cout << "Export to Matlab at currentTime = " << currentTime << endl;
                std::string phi_fileName   = phiFolder   + "phi_"   + std::to_string(iExport) + ".dat";
                std::string uStar_fileName = uStarFolder + "uStar_" + std::to_string(iExport) + ".dat";

                exportToMatlab(grid, phi_n, phi_fileName);
                exportToMatlab(grid, uStar, uStar_fileName);

                // uStarData.push_back(uStar);
                // cout << "dt: " << dt << endl;
                // cout << "phi_" + std::to_string(iExport) << ":" << phi_n(iMin, jMin) << std::endl;
                // phiMinData.push_back(phi_n(iMin, jMin));
                // cout << "uStar: " << uStar(iMin, jMin) << std::endl;
                iExport++;
            }
        }

        if (currentTime + dt > tFinal) {
            dt = tFinal - currentTime;
            currentTime = tFinal;

            cout << "Export to Matlab at currentTime = " << currentTime << endl;
            std::string phi_fileName   = phiFolder   + "phi_"   + std::to_string(iExport) + ".dat";
            std::string uStar_fileName = uStarFolder + "uStar_" + std::to_string(iExport) + ".dat";

            exportToMatlab(grid, phi_n, phi_fileName);
            exportToMatlab(grid, uStar, uStar_fileName);

            uStarData.push_back(uStar);
            cout << "dt: " << dt << endl;
            cout << "phi_" + std::to_string(iExport) << ":" << phi_n(iMin, jMin) << std::endl;
            phiMinData.push_back(phi_n(iMin, jMin));
            cout << "uStar: " << uStar(iMin, jMin) << std::endl;
            iExport++;
        }

        // TVD-RK3
        phi_n.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_n, phi_np1);

        phi_np1.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_np1, phi_np2);

        phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;

        phi_np12.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phi_np12, phi_np32);

        phi_np1 = (1.0 / 3.0) * phi_n + (2.0 / 3.0) * phi_np32;

        // Sanity check after update
        for (int ii = 1; ii <= nX; ++ii) {
            for (int jj = 1; jj <= nY; ++jj) {
                if (!std::isfinite(phi_np1(ii, jj))) {
                    std::cerr << "Nonfinite phi after RK update at step " << stepCounter
                              << " currentTime=" << currentTime
                              << " (i,j)=(" << ii << "," << jj << ")"
                              << " x=" << grid.x(ii)
                              << " y=" << grid.y(jj)
                              << " phi=" << phi_np1(ii, jj) << "\n";
                    return 1;
                }
            }
        }

        currentTime += dt;
        phi_n = phi_np1;
        stepCounter++;
        exportCounter++;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto sim_tim  = std::chrono::duration_cast<std::chrono::minutes>(end_time - start_time);

    std::cout << "Time for the whole simulation: "
              << std::setprecision(4) << sim_tim.count()
              << " minutes" << std::endl;

    return 0;
}

void save2DTRUE2D(const CaslGrid2D& grid, const CaslArray2D<double>& phi) {
    ofstream out("2DTRUE2D");
    if (out.fail()) {
        cout << "Error opening 2DTRUE2D. ABORT." << endl;
        exit(1);
    }
    for (int i = 1; i <= grid.nX(); i++) {
        for (int j = 1; j <= grid.nY(); j++) {
            out << phi(i, j) << " ";
        }
        out << endl;
    }
}

double calculateIC(double x, double y, double gama, double xTarg, double yTarg, double sigma) {
    const double r2 = (pow(x - xTarg, 2) + pow(y - yTarg, 2)) / sigma;
    return gama * (1.0 - exp(-r2));
}

void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName) {
    std::ofstream ofStream(fileName);
    if (ofStream.is_open()) {
        for (int i = 1; i <= un.nX(); ++i) {
            for (int j = 1; j <= un.nY(); ++j) {
                ofStream << std::setprecision(16) << un(i, j) << " ";
            }
            ofStream << std::endl;
        }
        ofStream.close();
        return;
    }

    std::cout << "File " << fileName << " could not be opened. EXITING." << std::endl;
    exit(1);
}

double minCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un) {
    double minVal = std::numeric_limits<double>::infinity();
    for (int i = 1; i <= grid.nX(); ++i) {
        for (int j = 1; j <= grid.nY(); ++j) {
            minVal = std::min(minVal, un(i, j));
        }
    }
    return minVal;
}

double maxCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un) {
    double maxVal = -std::numeric_limits<double>::infinity();
    for (int i = 1; i <= grid.nX(); ++i) {
        for (int j = 1; j <= grid.nY(); ++j) {
            maxVal = std::max(maxVal, un(i, j));
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