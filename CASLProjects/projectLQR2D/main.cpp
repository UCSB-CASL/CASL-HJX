//
// Created by Faranak Rajabi on 12/24/23.
//

#define CASL_LQR2D
#ifdef  CASL_LQR2D

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <filesystem>
#include <chrono>

// Modified paths to use CASLCommonLibrary instead of CASLUniformLibrary
#include "../CASLCommonLibrary/CaslArray2D.h"
#include "../CASLCommonLibrary/CaslGrid2D.h"
#include "../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "projectLQR2D_lib/CaslLQRSystemDynamics.h"
#include "projectLQR2D_lib/CaslHamiltonianLQR2D.h"
#include "../CASLCommonLibrary/CaslCppToMATLAB2D.h"

using namespace std;
namespace fs = std::filesystem;

double costFunction(double x, double xDot);
void findMinMax2D(const CaslGrid2D& grid, const CaslArray2D<double>& phi, double& minVal, double& maxVal,
                  int& iMin, int& jMin);
void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName);
void ensureDirectoryExists(const std::string& directoryPath);
void exportTimeStep(double currentTime, double dt, int& i,
                    const std::string& phiFolder,
                    const std::string& testCase,
                    const CaslGrid2D& grid, const CaslArray2D<double>& phi_n,
                    CaslCppToMATLAB2D& toMATLAB,
                    int iMin, int jMin,
                    double maxError = 0.0, double l2Error = 0.0, double hjbError = 0.0);
double computeL2Error(const CaslArray2D<double>& numerical, const CaslGrid2D& grid,
                      const CaslLQRSystemDynamics& lqr, double currentTime);
// Function to compute HJB residual at a point (optional helper function)
double computeHJBResidual(const CaslGrid2D& grid, const CaslArray2D<double>& phi_n,
                          const CaslLQRSystemDynamics& lqr, double currentTime,
                          int i, int j);

struct SimulationResult {
    int gridSize;
    double dx;
    double maxError;
    double l2Error;
    double simulationTime;  // in seconds
};

class ConvergenceAnalysis {
public:
    void addResult(int gridSize, double dx, double maxError, double l2Error, double hjbError, double simTime) {
        results.push_back({gridSize, dx, maxError, l2Error, hjbError, simTime});
    }

    void computeConvergenceRates() {
        std::cout << "\n=============== Convergence Analysis for LQR (T=5) ===============" << std::endl;
        std::cout << std::setw(10) << "Grid"
                  << std::setw(15) << "dx"
                  << std::setw(15) << "dt"
                  << std::setw(15) << "Max Error"
                  << std::setw(15) << "L2 Error"
                  << std::setw(15) << "HJB Error"
                  << std::setw(20) << "Conv. Rate"
                  << std::setw(15) << "Time (s)" << std::endl;
        std::cout << std::string(105, '-') << std::endl;

        // Print first result
        double dt = 0.5 * std::pow(results[0].dx, 2.0);
        std::cout << std::setw(10) << results[0].gridSize
                  << std::scientific << std::setprecision(4)
                  << std::setw(15) << results[0].dx
                  << std::setw(15) << dt
                  << std::setw(15) << results[0].maxError
                  << std::setw(15) << results[0].l2Error
                  << std::setw(15) << results[0].hjbError
                  << std::setw(20) << "-"
                  << std::fixed << std::setw(15) << results[0].simulationTime << std::endl;

        // Compute and print rates for subsequent results
        for (size_t i = 1; i < results.size(); ++i) {
            dt = std::pow(results[i].dx, 2.0);  // Using dx² for time step

            // Compute convergence rates for different error measures
            double maxErrorRate = computeConvergenceRate(
                    results[i-1].dx, results[i].dx,
                    results[i-1].maxError, results[i].maxError
            );

            double l2ErrorRate = computeConvergenceRate(
                    results[i-1].dx, results[i].dx,
                    results[i-1].l2Error, results[i].l2Error
            );

            double hjbErrorRate = computeConvergenceRate(
                    results[i-1].dx, results[i].dx,
                    results[i-1].hjbError, results[i].hjbError
            );

            std::cout << std::setw(10) << results[i].gridSize
                      << std::scientific << std::setprecision(4)
                      << std::setw(15) << results[i].dx
                      << std::setw(15) << dt
                      << std::setw(15) << results[i].maxError
                      << std::setw(15) << results[i].l2Error
                      << std::setw(15) << results[i].hjbError
                      << std::setw(20) << maxErrorRate
                      << std::fixed << std::setw(15) << results[i].simulationTime << std::endl;

            // Print detailed convergence info
            std::cout << "   Convergence Rates:" << std::endl;
            std::cout << "   Max Error: " << maxErrorRate << std::endl;
            std::cout << "   L2 Error:  " << l2ErrorRate << std::endl;
            std::cout << "   HJB Error: " << hjbErrorRate << std::endl;
        }
        std::cout << std::string(105, '-') << std::endl;
    }

private:
    struct SimulationResult {
        int gridSize;
        double dx;
        double maxError;
        double l2Error;
        double hjbError;
        double simulationTime;
    };

    std::vector<SimulationResult> results;

    double computeConvergenceRate(double dx1, double dx2, double error1, double error2) {
        if (error1 <= 0 || error2 <= 0) {
            return 0.0;  // Handle zero or negative errors
        }
        return std::log(error1/error2) / std::log(dx1/dx2);
    }
};

int main() {
    ConvergenceAnalysis convergence;
    std::vector<int> gridSizes = {20, 40, 80};

    for (int gridSize : gridSizes) {
        // Start timing
        auto startTime = std::chrono::high_resolution_clock::now();
        // 4D Grid Set-up:
        /* x1: x    (position in x-direction)
         * x2: xdot (Velocity in x-direction)
         */
        const double x1Min = -1, x1Max = 1;
        const double x2Min = -1, x2Max = 1;
        const int nX = gridSize, nY = gridSize;
        CaslGrid2D grid(x1Min, x1Max, x2Min, x2Max, nX, nY);
        const int nPadsX = 3, nPadsY = 3;
        cout << "dx: " << grid.dx() << endl;

        CaslArray2D<double> phi_n   (grid.nX(), grid.nY(), nPadsX, nPadsY);
        CaslArray2D<double> phi_np1 (grid.nX(), grid.nY(), nPadsX, nPadsY);
        CaslArray2D<double> phi_np2 (grid.nX(), grid.nY(), nPadsX, nPadsY);
        CaslArray2D<double> phi_np12(grid.nX(), grid.nY(), nPadsX, nPadsY);
        CaslArray2D<double> phi_np32(grid.nX(), grid.nY(), nPadsX, nPadsY);

        // Set Directories
        //    CASLHJB2D/
        //    └── CASLProjects/
        //      └── projectLQR2D/
        //          └── __Output/
        //              ├── LQR2D_[gridsize]/
        //                  └── phi/
        string projectPath = fs::current_path().string() + "/";  // Gets current working directory
        string outputFolder = projectPath + "__Output/";
        string testCase = outputFolder + "LQR2D_" + to_string(grid.nX()) + "/";
        string phiFolder = testCase + "phi/";

        // Create directories
        ensureDirectoryExists(outputFolder);
        ensureDirectoryExists(testCase);
        ensureDirectoryExists(phiFolder);

        // Time Set Up:
        double tInitial = 0, tFinal = 5, currentTime = tInitial;
        double dt = std::pow(grid.dx(), 2.0);  // Match what's in convergence analysis

        // Final Cost(Initial Function Set-up):
        for (int i = 1; i <= grid.nX(); i++) {
            for (int j = 1; j <= grid.nY(); j++) {
                double xi = grid.x(i);
                double yj = grid.y(j);
                double cost_ij = costFunction(xi, yj);
                phi_n(i, j) = cost_ij;
            }
        }

        //        std::ofstream finalCost("phi0(EndCost)_LQR2D.txt");
        //        for (int i = 1; i <= grid.nX(); i++) {
        //            for (int j = 1; j <= grid.nY(); j++) {
        //                finalCost << i << " " << j << " " << phi_n(i, j) << endl;
        //            }
        //        }

        // finalCost.close();

        double minVal, maxVal;
        int iMin, jMin;
        findMinMax2D(grid, phi_n, minVal, maxVal, iMin, jMin);
        cout << "minVal: " << minVal << endl;
        cout << "(i, j)_min: (" << iMin << ", " << jMin << ")" << endl;

        CaslArray2D<double> fx1Vec(grid.nX(), grid.nY()),
                fx2Vec(grid.nX(), grid.nY());
        double fx1, fx2;
        CaslLQRSystemDynamics sysDynamics;
        for (int i = 1; i <= grid.nX(); i++) {
            for (int j = 1; j <= grid.nY(); j++) {
                double xi = grid.x(i);
                double yj = grid.y(j);

                sysDynamics.LQRDynamics(xi, yj, fx1, fx2);
                fx1Vec(i, j) = fx1;
                fx2Vec(i, j) = fx2;
            }
        }

        // Hamiltonian set up:
        CaslHamiltonianLQR2D hamiltonianLQR(grid, fx1Vec, fx2Vec);
        CaslOptionPaddingWith boundaryConditions = withLinearExtrapolation;
        CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

        // HJSolver Initializing:
        CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianLQR, dt, currentTime, firstDerivativeScheme);

        // Define when to export to Matlab next:
        int exportCounter = 0;
        int exportFrequency = 40;
        int i = 0;
        double CFLNum = 0.9;

        // Exporting object
        CaslCppToMATLAB2D toMATLAB;

        // Inside the time-stepping loop
        CaslLQRSystemDynamics lqr;  // Add this before the loop

        while (currentTime < tFinal) {
            // Compute dt using findCFL function
            // dt = (CFLNum * (1.0 /HJSolver.findCFL(phi_n, currentTime)));

            if (currentTime + dt <= tFinal) {
                if (exportCounter % exportFrequency == 0) {
                    // Compute errors
                    double maxError = lqr.computeMaxError(grid, phi_n, currentTime);
                    double l2Error = lqr.computeL2Error(grid, phi_n, currentTime);

                    // Compute HJB error
                    double hjbError = 0.0;
                    for(int i = 1; i <= std::min(5, grid.nX()); i++) {
                        for(int j = 1; j <= std::min(5, grid.nY()); j++) {
                            hjbError = std::max(hjbError,
                                                std::abs(lqr.verifyHJB(grid.x(i), grid.y(j), currentTime)));
                        }
                    }

                    // Print errors to console
                    std::cout << "----------------------------------------" << std::endl;
                    std::cout << "At time " << currentTime << ":" << std::endl;
                    std::cout << "Maximum Error: " << std::scientific << std::setprecision(6) << maxError << std::endl;
                    std::cout << "L2 Error: " << std::scientific << std::setprecision(6) << l2Error << std::endl;
                    std::cout << "HJB Error: " << std::scientific << std::setprecision(6) << hjbError << std::endl;
                    std::cout << "----------------------------------------" << std::endl;

                    // Export time step data with all errors
                    exportTimeStep(currentTime, dt, i, phiFolder, testCase, grid, phi_n,
                                   toMATLAB, iMin, jMin, maxError, l2Error, hjbError);
                }
            }

            if (currentTime + dt > tFinal) {
                dt = tFinal - currentTime;
                currentTime = tFinal;

                // Compute final errors
                double maxError = lqr.computeMaxError(grid, phi_n, currentTime);
                double l2Error = lqr.computeL2Error(grid, phi_n, currentTime);

                // Compute final HJB error
                double hjbError = 0.0;
                for(int i = 1; i <= std::min(5, grid.nX()); i++) {
                    for(int j = 1; j <= std::min(5, grid.nY()); j++) {
                        hjbError = std::max(hjbError,
                                            std::abs(lqr.verifyHJB(grid.x(i), grid.y(j), currentTime)));
                    }
                }

                // Print final errors
                std::cout << "========== FINAL ERRORS ==========" << std::endl;
                std::cout << "Maximum Error: " << std::scientific << std::setprecision(6) << maxError << std::endl;
                std::cout << "L2 Error: " << std::scientific << std::setprecision(6) << l2Error << std::endl;
                std::cout << "HJB Error: " << std::scientific << std::setprecision(6) << hjbError << std::endl;
                std::cout << "=================================" << std::endl;

                exportTimeStep(currentTime, dt, i, phiFolder, testCase, grid, phi_n,
                               toMATLAB, iMin, jMin, maxError, l2Error, hjbError);
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
            exportCounter++;
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> simTime = endTime - startTime;

        double maxError = lqr.computeMaxError(grid, phi_n, currentTime);
        double l2Error = computeL2Error(phi_n, grid, lqr, currentTime);

        double hjbError = 0.0;
        for(int i = 1; i <= std::min(5, grid.nX()); i++) {
            for(int j = 1; j <= std::min(5, grid.nY()); j++) {
                hjbError = std::max(hjbError,
                                    std::abs(lqr.verifyHJB(grid.x(i), grid.y(j), currentTime)));
            }
        }

        // Add results including HJB error
        convergence.addResult(gridSize, grid.dx(), maxError, l2Error, hjbError, simTime.count());

    }


    // Display results
    convergence.computeConvergenceRates();

    return EXIT_SUCCESS;
}

double costFunction(double x, double xDot) {
    double cost = (1.0 / 2.0) * (pow(x, 2) + pow(xDot, 2));
    return cost;
}

void findMinMax2D(const CaslGrid2D& grid, const CaslArray2D<double>& phi,
                  double& minVal, double& maxVal,
                  int& iMin, int& jMin) {

    int nX, nY;
    nX = grid.nX(); nY = grid.nY();

    minVal = 100000000;
    maxVal = -100000000;

    // Find min and max values
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            double currentVal = phi(i, j);
            if (currentVal < minVal) {
                minVal = currentVal;
            }

            if (currentVal >= maxVal) {
                maxVal = currentVal;
            }
        }
    }

    // Find the first location of the minimum value
    bool isMinFound = false;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            if (phi(i, j) == minVal) {
                iMin = i;
                jMin = j;
                isMinFound = true;
                break;
            }
        }
        if (isMinFound) {
            break;
        }
    }
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

// Helper function for exporting
void exportTimeStep(double currentTime, double dt, int& i,
                    const std::string& phiFolder,
                    const std::string& testCase,
                    const CaslGrid2D& grid, const CaslArray2D<double>& phi_n,
                    CaslCppToMATLAB2D& toMATLAB,
                    int iMin, int jMin,
                    double maxError, double l2Error, double hjbError) {
    std::string dataFileName = phiFolder + "phi_" + std::to_string(i) + ".dat";
    std::string infoFileName = phiFolder + "simulation_info.txt";

    // Export data and info
    toMATLAB.exportDataToMatlab(grid, phi_n, dataFileName);

    // Open info file in append mode
    std::ofstream infoFile(infoFileName, std::ios::app);
    if (infoFile.is_open()) {
        infoFile << std::fixed;
        infoFile << "Export to Matlab at currentTime = " << currentTime << std::endl;
        infoFile << "dt: " << std::setprecision(16) << dt << std::endl;
        infoFile << "phi_" << i << ": " << std::setprecision(16) << phi_n(iMin, jMin) << std::endl;
        infoFile << "Maximum Error: " << std::scientific << std::setprecision(6) << maxError << std::endl;
        infoFile << "L2 Error: " << std::scientific << std::setprecision(6) << l2Error << std::endl;
        infoFile << "HJB Error: " << std::scientific << std::setprecision(6) << hjbError << std::endl;
        infoFile << "----------------------------------------" << std::endl;
        infoFile.close();
    }

    // Console output
    std::cout << "Export to Matlab at currentTime = " << currentTime << std::endl;
    std::cout << "dt: " << std::setprecision(16) << dt << std::endl;
    std::cout << "phi_" << i << ": " << std::setprecision(16) << phi_n(iMin, jMin) << std::endl;

    i++;
}

// Function to compute L2 error
double computeL2Error(const CaslArray2D<double>& numerical, const CaslGrid2D& grid,
                      const CaslLQRSystemDynamics& lqr, double currentTime) {
    double l2Error = 0.0;

    for(int i = 1; i <= grid.nX(); i++) {
        for(int j = 1; j <= grid.nY(); j++) {
            double x1 = grid.x(i);
            double x2 = grid.y(j);

            // Get exact solution at this point
            double exactVal = lqr.exactSolution(x1, x2, currentTime);

            // Compute squared error
            double error = pow(exactVal - numerical(i,j), 2);

            // Add to L2 error
            l2Error += error * grid.dx() * grid.dy(); // Include grid cell area
        }
    }

    return sqrt(l2Error); // Return square root for L2 norm
}

// Function to compute HJB residual at a point (optional helper function)
double computeHJBResidual(const CaslGrid2D& grid, const CaslArray2D<double>& phi_n,
                          const CaslLQRSystemDynamics& lqr, double currentTime,
                          int i, int j) {
    double x1 = grid.x(i);
    double x2 = grid.y(j);
    return std::abs(lqr.verifyHJB(x1, x2, currentTime));
}
#endif // CASL_LQR2D