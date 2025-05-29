//
// Created by Faranak Rajabi at 01/17/24
// Revised for proper output directory creation
//

#include <cmath>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <filesystem>  // For directory creation (C++17)
#include <sys/stat.h>  // Alternative for older systems

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"
#include "../../CASLCommonLibrary/DPMatrix2D.h"

#include "projectAdvectionDiffusion_lib/CaslHamiltonianDiffusion2D.h"

using namespace std;

double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);
void createOutputDirectory(const string& path);

int main()
{
    cout << "CASL-HJX: Advection-Diffusion Solver" << endl;
    cout << "=======================================" << endl;

    // Get current working directory
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        cout << "Current directory: " << cwd << endl;
    }

    // Create output directory in source folder - use absolute path
    string sourceDir = "/Users/faranakrajabi/CLionProjects/CASLHJB2D/CASLProjects/projectAdvectionDiffusion";
    string outputDir = sourceDir + "/__Output";
    createOutputDirectory(outputDir);
    cout << "Output directory: " << outputDir << endl;

    createOutputDirectory(outputDir);
    cout << "Output directory ready: " << outputDir << endl;

    // Domain boundaries
    double xMin = -2, xMax = 2,
           yMin = -2, yMax = 2;

    // Grid set up:
    int nX = 250, nY = 250;

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);
    cout << "Grid initialized: " << nX << "×" << nY << endl;

    // Time parameters
    double initialT = 0.0, finalT = 1.0, current_time = initialT;
    double dt = 0.02;
    cout << "Time setup: t ∈ [" << initialT << ", " << finalT << "], dt = " << dt << endl;

    // Phi arrays set up with padding:
    const int nPadsX = 3;
    const int nPadsY = 3;

    CaslArray2D<double> phi_n(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    // Physical parameters
    double diffusionC = 0.25;  // Diffusion coefficient
    double L = 1.0;
    int dimension = 1; // X-direction

    // Advection coefficients (can be modified for different problems)
    CaslArray2D<double> c1(nX, nY);  // x-velocity
    CaslArray2D<double> c2(nX, nY);  // y-velocity

    // Initialize advection field (modify as needed)
    for (int i = 1; i <= nX; i++)
    {
        for (int j = 1; j <= nY; j++)
        {
            // Example: constant advection
            c1(i, j) = 1.0;  // Advection in x-direction
            c2(i, j) = 0.0;  // No advection in y-direction

            // Alternative: rotating flow field
            // c1(i, j) = -grid.y(j);
            // c2(i, j) = grid.x(i);
        }
    }

    // Initial condition: circular blob
    cout << "Setting initial condition..." << endl;
    for (int i = 1; i <= nX; ++i)
    {
        for (int j = 1; j <= nY; ++j)
        {
            double radius_squared = pow(grid.x(i), 2) + pow(grid.y(j), 2);
            if (radius_squared < 1.0)
            {
                phi_n(i, j) = 1.0;
            }
            else
            {
                phi_n(i, j) = 0.0;
            }
        }
    }

    // Export initial condition
    CaslCppToMATLAB2D cppToMatlab;
    string initialFile = outputDir + "/phi_t0.dat";
    cppToMatlab.exportDataToMatlab(grid, phi_n, initialFile);
    cout << "Initial condition exported: " << initialFile << endl;

    // Solver configuration
    CaslOptionPaddingWith boundaryCondition = withConstantExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
    CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
    CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

    // Initialize solvers
    CaslHamiltonianDiffusion2D hamiltonianForDiffusion(grid, c1, c2);
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianForDiffusion, dt, current_time, firstDerivativeScheme);
    CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, current_time, diffusionC,
                                                       secondOrderTermDirection, secondDerivativeScheme,
                                                       boundaryCondition, HJSolver, heatSourceFunction);

    cout << "Solvers initialized successfully" << endl;
    cout << "Starting time integration..." << endl;

    int timestepNumber = 0;
    int exportFrequency = 5; // Export every 5 time steps

    // Time integration loop
    while (current_time < finalT)
    {
        // Adjust final time step
        if (current_time + dt > finalT)
        {
            dt = finalT - current_time;
        }

        // TVD-RK3 scheme for advection (first-order terms)
        phi_n.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_n, phi_np1);

        phi_np1.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_np1, phi_np2);

        phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;
        phi_np12.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_np12, phi_np32);

        phi_np1 = (1.0/3.0) * phi_n + (2.0/3.0) * phi_np32;

        // Implicit treatment for diffusion (second-order terms)
        LaplacianSolver.backwardTimeCentralSpacing(phi_n, phi_np1, 1e-30, 1e6);
        LaplacianSolver.backwardTimeCentralSpacing(phi_np1, phi_np2, 1e-30, 1e6);

        phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;
        LaplacianSolver.backwardTimeCentralSpacing(phi_np12, phi_np32, 1e-30, 1e6);

        phi_np1 = (1.0/3.0) * phi_n + (2.0/3.0) * phi_np32;

        // Update for next time step
        current_time += dt;
        phi_n = phi_np1;
        timestepNumber++;

        // Progress output
        if (timestepNumber % 10 == 0)
        {
            cout << "t = " << fixed << setprecision(3) << current_time
                 << " (step " << timestepNumber << ")" << endl;
        }

        // Export data periodically
        if (timestepNumber % exportFrequency == 0 || current_time >= finalT)
        {
            string filename = outputDir + "/phi_t" + to_string(static_cast<int>(current_time*10)) + ".dat";
            cppToMatlab.exportDataToMatlab(grid, phi_n, filename);
        }
    }

    // Export final solution
    string finalFile = outputDir + "/phi_final.dat";
    cppToMatlab.exportDataToMatlab(grid, phi_n, finalFile);

    cout << "Simulation completed successfully!" << endl;
    cout << "Results exported to: " << outputDir << "/" << endl;
    cout << "Total time steps: " << timestepNumber << endl;
    cout << "Final time: " << current_time << endl;

    return EXIT_SUCCESS;
}

void createOutputDirectory(const string& path)
{
    // Method 1: Using C++17 filesystem (preferred if available)
    #if __cplusplus >= 201703L
        try {
            std::filesystem::create_directories(path);
        }
        catch (const std::exception& e) {
            cout << "Warning: Could not create directory using filesystem: " << e.what() << endl;
            // Fall back to system call
            system(("mkdir -p " + path).c_str());
        }
    #else
        // Method 2: Using system call (more portable)
        system(("mkdir -p " + path).c_str());
    #endif

    // Method 3: Check if directory exists (verification)
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
    {
        cout << "Warning: Could not access " << path << endl;
    }
    else if (!(info.st_mode & S_IFDIR))
    {
        cout << "Warning: " << path << " is not a directory" << endl;
    }
}

double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension)
{
    // Homogeneous case (no source term)
    return 0.0;

    // Example non-homogeneous source term:
    // double x = (dimension == 1) ? grid.x(i) : grid.y(j);
    // return x * exp(-time);  // Decaying source

    // Example Gaussian source:
    // double x = grid.x(i);
    // double y = grid.y(j);
    // return exp(-(x*x + y*y)) * sin(time);
}