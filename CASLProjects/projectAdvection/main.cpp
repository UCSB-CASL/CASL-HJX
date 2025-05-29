/**
 * @file main.cpp
 * @brief Linear Advection Equation Solver using CASL-HJX Framework
 *
 * MATHEMATICAL MODEL:
 * ∂φ/∂t + c₁∂φ/∂x + c₂∂φ/∂y = 0
 *
 * Where:
 * - φ(x,y,t) is the transported scalar field
 * - c₁, c₂ are the advection velocity components (constant)
 * - Domain: [0, 2π] × [0, 2π] with periodic boundary conditions
 * - Initial condition: φ(x,y,0) = cos(x) + sin(y)
 * - Exact solution: φ(x,y,t) = cos(x - c₁t) + sin(y - c₂t)
 *
 * NUMERICAL METHODS:
 * - Spatial discretization: WENO5 (Weighted Essentially Non-Oscillatory, 5th order)
 * - Time integration: TVD-RK3 (Total Variation Diminishing Runge-Kutta, 3rd order)
 * - Boundary conditions: Periodic (suitable for transport on torus)
 * - CFL condition: dt ≤ 0.4 * min(dx/|c₁|, dy/|c₂|, dx, dy)
 *
 * PHYSICAL SIGNIFICANCE:
 * Linear advection represents the transport of a scalar quantity (e.g., temperature,
 * concentration, or potential) by a prescribed velocity field. This is fundamental
 * in atmospheric modeling, oceanography, and many engineering applications.
 *
 * OUTPUT FORMAT:
 * - phi_t*.dat: Scalar field φ at different time instances
 * - Convergence analysis with exact solution comparison
 * - ASCII format suitable for MATLAB/Python visualization
 *
 * @author Faranak Rajabi
 * @date 2024-01-17
 * @version 2.0 - Professional version with comprehensive documentation and output
 *
 * USAGE:
 * 1. Compile: mkdir build && cd build && cmake .. && make
 * 2. Run: ./projectAdvection
 * 3. Results: Check __Output/ folder for phi_t*.dat files
 * 4. Visualize: Use MATLAB/Python scripts to create animations
 *
 * DEPENDENCIES:
 * - CASL-HJX Common Library (CaslGrid2D, CaslArray2D, etc.)
 * - C++11 or later
 * - CMake build system
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <unistd.h>
#include <sys/stat.h>

// CASL-HJX Core Libraries
#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"
#include "../../CASLCommonLibrary/DPMatrix2D.h"

// Project-Specific Library
#include "projectAdvection_lib/CaslHamiltonianAdvection2D.h"

using namespace std;

/**
 * @brief Creates output directory with proper error handling
 * @param dirPath Path to the directory to create
 * @return true if directory exists or was created successfully
 */
bool createOutputDirectory(const string& dirPath) {
    struct stat st = {0};
    if (stat(dirPath.c_str(), &st) == -1) {
        #ifdef _WIN32
            return _mkdir(dirPath.c_str()) == 0;
        #else
            return mkdir(dirPath.c_str(), 0755) == 0;
        #endif
    }
    return true; // Directory already exists
}

/**
 * @brief Generates timestep filename with proper formatting
 * @param outputDir Base output directory path
 * @param timeStep Current time step number
 * @param currentTime Current simulation time
 * @return Formatted filename string
 */
string generateFilename(const string& outputDir, int timeStep, double currentTime) {
    string filename = outputDir + "/phi_t";
    if (timeStep == 0) {
        filename += "0.dat";
    } else {
        // Format time as t1p234 for t=1.234 (replace decimal with 'p')
        string timeStr = to_string(currentTime);
        size_t dotPos = timeStr.find('.');
        if (dotPos != string::npos) {
            timeStr[dotPos] = 'p';
        }
        // Remove trailing zeros
        timeStr.erase(timeStr.find_last_not_of('0') + 1, string::npos);
        if (timeStr.back() == 'p') timeStr.pop_back();

        filename += timeStr + ".dat";
    }
    return filename;
}

/**
 * @brief Main solver routine for linear advection equation
 *
 * Implements the complete solution pipeline:
 * 1. Grid and domain setup
 * 2. Initial condition specification
 * 3. Velocity field definition
 * 4. Numerical solver configuration
 * 5. Time marching with TVD-RK3
 * 6. Data export and error analysis
 *
 * @return EXIT_SUCCESS on successful completion
 */
int main()
{
    // ====================================================================
    // SECTION 1: SIMULATION PARAMETERS AND DOMAIN SETUP
    // ====================================================================

    cout << "CASL-HJX: Linear Advection Equation Solver" << endl;
    cout << "==========================================" << endl;
    cout << "Mathematical Model: ∂φ/∂t + c₁∂φ/∂x + c₂∂φ/∂y = 0" << endl;
    cout << "Numerical Methods: WENO5 + TVD-RK3 + Periodic BC" << endl;

    // Physical domain: [0, 2π] × [0, 2π] (periodic torus)
    const double xMin = 0.0, xMax = 2.0 * M_PI;
    const double yMin = 0.0, yMax = 2.0 * M_PI;
    const int nX = 80, nY = 80;  // Grid resolution (80×80 = 6400 points)

    // Temporal domain
    const double tInitial = 0.0, tFinal = 1.0;
    const int nOutput = 50;  // Export frequency (20 ms intervals)

    cout << "Domain: [" << xMin << ", " << xMax << "] × [" << yMin << ", " << yMax << "]" << endl;
    cout << "Grid: " << nX << "×" << nY << " points" << endl;
    cout << "Time: t ∈ [" << tInitial << ", " << tFinal << "]" << endl;

    // ====================================================================
    // SECTION 2: DIRECTORY DETECTION AND OUTPUT SETUP
    // ====================================================================

    // Smart directory detection for cross-platform compatibility
    char cwd[1024];
    string outputDir;

    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        cout << "Current working directory: " << cwd << endl;

        string currentPath(cwd);
        if (currentPath.find("cmake-build-release") != string::npos) {
            cout << "Build environment: CMake build directory detected" << endl;
            outputDir = "../../../CASLProjects/projectAdvection/__Output";
        } else if (currentPath.find("build") != string::npos) {
            cout << "Build environment: Standard build directory" << endl;
            outputDir = "../__Output";
        } else {
            cout << "Build environment: Source directory" << endl;
            outputDir = "__Output";
        }
    } else {
        cout << "Build environment: Default (source directory assumed)" << endl;
        outputDir = "__Output";
    }

    cout << "Output directory: " << outputDir << endl;

    // Create output directory structure
    if (createOutputDirectory(outputDir)) {
        cout << "Output directory ready: " << outputDir << endl;
    } else {
        cerr << "Error: Could not create output directory: " << outputDir << endl;
        return EXIT_FAILURE;
    }

    // ====================================================================
    // SECTION 3: GRID INITIALIZATION AND ARRAY ALLOCATION
    // ====================================================================

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);
    double dx = grid.dx(), dy = grid.dy();

    cout << "Grid spacing: dx = " << dx << ", dy = " << dy << endl;

    // Padding for high-order WENO5 scheme (requires 3 ghost points)
    const int nPads = 3;

    // Solution arrays for TVD-RK3 time stepping
    CaslArray2D<double> phiN(nX, nY, nPads);      // Current solution φⁿ
    CaslArray2D<double> phiNp1(nX, nY, nPads);    // Intermediate solution φⁿ⁺¹
    CaslArray2D<double> phiNp2(nX, nY, nPads);    // Second RK stage
    CaslArray2D<double> phiNp12(nX, nY, nPads);   // Half-step combination
    CaslArray2D<double> phiNp32(nX, nY, nPads);   // Three-halves step

    // ====================================================================
    // SECTION 4: VELOCITY FIELD AND INITIAL CONDITIONS
    // ====================================================================

    // Define constant advection velocity field
    // c₁ = 0 (no x-direction transport), c₂ = 1 (unit velocity in y-direction)
    CaslArray2D<double> c1(nX, nY), c2(nX, nY);
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            c1(i, j) = 0.0;  // No x-advection
            c2(i, j) = 1.0;  // Unit y-advection
        }
    }

    cout << "Velocity field: c₁ = " << c1(1,1) << ", c₂ = " << c2(1,1) << " (uniform)" << endl;

    // Set initial condition: φ(x,y,0) = cos(x) + sin(y)
    // This creates a smooth, periodic pattern suitable for transport analysis
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            phiN(i, j) = cos(grid.x(i)) + sin(grid.y(j));
        }
    }

    cout << "Initial condition: φ(x,y,0) = cos(x) + sin(y)" << endl;

    // ====================================================================
    // SECTION 5: NUMERICAL SOLVER CONFIGURATION
    // ====================================================================

    // Configure Hamiltonian for linear advection
    CaslHamiltonianAdvection2D hamiltonian(grid, c1, c2);

    // Numerical method selections
    CaslOptionPaddingWith boundaryConditions = withPeriodicCondition;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

    // CFL-limited time step for numerical stability
    // CFL = max(|c₁|Δt/Δx + |c₂|Δt/Δy) ≤ 0.4 (safety factor)
    double dt = 0.4 * min({dx / max(c1.maxAbs(), 1e-12),
                           dy / max(c2.maxAbs(), 1e-12),
                           dx, dy});

    cout << "Time step: dt = " << scientific << setprecision(3) << dt
         << " (CFL-limited)" << endl;

    // Initialize Hamilton-Jacobi solver
    double currentTime = tInitial;
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, currentTime,
                                         firstDerivativeScheme);

    // ====================================================================
    // SECTION 6: DATA EXPORT INITIALIZATION
    // ====================================================================

    CaslCppToMATLAB2D cppToMatlab;

    // Export initial condition
    string initialFile = generateFilename(outputDir, 0, currentTime);
    cout << "Exporting initial condition: " << initialFile << endl;
    cppToMatlab.exportDataToMatlab(grid, phiN, initialFile);

    // Setup for time-based exports
    double saveIncrement = tFinal / static_cast<double>(nOutput);
    double nextOutputTime = saveIncrement;
    int fileIndex = 1;

    cout << "Export frequency: every " << saveIncrement << " time units" << endl;
    cout << "Starting time integration..." << endl;

    // ====================================================================
    // SECTION 7: MAIN TIME INTEGRATION LOOP (TVD-RK3)
    // ====================================================================

    /**
     * Third-order Total Variation Diminishing Runge-Kutta scheme:
     *
     * Stage 1: φ⁽¹⁾ = φⁿ + Δt L(φⁿ)
     * Stage 2: φ⁽²⁾ = φ⁽¹⁾ + Δt L(φ⁽¹⁾)
     * Stage 3: φ̃ = (3/4)φⁿ + (1/4)φ⁽²⁾
     * Stage 4: φ⁽³⁾ = φ̃ + Δt L(φ̃)
     * Final:   φⁿ⁺¹ = (1/3)φⁿ + (2/3)φ⁽³⁾
     *
     * Where L(φ) = -c₁∂φ/∂x - c₂∂φ/∂y (advection operator)
     */

    double originalDt = dt;
    bool shouldExport = false;
    int stepCount = 0;

    while (currentTime < tFinal) {
        // Adaptive time step near output times and final time
        if (currentTime + dt >= nextOutputTime || currentTime + dt >= tFinal) {
            dt = min(nextOutputTime - currentTime, tFinal - currentTime);
            shouldExport = true;
        }

        // TVD-RK3 Stage 1: Forward Euler step
        phiN.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phiN, phiNp1);

        // TVD-RK3 Stage 2: Second Euler step
        phiNp1.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phiNp1, phiNp2);

        // TVD-RK3 Stage 3: Convex combination
        phiNp12 = 0.75 * phiN + 0.25 * phiNp2;

        // TVD-RK3 Stage 4: Third Euler step
        phiNp12.fillPaddingPoints(boundaryConditions);
        HJSolver.eulerStep(phiNp12, phiNp32);

        // TVD-RK3 Final: Second convex combination
        phiNp1 = (1.0/3.0) * phiN + (2.0/3.0) * phiNp32;

        // Update solution and time
        currentTime += dt;
        phiN = phiNp1;
        stepCount++;

        // Progress reporting (every 10 steps or at export times)
        if (stepCount % 10 == 0 || shouldExport) {
            cout << "t = " << fixed << setprecision(3) << currentTime
                 << " (step " << stepCount << ")" << endl;
        }

        // Export data at specified intervals
        if (shouldExport) {
            string filename = generateFilename(outputDir, fileIndex, currentTime);
            cout << "Exporting solution: " << filename << endl;
            cppToMatlab.exportDataToMatlab(grid, phiN, filename);

            fileIndex++;
            nextOutputTime += saveIncrement;
            shouldExport = false;
            dt = originalDt;  // Restore original time step
        }
    }

    // ====================================================================
    // SECTION 8: ERROR ANALYSIS AND VALIDATION
    // ====================================================================

    cout << "\nPerforming convergence analysis..." << endl;

    // Calculate L∞ error against analytical solution
    // Exact solution: φ_exact(x,y,t) = cos(x - c₁t) + sin(y - c₂t)
    double maxError = 0.0;
    double l2Error = 0.0;
    int totalPoints = 0;

    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            double xShifted = grid.x(i) - c1(i,j) * currentTime;
            double yShifted = grid.y(j) - c2(i,j) * currentTime;
            double exactSolution = cos(xShifted) + sin(yShifted);
            double pointError = abs(phiN(i,j) - exactSolution);

            maxError = max(maxError, pointError);
            l2Error += pointError * pointError;
            totalPoints++;
        }
    }

    l2Error = sqrt(l2Error / totalPoints);

    // ====================================================================
    // SECTION 9: RESULTS SUMMARY AND OUTPUT
    // ====================================================================

    cout << "\nSimulation completed successfully!" << endl;
    cout << "=================================" << endl;
    cout << "Grid resolution: " << nX << "×" << nY << " points" << endl;
    cout << "Total time steps: " << stepCount << endl;
    cout << "Final time: " << fixed << setprecision(6) << currentTime << endl;
    cout << "Files exported: " << fileIndex << " snapshots" << endl;
    cout << "\nConvergence Analysis:" << endl;
    cout << "L∞ error: " << scientific << setprecision(6) << maxError << endl;
    cout << "L² error: " << scientific << setprecision(6) << l2Error << endl;
    cout << "\nResults exported to: " << outputDir << endl;

    // Verify file creation
    cout << "\nCreated files:" << endl;
    string listCommand = "ls -la " + outputDir + "/";
    int result = system(listCommand.c_str());

    if (result != 0) {
        cout << "Note: Could not list output files (ls command failed)" << endl;
        cout << "Please check the output directory manually: " << outputDir << endl;
    }

    return EXIT_SUCCESS;
}