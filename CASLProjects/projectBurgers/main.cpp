/**
 * @file main.cpp
 * @brief CASL-HJX Burgers Equation Solver
 *
 * @details This program solves the viscous Burgers equation in 2D using high-order
 * numerical methods within the CASL-HJX framework. The Burgers equation serves as
 * a fundamental model for nonlinear wave propagation and shock formation in fluid
 * dynamics.
 *
 * Mathematical Model:
 * ∂u/∂t + u·∇u = ν∇²u
 *
 * where:
 * - u(x,y,t) is the velocity field
 * - ν is the kinematic viscosity coefficient
 * - The nonlinear term u·∇u represents convective acceleration
 * - The diffusion term ν∇²u provides viscous dissipation
 *
 * Numerical Methods:
 * - Spatial discretization: WENO5 (5th-order Weighted Essentially Non-Oscillatory)
 * - Time integration: TVD-RK3 (3rd-order Total Variation Diminishing Runge-Kutta)
 * - Viscous terms: Backward Euler with Central Differences (BTCS)
 * - Boundary conditions: Constant extrapolation
 *
 * Applications:
 * - Shock wave formation and propagation
 * - Turbulence modeling (simplified)
 * - Nonlinear wave dynamics
 * - Algorithm validation for hyperbolic PDEs
 *
 * @author Andrew Wang
 * @date 2025-01-14
 * @version 1.0
 *
 * @copyright UC Santa Barbara - Computational Applied Sciences Laboratory
 *
 * Usage:
 *   ./projectBurgers
 *
 * Output:
 *   Results are exported to __Output/projectBurgers/ directory as .dat files
 *   compatible with MATLAB visualization tools.
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>      // String operations
#include <filesystem>  // Directory creation (C++17)
#include <sys/stat.h>  // System file operations
#include <iomanip>     // Output formatting
#include <cstdio>      // C-style I/O
#include <algorithm>   // STL algorithms
#include <unistd.h>    // For getcwd function
#include <cstdlib>     // For realpath function

// CASL-HJX Core Libraries
#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"
#include "../../CASLCommonLibrary/DPMatrix2D.h"

// Project-specific Hamiltonian
#include "projectBurgers_lib/CaslHamiltonianBurgers2D.h"

using namespace std;

/**
 * @brief Source term function for the Burgers equation
 *
 * @details This function defines any external forcing or source terms in the
 * Burgers equation. For the standard homogeneous case, this returns zero.
 * Can be modified to include:
 * - External forcing terms
 * - Source/sink terms
 * - Time-dependent driving forces
 *
 * @param grid Reference to the computational grid
 * @param i Grid index in x-direction
 * @param j Grid index in y-direction
 * @param time Current simulation time
 * @param dimension Spatial dimension (1=x, 2=y)
 * @return Source term value at grid point (i,j) and time t
 */
double homogeneousFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);

/**
 * @brief Creates output directory with cross-platform compatibility
 *
 * @details Attempts to create the specified directory using multiple methods
 * for maximum compatibility across different operating systems and C++ standards.
 *
 * @param path Directory path to create
 */
void createOutputDirectory(const string& path);

/**
 * @brief Main simulation function for Burgers equation solver
 *
 * @details Implements the complete numerical solution of the 2D viscous Burgers
 * equation using operator splitting and high-order methods:
 *
 * 1. Grid Setup: Uniform Cartesian grid on [-π/2, π/2]²
 * 2. Initial Conditions: Gaussian profile exp(-(x² + y²))
 * 3. Time Integration: Adaptive timestep with CFL condition
 * 4. Spatial Discretization: WENO5 for hyperbolic terms
 * 5. Viscous Treatment: Implicit backward Euler
 * 6. Data Export: Regular output for visualization
 *
 * @return EXIT_SUCCESS on successful completion
 */
int main()
{
    // ========================================================================
    // SIMULATION HEADER AND SETUP
    // ========================================================================

    cout << "CASL-HJX: Burgers Equation Solver" << endl;
    cout << "==================================" << endl;
    cout << "Mathematical Model: ∂u/∂t + u·∇u = ν∇²u" << endl;
    cout << "Numerical Methods: WENO5 + TVD-RK3 + BTCS" << endl;
    cout << endl;

    // ========================================================================
    // OUTPUT DIRECTORY MANAGEMENT
    // ========================================================================

    string projectName = "projectBurgers";

    // Get current working directory for debugging
    char cwd[1024];
    string currentPath;
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        currentPath = string(cwd);
        cout << "Current working directory: " << currentPath << endl;
    }

    // Robust detection: check if we're in cmake-build-release path
    string fullOutputDir;

    if (currentPath.find("cmake-build-release") != string::npos) {
        // We're definitely in a build directory
        // Go up to CASLHJB2D level, then down to CASLProjects/projectBurgers
        fullOutputDir = "../../../CASLProjects/projectBurgers/__Output";
        cout << "Build environment: CMake build directory detected" << endl;
    } else if (currentPath.find("build") != string::npos) {
        // Generic build directory
        fullOutputDir = "../../__Output";
        cout << "Build environment: Generic build directory" << endl;
    } else {
        // Assume we're in source directory
        fullOutputDir = "__Output";
        cout << "Build environment: Source directory" << endl;
    }

    createOutputDirectory(fullOutputDir);
    cout << "Output directory: " << fullOutputDir << endl;
    cout << "Resolved to: ";

    // Show the absolute path of output directory
    char resolvedPath[1024];
    if (realpath(fullOutputDir.c_str(), resolvedPath) != NULL) {
        cout << resolvedPath << endl;
    } else {
        cout << "Path resolution failed" << endl;
    }
    cout << endl;

    // ========================================================================
    // COMPUTATIONAL DOMAIN SETUP
    // ========================================================================

    // Domain: [-π/2, π/2] × [-π/2, π/2]
    // This domain is chosen to demonstrate shock formation while maintaining
    // computational efficiency and avoiding boundary effects
    double xMin = -M_PI_2, xMax = M_PI_2;
    double yMin = -M_PI_2, yMax = M_PI_2;
    int nX = 250, nY = 250;  // High resolution for accurate shock capture

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    cout << "Computational Domain:" << endl;
    cout << "  x ∈ [" << xMin << ", " << xMax << "]" << endl;
    cout << "  y ∈ [" << yMin << ", " << yMax << "]" << endl;
    cout << "  Resolution: " << nX << " × " << nY << " points" << endl;
    cout << "  Grid spacing: Δx = " << grid.dx() << ", Δy = " << grid.dy() << endl;
    cout << endl;

    // ========================================================================
    // TIME INTEGRATION PARAMETERS
    // ========================================================================

    double tInitial = 0.0;    // Start time
    double tFinal = 4.0;      // End time (sufficient for shock development)
    double tCurrent = tInitial;
    double dt;                // Adaptive time step
    double nOutput = 200;     // Number of output files

    cout << "Time Integration:" << endl;
    cout << "  Time interval: [" << tInitial << ", " << tFinal << "]" << endl;
    cout << "  Output frequency: " << nOutput << " snapshots" << endl;
    cout << "  Time step: Adaptive (CFL-based)" << endl;
    cout << endl;

    // ========================================================================
    // SOLUTION ARRAYS WITH PADDING
    // ========================================================================

    // Padding for high-order stencils and boundary conditions
    const int nPadsX = 3;  // Required for WENO5 scheme
    const int nPadsY = 3;

    // Solution arrays for TVD-RK3 time integration
    CaslArray2D<double> un(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);      // u^n
    CaslArray2D<double> unp1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);    // u^{n+1}
    CaslArray2D<double> unp2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);    // u^{n+2}
    CaslArray2D<double> unp12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);   // u^{n+1/2}
    CaslArray2D<double> unp32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);   // u^{n+3/2}

    // ========================================================================
    // INITIAL CONDITIONS
    // ========================================================================

    cout << "Setting initial conditions..." << endl;

    // Gaussian initial profile: u(x,y,0) = exp(-(x² + y²))
    // This smooth initial condition will develop into shocks due to nonlinearity
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            double x = grid.x(i);
            double y = grid.y(j);
            un(i, j) = exp(-1.0 * (x * x + y * y));
        }
    }

    cout << "  Profile: Gaussian exp(-(x² + y²))" << endl;
    cout << "  Maximum value: " << un.maxAbs() << endl;
    cout << endl;

    // ========================================================================
    // NUMERICAL METHOD CONFIGURATION
    // ========================================================================

    // Physical parameters
    double diffusionC = 0.0;  // Inviscid Burgers (set to small value for viscous)

    // CFL-based time step calculation
    double dx = grid.dx(), dy = grid.dy();
    dt = 0.4 * std::min({dx / (un.maxAbs() + 1e-12), dy / (un.maxAbs() + 1e-12), dx, dy});

    cout << "Numerical Configuration:" << endl;
    cout << "  Viscosity coefficient: ν = " << diffusionC << endl;
    cout << "  Initial time step: Δt = " << dt << endl;
    cout << "  CFL safety factor: 0.4" << endl;
    cout << endl;

    // Solver options
    CaslOptionPaddingWith boundaryCondition = withConstantExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
    CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
    CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

    cout << "Numerical Schemes:" << endl;
    cout << "  Spatial discretization: WENO5 (5th-order)" << endl;
    cout << "  Time integration: TVD-RK3 (3rd-order)" << endl;
    cout << "  Viscous terms: BTCS (Backward Euler)" << endl;
    cout << "  Boundary conditions: Constant extrapolation" << endl;
    cout << endl;

    // ========================================================================
    // SOLVER INITIALIZATION
    // ========================================================================

    // Hamiltonian for Burgers equation: H = u·∇u
    CaslHamiltonianBurgers2D hamiltonian(grid, un, un);

    // Hamilton-Jacobi solver for hyperbolic terms
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, tCurrent, firstDerivativeScheme);

    // Laplacian solver for viscous terms
    CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, tCurrent, diffusionC,
                                                        secondOrderTermDirection, secondDerivativeScheme,
                                                        boundaryCondition, HJSolver, homogeneousFunction);

    cout << "Solvers initialized successfully." << endl;
    cout << endl;

    // ========================================================================
    // DATA EXPORT SETUP
    // ========================================================================

    CaslCppToMATLAB2D cppToMatlab;

    // Export initial condition
    string initialFile = fullOutputDir + "/un_0.dat";
    cppToMatlab.exportDataToMatlab(grid, un, initialFile);
    cout << "Initial condition exported: " << initialFile << endl;

    // Output management
    double dtGlobal = dt;                    // Store original time step
    bool isOutput = false;                   // Output flag
    int nextFileIdx = 1;                     // File counter
    double saveIncrement = tFinal / nOutput; // Time between saves
    double nextOutputTime = saveIncrement;   // Next output time

    cout << "Output interval: Δt_output = " << saveIncrement << endl;
    cout << endl;

    // ========================================================================
    // MAIN TIME INTEGRATION LOOP
    // ========================================================================

    cout << "Starting time integration..." << endl;
    cout << "Progress: ";

    int progressCounter = 0;
    double progressIncrement = tFinal / 50;  // Progress updates
    double nextProgressTime = progressIncrement;

    while (tCurrent < tFinal) {

        // ====================================================================
        // ADAPTIVE TIME STEPPING
        // ====================================================================

        // Adjust time step for output synchronization
        if (tCurrent + dt >= nextOutputTime) {
            dt = nextOutputTime - tCurrent;
            isOutput = true;
        }

        // Adjust time step for final time
        if (tCurrent + dt > tFinal) {
            dt = tFinal - tCurrent;
            isOutput = true;
        }

        // Update CFL constraint (solution-dependent)
        double maxSpeed = un.maxAbs();
        double dtCFL = 0.4 * std::min({dx / (maxSpeed + 1e-12), dy / (maxSpeed + 1e-12)});
        dt = std::min(dt, dtCFL);

        // ====================================================================
        // TVD-RK3 TIME INTEGRATION SCHEME
        // ====================================================================

        // Step 1: u^{n+1} = u^n + Δt·L(u^n)
        LaplacianSolver.backwardTimeCentralSpacing(un, unp1);

        // Step 2: u^{n+2} = u^{n+1} + Δt·L(u^{n+1})
        LaplacianSolver.backwardTimeCentralSpacing(unp1, unp2);

        // Step 3: u^{n+1/2} = 3/4·u^n + 1/4·u^{n+2}
        unp12 = 0.75 * un + 0.25 * unp2;
        LaplacianSolver.backwardTimeCentralSpacing(unp12, unp32);

        // Step 4: u^{n+1} = 1/3·u^n + 2/3·u^{n+3/2}
        unp1 = (1.0/3.0) * un + (2.0/3.0) * unp32;

        // ====================================================================
        // TIME ADVANCEMENT
        // ====================================================================

        tCurrent += dt;
        un = unp1;

        // ====================================================================
        // PROGRESS REPORTING
        // ====================================================================

        if (tCurrent >= nextProgressTime) {
            cout << ".";
            cout.flush();
            nextProgressTime += progressIncrement;
            progressCounter++;
        }

        // ====================================================================
        // DATA EXPORT
        // ====================================================================

        if (isOutput) {
            // Generate filename with proper time formatting
            char timeStr[32];
            sprintf(timeStr, "%.3f", tCurrent);
            string timeString(timeStr);
            size_t pos = timeString.find('.');
            if (pos != string::npos) {
                timeString[pos] = 'p';
            }

            string filename = fullOutputDir + "/un_" + to_string(nextFileIdx) + "_t" + timeString + ".dat";
            cppToMatlab.exportDataToMatlab(grid, un, filename);

            ++nextFileIdx;
            isOutput = false;
            dt = dtGlobal;  // Reset to global time step
            nextOutputTime += saveIncrement;
        }
    }

    cout << " Complete!" << endl;
    cout << endl;

    // ========================================================================
    // SIMULATION SUMMARY
    // ========================================================================

    cout << "Simulation completed successfully!" << endl;
    cout << "Summary:" << endl;
    cout << "  Final time: t = " << tCurrent << endl;
    cout << "  Total output files: " << nextFileIdx << endl;
    cout << "  Final maximum value: " << un.maxAbs() << endl;
    cout << "  Results directory: " << fullOutputDir << endl;

    // Verify files were created
    cout << "Verifying output files..." << endl;
    system(("ls -la " + fullOutputDir + "/").c_str());

    cout << endl;

    cout << "Post-processing suggestions:" << endl;
    cout << "  1. Visualize shock formation with MATLAB/Python" << endl;
    cout << "  2. Analyze conservation properties" << endl;
    cout << "  3. Study nonlinear wave steepening" << endl;
    cout << "  4. Compare with analytical solutions (if available)" << endl;

    return EXIT_SUCCESS;
}

// ============================================================================
// FUNCTION IMPLEMENTATIONS
// ============================================================================

double homogeneousFunction(CaslGrid2D &grid, int i, int j, double time, int dimension)
{
    /**
     * @brief Homogeneous source term (zero forcing)
     *
     * For the standard Burgers equation, no external forcing is applied.
     * This function can be modified to include:
     * - Spatially varying source terms
     * - Time-dependent forcing
     * - Stochastic perturbations
     */
    return 0.0;

    // Example alternative source terms:
    //
    // Gaussian source:
    // double x = grid.x(i);
    // double y = grid.y(j);
    // return 0.1 * exp(-(x*x + y*y)) * sin(2*M_PI*time);
    //
    // Linear source:
    // return 0.01 * grid.x(i) * exp(-time);
}

void createOutputDirectory(const string& path)
{
    /**
     * @brief Cross-platform directory creation with fallback methods
     *
     * Attempts multiple methods for maximum compatibility:
     * 1. C++17 std::filesystem (if available)
     * 2. System-specific calls (mkdir)
     * 3. Verification using stat()
     */

    #if __cplusplus >= 201703L
        // Method 1: C++17 filesystem (preferred)
        try {
            std::filesystem::create_directories(path);
        }
        catch (const std::exception& e) {
            cout << "Warning: Filesystem method failed: " << e.what() << endl;
            // Fall back to system call
            system(("mkdir -p " + path).c_str());
        }
    #else
        // Method 2: System call (portable)
        system(("mkdir -p " + path).c_str());
    #endif

    // Method 3: Verification
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        cout << "Warning: Could not access directory: " << path << endl;
    } else if (!(info.st_mode & S_IFDIR)) {
        cout << "Warning: Path exists but is not a directory: " << path << endl;
    }
}