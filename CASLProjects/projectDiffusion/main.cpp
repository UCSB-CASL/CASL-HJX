/**
 * @file main.cpp
 * @brief CASL-HJX Pure Diffusion Equation Solver (Heat Equation)
 *
 * @details This program solves the linear diffusion equation (heat equation) in 2D
 * using high-order numerical methods within the CASL-HJX framework. The diffusion
 * equation is fundamental to modeling heat conduction, mass transport, and various
 * other physical phenomena.
 *
 * Mathematical Model:
 * ∂φ/∂t = D∇²φ + f(x,y,t)
 *
 * where:
 * - φ(x,y,t) is the scalar field (temperature, concentration, etc.)
 * - D is the diffusion coefficient (thermal diffusivity, etc.)
 * - ∇²φ is the Laplacian operator (second-order spatial derivatives)
 * - f(x,y,t) is an optional source/sink term
 *
 * Numerical Methods:
 * - Spatial discretization: Central differences for Laplacian
 * - Time integration: Backward Euler (BTCS) - unconditionally stable
 * - Boundary conditions: Constant extrapolation
 * - Initial condition: Unit circle characteristic function
 *
 * Applications:
 * - Heat conduction in solids
 * - Mass diffusion in fluids
 * - Concentration spreading in chemical systems
 * - Image processing (diffusion filtering)
 * - Financial mathematics (option pricing)
 *
 * @author Faranak Rajabi
 * @date 2025-01-17
 * @version 1.0
 *
 * @copyright UC Santa Barbara - Computational Applied Sciences Laboratory
 *
 * Usage:
 *   ./projectDiffusion
 *
 * Output:
 *   Results are exported to __Output/ directory as .dat files
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

// Project-specific Hamiltonian (used for operator splitting framework)
#include "projectDiffusion_lib/CaslHamiltonianDiffusion2D.h"

using namespace std;

/**
 * @brief Heat source/sink function for the diffusion equation
 *
 * @details This function defines external heat sources or sinks in the diffusion
 * equation. For the homogeneous heat equation, this returns zero. Can be modified
 * to include:
 * - Spatially varying heat sources
 * - Time-dependent heating/cooling
 * - Nonlinear source terms
 * - Reaction terms for reaction-diffusion systems
 *
 * @param grid Reference to the computational grid
 * @param i Grid index in x-direction
 * @param j Grid index in y-direction
 * @param time Current simulation time
 * @param dimension Spatial dimension (1=x, 2=y)
 * @return Source term value at grid point (i,j) and time t
 */
double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);

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
 * @brief Main simulation function for diffusion equation solver
 *
 * @details Implements the complete numerical solution of the 2D diffusion equation
 * using implicit methods for unconditional stability:
 *
 * 1. Grid Setup: Uniform Cartesian grid on [-2, 2]²
 * 2. Initial Conditions: Unit circle characteristic function
 * 3. Time Integration: Fixed timestep with BTCS scheme
 * 4. Spatial Discretization: Central differences for Laplacian
 * 5. Boundary Treatment: Constant extrapolation
 * 6. Data Export: Every timestep for animation/analysis
 *
 * @return EXIT_SUCCESS on successful completion
 */
int main()
{
    // ========================================================================
    // SIMULATION HEADER AND SETUP
    // ========================================================================

    cout << "CASL-HJX: Pure Diffusion Equation Solver" << endl;
    cout << "=========================================" << endl;
    cout << "Mathematical Model: ∂φ/∂t = D∇²φ + f(x,y,t)" << endl;
    cout << "Numerical Method: Backward Euler (BTCS) - Unconditionally Stable" << endl;
    cout << endl;

    // ========================================================================
    // OUTPUT DIRECTORY MANAGEMENT
    // ========================================================================

    string projectName = "projectDiffusion";

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
        // Go up to CASLHJB2D level, then down to CASLProjects/projectDiffusion
        fullOutputDir = "../../../CASLProjects/projectDiffusion/__Output";
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

    // Domain: [-2, 2] × [-2, 2]
    // This domain is chosen to contain the initial circular region while
    // providing sufficient space to observe diffusion evolution
    int nX = 250, nY = 250;  // High resolution for smooth visualization
    double xMin = -2, xMax = 2;
    double yMin = -2, yMax = 2;

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

    double initialT = 0.0;    // Start time
    double finalT = 1.0;      // End time (sufficient for diffusion spreading)
    double current_time = initialT;
    double dt = 0.02;         // Fixed time step (stable for BTCS)

    cout << "Time Integration:" << endl;
    cout << "  Time interval: [" << initialT << ", " << finalT << "]" << endl;
    cout << "  Time step: Δt = " << dt << " (fixed)" << endl;
    cout << "  Total steps: " << static_cast<int>(finalT / dt) << endl;
    cout << "  Scheme: BTCS (unconditionally stable)" << endl;
    cout << endl;

    // ========================================================================
    // SOLUTION ARRAYS WITH PADDING
    // ========================================================================

    // Padding for boundary conditions and stencil operations
    const int nPadsX = 3;
    const int nPadsY = 3;

    // Solution arrays for time integration
    CaslArray2D<double> phi_n(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);      // φ^n
    CaslArray2D<double> phi_np1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);    // φ^{n+1}
    CaslArray2D<double> phi_np2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);    // φ^{n+2}
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);   // φ^{n+1/2}
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);   // φ^{n+3/2}

    // ========================================================================
    // PHYSICAL PARAMETERS
    // ========================================================================

    double diffusionC = 0.25; // Diffusion coefficient D
    double L = 1.0;           // Characteristic length scale
    int dimension = 1;        // X-direction (for source function)

    cout << "Physical Parameters:" << endl;
    cout << "  Diffusion coefficient: D = " << diffusionC << endl;
    cout << "  Characteristic length: L = " << L << endl;
    cout << "  Diffusion time scale: τ = L²/D = " << (L*L)/diffusionC << endl;
    cout << "  Dimensionless time: t* = t·D/L² = t·" << diffusionC/(L*L) << endl;
    cout << endl;

    // ========================================================================
    // ADVECTION COEFFICIENTS (ZERO FOR PURE DIFFUSION)
    // ========================================================================

    // Initialize advection coefficients (zero for pure diffusion)
    // These are required by the operator splitting framework
    CaslArray2D<double> c1(nX, nY);  // x-velocity (zero)
    CaslArray2D<double> c2(nX, nY);  // y-velocity (zero)

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            c1(i, j) = 0.0;  // No advection in x-direction
            c2(i, j) = 0.0;  // No advection in y-direction
        }
    }

    // ========================================================================
    // INITIAL CONDITIONS
    // ========================================================================

    cout << "Setting initial conditions..." << endl;

    // Unit circle characteristic function: φ(x,y,0) = χ_{x²+y²<1}
    // This represents an initial "hot spot" or concentration within unit circle
    double initialMass = 0.0;  // For conservation checking

    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            double x = grid.x(i);
            double y = grid.y(j);
            double radius_squared = x*x + y*y;

            if (radius_squared < 1.0) {
                phi_n(i, j) = 1.0;
                initialMass += 1.0;
            } else {
                phi_n(i, j) = 0.0;
            }
        }
    }

    cout << "  Profile: Unit circle χ_{x²+y²<1}" << endl;
    cout << "  Initial total mass: " << initialMass * grid.dx() * grid.dy() << endl;
    cout << "  Expected diffusion radius: r(t) ≈ √(4Dt)" << endl;
    cout << endl;

    // ========================================================================
    // NUMERICAL METHOD CONFIGURATION
    // ========================================================================

    // Solver options for maximum stability and accuracy
    CaslOptionPaddingWith boundaryCondition = withConstantExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
    CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
    CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

    cout << "Numerical Configuration:" << endl;
    cout << "  Second-order discretization: Central differences" << endl;
    cout << "  Time integration: BTCS (Backward Euler)" << endl;
    cout << "  Boundary conditions: Constant extrapolation" << endl;
    cout << "  Stability: Unconditionally stable (implicit method)" << endl;
    cout << endl;

    // ========================================================================
    // SOLVER INITIALIZATION
    // ========================================================================

    // Hamiltonian for operator splitting framework (zero for pure diffusion)
    CaslHamiltonianDiffusion2D hamiltonianForDiffusion(grid, c1, c2);

    // Hamilton-Jacobi solver (not used for pure diffusion, but required)
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianForDiffusion, dt, current_time, firstDerivativeScheme);

    // Laplacian solver for diffusion terms (main solver)
    CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, current_time, diffusionC,
                                                        secondOrderTermDirection, secondDerivativeScheme,
                                                        boundaryCondition, HJSolver, heatSourceFunction);

    cout << "Solvers initialized successfully." << endl;
    cout << endl;

    // ========================================================================
    // DATA EXPORT SETUP
    // ========================================================================

    CaslCppToMATLAB2D cppToMatlab;

    // Export initial condition
    string initialFile = fullOutputDir + "/phi_t0p00.dat";
    cppToMatlab.exportDataToMatlab(grid, phi_n, initialFile);
    cout << "Initial condition exported: " << initialFile << endl;
    cout << endl;

    // ========================================================================
    // MAIN TIME INTEGRATION LOOP
    // ========================================================================

    cout << "Starting time integration..." << endl;
    cout << "Progress: ";

    int timestepNumber = 0;
    int progressCounter = 0;
    double progressIncrement = finalT / 50;  // Progress updates
    double nextProgressTime = progressIncrement;

    while (current_time < finalT) {

        // ====================================================================
        // TIME STEP ADJUSTMENT
        // ====================================================================

        // Adjust final time step to hit final time exactly
        if (current_time + dt > finalT) {
            dt = finalT - current_time;
        }

        // ====================================================================
        // TVD-RK3 SCHEME (Applied to diffusion equation)
        // ====================================================================

        // Note: For pure diffusion, the TVD-RK3 is applied to the implicit
        // diffusion operator, providing high-order temporal accuracy

        // Step 1: φ^{n+1} = φ^n + Δt·D∇²φ^{n+1} (implicit)
        phi_n.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_n, phi_np1);  // Explicit part (zero for pure diffusion)

        phi_np1.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_np1, phi_np2);

        phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;
        phi_np12.fillPaddingPoints(boundaryCondition);
        HJSolver.eulerStep(phi_np12, phi_np32);

        phi_np1 = (1.0/3.0) * phi_n + (2.0/3.0) * phi_np32;

        // Apply implicit diffusion operator (main computation)
        LaplacianSolver.backwardTimeCentralSpacing(phi_n, phi_np1, 1e-30, 1e6);
        LaplacianSolver.backwardTimeCentralSpacing(phi_np1, phi_np2, 1e-30, 1e6);

        phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;
        LaplacianSolver.backwardTimeCentralSpacing(phi_np12, phi_np32, 1e-30, 1e6);

        phi_np1 = (1.0/3.0) * phi_n + (2.0/3.0) * phi_np32;

        // ====================================================================
        // TIME ADVANCEMENT
        // ====================================================================

        current_time += dt;
        phi_n = phi_np1;
        timestepNumber++;

        // ====================================================================
        // PROGRESS REPORTING
        // ====================================================================

        if (current_time >= nextProgressTime) {
            cout << ".";
            cout.flush();
            nextProgressTime += progressIncrement;
            progressCounter++;
        }

        // ====================================================================
        // DATA EXPORT (EVERY TIMESTEP)
        // ====================================================================

        // Generate filename with proper time formatting
        char timeStr[32];
        sprintf(timeStr, "%.2f", current_time);
        string timeString(timeStr);
        size_t pos = timeString.find('.');
        if (pos != string::npos) {
            timeString[pos] = 'p';
        }

        string filename = fullOutputDir + "/phi_t" + timeString + ".dat";
        cppToMatlab.exportDataToMatlab(grid, phi_n, filename);
    }

    cout << " Complete!" << endl;
    cout << endl;

    // ========================================================================
    // SIMULATION SUMMARY AND VERIFICATION
    // ========================================================================

    // Calculate final mass for conservation check
    double finalMass = 0.0;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            finalMass += phi_n(i, j);
        }
    }
    finalMass *= grid.dx() * grid.dy();

    cout << "Simulation completed successfully!" << endl;
    cout << "Summary:" << endl;
    cout << "  Final time: t = " << current_time << endl;
    cout << "  Total timesteps: " << timestepNumber << endl;
    cout << "  Initial mass: " << initialMass * grid.dx() * grid.dy() << endl;
    cout << "  Final mass: " << finalMass << endl;
    cout << "  Mass conservation error: " << abs(finalMass - initialMass * grid.dx() * grid.dy()) << endl;
    cout << "  Results directory: " << fullOutputDir << endl;

    // Theoretical diffusion radius at final time
    double theoreticalRadius = sqrt(4.0 * diffusionC * finalT);
    cout << "  Theoretical diffusion radius: r = " << theoreticalRadius << endl;

    // Verify files were created
    cout << endl;
    cout << "Verifying output files..." << endl;
    system(("ls -la " + fullOutputDir + "/ | head -10 && echo '...' && ls -la " + fullOutputDir + "/ | tail -5").c_str());

    cout << endl;
    cout << "Post-processing suggestions:" << endl;
    cout << "  1. Create diffusion animation with MATLAB/Python" << endl;
    cout << "  2. Verify mass conservation over time" << endl;
    cout << "  3. Compare with analytical Gaussian solution" << endl;
    cout << "  4. Study diffusion radius growth: r(t) ∝ √t" << endl;
    cout << "  5. Analyze heat equation eigenmodes" << endl;

    return EXIT_SUCCESS;
}

// ============================================================================
// FUNCTION IMPLEMENTATIONS
// ============================================================================

double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension)
{
    /**
     * @brief Homogeneous heat equation (no source terms)
     *
     * For the standard heat equation ∂φ/∂t = D∇²φ, no external sources
     * are applied. This function can be modified to include:
     * - Spatially varying heat sources: f(x,y)
     * - Time-dependent heating: f(t)
     * - Nonlinear sources: f(φ) for reaction-diffusion
     * - Gaussian sources: A·exp(-(x²+y²)/σ²)
     */
    return 0.0;

    // Example alternative source terms:
    //
    // Constant source everywhere:
    // return 0.1;
    //
    // Point source at origin:
    // double x = grid.x(i);
    // double y = grid.y(j);
    // return (abs(x) < 0.1 && abs(y) < 0.1) ? 10.0 : 0.0;
    //
    // Time-decaying Gaussian source:
    // double x = grid.x(i);
    // double y = grid.y(j);
    // return exp(-(x*x + y*y)/0.25) * exp(-time);
    //
    // Linear spatial variation:
    // double x = (dimension == 1) ? grid.x(i) : grid.y(j);
    // return x;
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