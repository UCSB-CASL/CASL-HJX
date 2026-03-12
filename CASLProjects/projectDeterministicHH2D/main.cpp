// main_cpu_hh2D_validate.cpp
// CPU 2D HJB solver — modified to match GPU 4D print format exactly.
//
// Changes from original:
//   1. sigma = 0.01  (was 0.001 — must be > dx² for G=40)
//   2. nX=nY=40  (was 40 in original, kept)
//   3. Print format matches GPU NaN-check lines so you can compare directly
//   4. Terminal cost min/max printed immediately after IC setup
//   5. Per-step prints at steps 1,2,3,5,10,25,50,100,250,500,1000 matching GPU
//
// Expected relationship to GPU 4D output:
//   phi_4D(x_pl,y_pl,x_pl,y_pl) = 2 * phi_2D(x_pl,y_pl)
//   because phaseless terminal cost is additive over both neurons.
//   So GPU min=18.28 → CPU min should be ≈ 9.14.

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
#include <set>

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
void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName);
void ensureDirectoryExists(const std::string& directoryPath);
void printPhiStats(const CaslGrid2D& grid, const CaslArray2D<double>& phi,
                   const std::string& label, int nX, int nY, int scaleFac);

int main() {
    auto start_time = chrono::high_resolution_clock::now();

    string projectPath  = fs::current_path().string() + "/";
    string outputFolder = projectPath + "__Output/";
    string phiFolder    = outputFolder + "phi/";
    string uStarFolder  = outputFolder + "uStar/";

    ensureDirectoryExists(outputFolder);
    ensureDirectoryExists(phiFolder);
    ensureDirectoryExists(uStarFolder);

    // ================================================================
    // Parameters — matched to GPU 4D validation run
    // ================================================================
    const int    scaleFac = Ks;
    const double uMax     = 10.0;   // = Ib
    const double vSpike   = 44.8;
    const double nSpike   = 0.459;

    // Grid — G=40 matching GPU
    const double vMin = -100.0, vMax = 100.0;
    const double xMin = (1.0 / scaleFac) * vMin;   // -1.0
    const double xMax = (1.0 / scaleFac) * vMax;   //  1.0
    const double yMin = 0.0, yMax = 1.0;
    const int nX = 40, nY = 40;

    CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

    printf("============================================================\n");
    printf("  CPU 2D HJB  G=%d  uMax=%.0f  Tf=7 ms  sigma2=0.01  cfl=0.9  bc=linear\n", nX, uMax);
    printf("  domain: xMin=%.0f xMax=%.0f yMin=%.0f yMax=%.0f\n", xMin, xMax, yMin, yMax);
    printf("  phaseless: V_pl=-59.6 mV  n_pl=0.403\n");
    printf("  dx=%.5f  dy=%.5f\n", grid.dx(), grid.dy());
    printf("============================================================\n\n");

    const int nPadsX = 3, nPadsY = 3;
    CaslArray2D<double> phi_n   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_1   (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np1 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np2 (nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
    CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

    double tInitial = 0.0, tFinal = 7.0, currentTime = tInitial;

    // ================================================================
    // Terminal cost
    // ================================================================
    const double v_pl  = -59.6;
    const double n_pl  = 0.403;
    const double gama  = 1000.0;
    const double xTarg = (1.0 / scaleFac) * v_pl;   // -0.596
    const double yTarg = n_pl;
    const double sigma = 0.01;   // sigma2 — must be > dx² for G=40

    for (int i = 1; i <= nX; i++)
        for (int j = 1; j <= nY; j++)
            phi_n(i, j) = calculateIC(grid.x(i), grid.y(j), gama, xTarg, yTarg, sigma);

    // Print matching GPU "0-terminal-cost" line
    printPhiStats(grid, phi_n, "0-terminal-cost", nX, nY, scaleFac);

    // Find min cell (for tracking throughout solve)
    double minPhi_0 = minCalculator(grid, phi_n);
    int iMin=1, jMin=1;
    for (int i=1;i<=nX;i++)
        for (int j=1;j<=nY;j++)
            if (phi_n(i,j) == minPhi_0){ iMin=i; jMin=j; }

    // ================================================================
    // Dynamics
    // ================================================================
    CaslArray2D<double> fxVec(grid.nX(), grid.nY());
    CaslArray2D<double> fyVec(grid.nX(), grid.nY());

    for (int i = 1; i <= grid.nX(); i++)
        for (int j = 1; j <= grid.nY(); j++){
            const double v = scaleFac * grid.x(i);
            const double n = grid.y(j);
            fxVec(i,j) = (1.0/scaleFac) * fv(v,n);
            fyVec(i,j) = fn(v,n);
        }

    CaslHamiltonianHHModel hamiltonian(grid, fxVec, fyVec);
    CaslOptionPaddingWith boundaryConditions = withLinearExtrapolation;
    CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

    double CFLNum = 0.9;
    double dt = 0.0;
    CaslArray2D<double> uStar(nX, nY);
    CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, currentTime, firstDerivativeScheme);

    // Steps at which to print (matching GPU NaN-check schedule)
    set<int> print_steps = {1,2,3,5,10,25,50,100,250,500,1000,2000,5000,10000,15000,20000};

    int exportCounter = 0;
    int iExport = 0;
    int stepCounter = 1;

    printf("\nStarting backward HJB solve...\n");
    printf("    step     t(ms)   %%done   wall(s)\n");
    printf("------------------------------------------\n");

    auto t_solve0 = chrono::high_resolution_clock::now();

    while (currentTime < tFinal) {
        // const double cflVal = HJSolver.findCFL(phi_n, currentTime);
        // dt = CFLNum * (1.0 / cflVal);
        dt = 4.17e-04;

        if (currentTime + dt > tFinal) dt = tFinal - currentTime;

        HJSolver.computeOptimalControl(phi_n, scaleFac, uMax, uStar);

        // Export every 120 steps
        if (exportCounter % 120 == 0 && currentTime + dt <= tFinal){
            string phi_fn   = phiFolder   + "phi_"   + to_string(iExport) + ".dat";
            string ustar_fn = uStarFolder + "uStar_" + to_string(iExport) + ".dat";
            exportToMatlab(grid, phi_n,  phi_fn);
            exportToMatlab(grid, uStar, ustar_fn);
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

        phi_np1 = (1.0/3.0) * phi_n + (2.0/3.0) * phi_np32;

        currentTime += dt;
        phi_n = phi_np1;

        // Print matching GPU format
        if (print_steps.count(stepCounter)){
            double wall = chrono::duration<double>(
                chrono::high_resolution_clock::now()-t_solve0).count();
            printPhiStats(grid, phi_n,
                "step-" + to_string(stepCounter) +
                " t=" + to_string(currentTime).substr(0,6) + "ms",
                nX, nY, scaleFac);
            printf("  %6d   %7.4f   %4.0f%%   %.4f   dt=%.2e\n",
                stepCounter, currentTime,
                100.0*currentTime/tFinal, wall, dt);
        }

        stepCounter++;
        exportCounter++;
    }

    // Final export
    {
        string phi_fn   = phiFolder   + "phi_"   + to_string(iExport) + ".dat";
        string ustar_fn = uStarFolder + "uStar_" + to_string(iExport) + ".dat";
        exportToMatlab(grid, phi_n,  phi_fn);
        exportToMatlab(grid, uStar, ustar_fn);
    }

    // ================================================================
    // Final diagnostics — matching GPU format
    // ================================================================
    printf("\n--- Final diagnostics ---\n");
    printPhiStats(grid, phi_n, "final", nX, nY, scaleFac);
    printf("dt: %.4e\n", dt);
    printf("uStar @ phaseless cell (%d,%d): %.4f\n", iMin-1, jMin-1, uStar(iMin,jMin));
    printf("phi @ phaseless cell:  %.4f\n", phi_n(iMin,jMin));
    printf("\nNOTE: GPU phi_4D at diagonal = 2 * this value (additive terminal cost)\n");
    printf("  GPU min expected: %.4f  (= 2 * %.4f)\n",
           2.0*phi_n(iMin,jMin), phi_n(iMin,jMin));

    auto end_time = chrono::high_resolution_clock::now();
    auto sim_tim  = chrono::duration_cast<chrono::minutes>(end_time - start_time);
    printf("\nTime for the whole simulation: %ld minutes\n", sim_tim.count());

    return 0;
}

// ================================================================
// printPhiStats — matches GPU NaN-check line format
// ================================================================
void printPhiStats(const CaslGrid2D& grid, const CaslArray2D<double>& phi,
                   const string& label, int nX, int nY, int scaleFac)
{
    double minVal =  1e99, maxVal = -1e99;
    int iMin=1,jMin=1,iMax=1,jMax=1;
    for(int i=1;i<=nX;i++) for(int j=1;j<=nY;j++){
        double v=phi(i,j);
        if(v<minVal){minVal=v;iMin=i;jMin=j;}
        if(v>maxVal){maxVal=v;iMax=i;jMax=j;}
    }
    printf("[CPU-check %-30s]\n", label.c_str());
    printf("   phi: min=+%.4g @ (%d,%d) V=%.1fmV n=%.3f\n",
        minVal, iMin-1, jMin-1,
        (double)scaleFac*grid.x(iMin), grid.y(jMin));
    printf("        max=+%.4g @ (%d,%d)   phi[mid]=%.4g\n",
        maxVal, iMax-1, jMax-1,
        phi(nX/2, nY/2));
}

// ================================================================
// Helpers
// ================================================================
double calculateIC(double x, double y, double gama, double xTarg, double yTarg, double sigma){
    const double r2 = (pow(x-xTarg,2) + pow(y-yTarg,2)) / sigma;
    return gama * (1.0 - exp(-r2));
}

double minCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un){
    double v = numeric_limits<double>::infinity();
    for(int i=1;i<=grid.nX();i++) for(int j=1;j<=grid.nY();j++) v=min(v,un(i,j));
    return v;
}

double maxCalculator(const CaslGrid2D& grid, const CaslArray2D<double>& un){
    double v = -numeric_limits<double>::infinity();
    for(int i=1;i<=grid.nX();i++) for(int j=1;j<=grid.nY();j++) v=max(v,un(i,j));
    return v;
}

void exportToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const string& fileName){
    ofstream f(fileName);
    for(int i=1;i<=un.nX();i++){
        for(int j=1;j<=un.nY();j++) f<<setprecision(16)<<un(i,j)<<" ";
        f<<"\n";
    }
}

void ensureDirectoryExists(const string& p){
    if(!fs::exists(p)) fs::create_directories(p);
}