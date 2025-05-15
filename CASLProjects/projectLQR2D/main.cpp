// main_lqr2d.cpp
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslOptions.h"
#include "projectLQR2D_lib/CaslHamiltonianLQR2D.h"

using namespace std;

int main() {
    std::vector<int> gridSizes = {20, 40, 80, 160, 320};
    const double tFinal = 5.0;
    std::vector<double> exportTimes = {0, 1, 2, 3, 4, 5};

    for (int N : gridSizes) {
        cout << "\n=== Grid Size: " << N << " ===" << endl;

        CaslGrid2D grid(-2.0, 2.0, -2.0, 2.0, N, N);
        double dx = grid.dx();
        double dt = 0.5 * dx * dx;

        double currentTime = 0.0;
        const int nPads = 3;

        CaslArray2D<double> phi_n(N, N, nPads);
        CaslArray2D<double> phi_np1(N, N, nPads);
        CaslArray2D<double> phi_np2(N, N, nPads);
        CaslArray2D<double> phi_np12(N, N, nPads);
        CaslArray2D<double> phi_np32(N, N, nPads);

        // Terminal cost
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j) {
                double x = grid.x(i), y = grid.y(j);
                phi_n(i, j) = 0.5 * (x * x + y * y);
            }

        // System dynamics
        CaslArray2D<double> fx1(N, N), fx2(N, N);
        for (int i = 1; i <= N; ++i)
            for (int j = 1; j <= N; ++j) {
                fx1(i, j) = grid.y(j); // x2
                fx2(i, j) = 0.0;
            }

        CaslHamiltonianLQR2D hamiltonian(grid, fx1, fx2);
        CaslHamiltonJacobi2D<double> solver(grid, hamiltonian, dt, currentTime, WENO5);
        CaslCppToMATLAB2D matlabExporter;

        std::string folder = "./LQR2D_Output/LQR2D_" + to_string(N) + "/phi/";
        system(("mkdir -p " + folder).c_str());

        size_t nextExport = 0;
        while (currentTime < tFinal) {
            // Export if current time ≈ target
            if (nextExport < exportTimes.size() &&
                std::abs(currentTime - exportTimes[nextExport]) < 1e-2) {
                std::string filename = folder + "phi_t" + std::to_string(int(exportTimes[nextExport])) + ".dat";
                matlabExporter.exportDataToMatlab(grid, phi_n, filename);
                cout << "Exported: " << filename << " at time: " << currentTime << endl;
                nextExport++;
            }

            // Adjust final time step
            if (currentTime + dt > tFinal)
                dt = tFinal - currentTime;

            // TVD-RK3
            phi_n.fillPaddingPoints(withQuadraticExtrapolation);
            solver.eulerStep(phi_n, phi_np1);

            phi_np1.fillPaddingPoints(withQuadraticExtrapolation);
            solver.eulerStep(phi_np1, phi_np2);

            phi_np12 = 0.75 * phi_n + 0.25 * phi_np2;
            phi_np12.fillPaddingPoints(withQuadraticExtrapolation);
            solver.eulerStep(phi_np12, phi_np32);

            phi_np1 = (1.0 / 3.0) * phi_n + (2.0 / 3.0) * phi_np32;
            phi_n = phi_np1;

            currentTime += dt;
        }

        // Export final solution
        matlabExporter.exportDataToMatlab(grid, phi_n, folder + "phi_final.dat");
        cout << "Exported final solution at time: " << currentTime << endl;
    }

    return EXIT_SUCCESS;
}
