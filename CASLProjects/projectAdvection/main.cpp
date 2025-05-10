#include <cmath>
#include <iostream>
#include <fstream>

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"
#include "../../CASLCommonLibrary/DPMatrix2D.h"

#include "projectAdvection_lib/CaslHamiltonianAdvection2D.h"

using namespace std;

int main()
{

  const double xMin = 0, xMax = 2 * M_PI;
  const int nX = 80;
  const double yMin = 0, yMax = 2 * M_PI;
  const int nY = 80;
  CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);
  double dx = grid.dx();
  double dy = grid.dy();

  double tInitial = 0, tFinal = 1, currentTime = tInitial;
  const int nOutput = 50;

  const int nPads = 3; // 3 pads all around since we will use at most ENO3 or WENO5.

  CaslArray2D<double> phiN(nX, nY, nPads);
  CaslArray2D<double> phiNp1(nX, nY, nPads);
  CaslArray2D<double> phiNp2(nX, nY, nPads);
  CaslArray2D<double> phiNp12(nX, nY, nPads);
  CaslArray2D<double> phiNp32(nX, nY, nPads);

  // velocity for advection HJ:
  CaslArray2D<double> c1(nX, nY), c2(nX, nY);
  for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      c1(i, j) = 0;
      c2(i, j) = 1;
    }

  for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      phiN(i, j) = cos(grid.x(i)) + sin(grid.y(j));
    }

  // Hamiltonian for linear advection:
  CaslHamiltonianAdvection2D hamiltonian(grid, c1, c2);
  CaslOptionPaddingWith boundaryConditions = withPeriodicCondition;
  CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;

  // time step for advection:
  double dt = 0.4 * std::min(std::min(std::min(dx / c1.maxAbs(), dy / c2.maxAbs()), dx), dy);

  CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, currentTime, firstDerivativeScheme);

  // March in time:
  while (currentTime < tFinal)
  {
    if (currentTime + dt > tFinal)
    {
      dt = tFinal - currentTime;
    }

    /*
     * TVD-RK3:
     */
    phiN.fillPaddingPoints(boundaryConditions);
    HJSolver.eulerStep(phiN, phiNp1);

    phiNp1.fillPaddingPoints(boundaryConditions);
    HJSolver.eulerStep(phiNp1, phiNp2);

    phiNp12 = .75 * phiN + .25 * phiNp2;

    phiNp12.fillPaddingPoints(boundaryConditions);
    HJSolver.eulerStep(phiNp12, phiNp32);

    phiNp1 = 1. / 3 * phiN + 2. / 3 * phiNp32;

    // Update for next time step:
    currentTime += dt;
    phiN = phiNp1;
    //cout << "currentTime = " << currentTime << endl;

    // Calculate the error:
    if (currentTime == tFinal)
    {
      double uExact;
      double maxError = 0.;
      for (int i = 1; i <= phiNp1.nX(); ++i)
      {
        for (int j = 1; j <= phiNp1.nY(); ++j)
        {
          uExact = cos(grid.x(i) - c1(i, j) * currentTime) + sin(grid.y(j) - c2(i, j) * currentTime);
          maxError = max(maxError, fabs(phiNp1(i, j) - uExact));
        }
      }
      cout << endl << "The grid size is " << nX << " by " << nY << endl;
      cout << "The maximum error is " << maxError << endl << endl;
    }
  }

  return EXIT_SUCCESS;
}