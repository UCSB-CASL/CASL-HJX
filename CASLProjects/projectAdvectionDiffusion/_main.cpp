#include <cmath>
#include <iostream>
#include <fstream>

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonian2D.h"
#include "../../CASLCommonLibrary/CaslSecondOrderDerivative2D_1D.h"

#include "projectDiffusion_lib/CaslHamiltonianDiffusion2D.h"


using namespace std;

double sourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);

int main()
{
  int nX = 250;
  int nY = nX;

  double xMin = -M_PI_2, xMax = M_PI_2, yMin = xMin, yMax = xMax;

  CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

  double dx = grid.dx();
  double dy = grid.dy();

  double tInitial = 0.0, tFinal = 1.0, currentTime = tInitial;
  double dt = 0.01;

  const int nPadsX = 3;
  const int nPadsY = 3;

  CaslArray2D<double> phiN(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phiNp1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phiNp2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phiNp12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phiNp32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

  for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      if (pow(grid.x(i) - 1, 2) + pow(grid.y(j), 2) < 0.25 * 0.25)
        phiN(i, j) = 1.0;
      else
        phiN(i, j) = 0;
    }

  CaslArray2D<double> c1(nX, nY), c2(nX, nY);
  for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      c1(i, j) = 0;
      c2(i, j) = 0;
    }

  double diffusionCoefficient = 0.5;
  CaslOptionPaddingWith boundaryCondition = withPeriodicCondition;
  CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
  CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
  CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

  CaslHamiltonianDiffusion2D hamiltonianForDiffusion(grid, c1, c2);
  CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianForDiffusion, dt, currentTime, firstDerivativeScheme);
  CaslSecondOrderDerivative2D<double> diffusionSolver(grid, dt, currentTime, diffusionCoefficient, secondOrderTermDirection, secondDerivativeScheme, boundaryCondition, HJSolver, sourceFunction);

  

  CaslCppToMATLAB2D cppToMatlab;
  cppToMatlab.exportDataToMatlab(grid, phiN, "phi_0.dat");


  int exportCount = 1;

  while (currentTime < tFinal)
  {
    if (currentTime + dt > tFinal)
    {
      dt = tFinal - currentTime;
      currentTime = tFinal;
    }

    for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      c1(i, j) = -cos(grid.x(i)) * sin(grid.y(j));
      c2(i, j) = sin(grid.x(i)) * cos(grid.y(j));
    }

    cout << "currentTime = " << currentTime << endl;

    diffusionSolver.backwardTimeCentralSpacing(phiN, phiNp1);
    diffusionSolver.backwardTimeCentralSpacing(phiNp1, phiNp2);

    phiNp12 = .75 * phiN + .25 * phiNp2;

    diffusionSolver.backwardTimeCentralSpacing(phiNp12, phiNp32);

    phiNp1 = 1. / 3 * phiN + 2. / 3 * phiNp32;

    cppToMatlab.exportDataToMatlab(grid, phiNp1, "phi_" + to_string(exportCount) + ".dat");

    ++exportCount;
    currentTime += dt;
    phiN = phiNp1;
  }
}

double sourceFunction(CaslGrid2D & grid, int i, int j, double time, int dimension) {
     return 0.0;
}