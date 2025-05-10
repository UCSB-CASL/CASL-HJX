// Created by Andrew Wang    2025-01-14

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

#include "projectBurgers_lib/CaslHamiltonianBurgers2D.h"

double homogeneousFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);

int main()
{
  double xMin = -M_PI_2, xMax = M_PI_2;
  double yMin = -M_PI_2, yMax = M_PI_2;
  int nX = 250, nY = 250;

  CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

  double tInitial = 0, tFinal = 4, tCurrent = tInitial;
  double dt;
  double nOutput = 200;

  const int nPadsX = 3;
  const int nPadsY = 3;

  CaslArray2D<double> un(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> unp1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> unp2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> unp12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> unp32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

  double dx = grid.dx(), dy = grid.dy();
  dt = 0.4 * std::min(std::min(std::min(dx / un.maxAbs(), dy / un.maxAbs()), dx), dy);

  for (int i = 1; i <= nX; i++)
    for (int j = 1; j <= nY; j++)
    {
      un(i, j) = exp(-1 * (grid.x(i) * grid.x(i) + grid.y(j) * grid.y(j)));
    }

  double diffusionC = 0.0;

  CaslOptionPaddingWith boundaryCondition = withConstantExtrapolation;
  CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
  CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
  CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

  CaslHamiltonianBurgers2D hamiltonian(grid, un, un);
  CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonian, dt, tCurrent, firstDerivativeScheme);
  CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, tCurrent, diffusionC, secondOrderTermDirection, secondDerivativeScheme, boundaryCondition, HJSolver, homogeneousFunction);

  string projectName = "projectBurgers";

  CaslCppToMATLAB2D cppToMatlab;

  string outputFolder = "../__Output/" + projectName + "/";
  string fileNameU = outputFolder + "un_";
  fileNameU += to_string(0) + ".dat";
  cout << "Exporting to MATLAB at time t = " << tCurrent << endl;
  cppToMatlab.exportDataToMatlab(grid, un, fileNameU);

  double dtGlobal = dt;
  bool isOutput = false;
  int nextFileIdx = 1;
  double saveIncrement = tFinal / nOutput;
  double nextOutputTime = saveIncrement;

  while (tCurrent < tFinal)
  {
    if (tCurrent + dt >= nextOutputTime)
    {
      dt = nextOutputTime - tCurrent;
      isOutput = true;
    }
    if (tCurrent + dt > tFinal)
    {
      dt = tFinal - tCurrent;
      isOutput = true;
    }

    LaplacianSolver.backwardTimeCentralSpacing(un, unp1);
    LaplacianSolver.backwardTimeCentralSpacing(unp1, unp2);
    unp12 = .75 * un + .25 * unp2;
    LaplacianSolver.backwardTimeCentralSpacing(unp12, unp32);
    unp1 = 1. / 3 * un + 2. / 3 * unp32;
    tCurrent += dt;
    un = unp1;
    cout << "t = " << tCurrent << endl;

    if (isOutput)
    {
      cout << "Exporting to MATLAB at time t = " << tCurrent << endl;
      fileNameU = outputFolder + "un_";
      fileNameU += to_string(nextFileIdx) + ".dat";
      cppToMatlab.exportDataToMatlab(grid, un, fileNameU);

      ++nextFileIdx;
      isOutput = false;
      dt = dtGlobal;
      nextOutputTime += saveIncrement;
    }
  }

  return EXIT_SUCCESS;
}

double homogeneousFunction(CaslGrid2D &grid, int i, int j, double time, int dimension)
{
  return 0.0;
}