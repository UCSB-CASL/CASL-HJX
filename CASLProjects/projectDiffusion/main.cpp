//
// Created by Faranak Rajabi at 01/17/24
//

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

#include "projectDiffusion_lib/CaslHamiltonianDiffusion2D.h"

using namespace std;

double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension);

int main()
{
  
  string projectName = "projectDiffusion";
  
  // Steps:
  // 1. Grid
  // 2. Padding
  // 3. Time set up
  // 4. Initial function definition
  // 5. CaslSecondOrder2D constructor
  //

  // Grid set up:
  int nX = 250, nY = 250;

  // For example 1 & 2:
  double xMin = -2, xMax = 2,
         yMin = -2, yMax = 2;

  CaslGrid2D grid(xMin, xMax, yMin, yMax, nX, nY);

  double initialT = 0.0, finalT = 1.0, current_time = initialT;
  double dt;

  dt = 0.02;

  // Phi arrays set up:
  const int nPadsX = 3;
  const int nPadsY = 3;

  CaslArray2D<double> phi_n(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phi_np1(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phi_np2(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phi_np12(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);
  CaslArray2D<double> phi_np32(nX, nY, nPadsX, nPadsX, nPadsY, nPadsY);

  // Initial function:
  // For example 1 & 2:
  // double L = 1;
  // For example 3:
  // double L = 2;
  // For example 1, 2, 3 & 5:
  double diffusionC = 0.25;
  //    double diffusionC = 1.0 / pow(M_PI, 2);
  // double diffusionC = 1;
  double L = 1.0;
  int dimension = 1.0; // X-direction

  // Inputs for the hamiltonian constructor, and initial functions if needed
  // For example 1, 2 & 3: Hamiltonian's all zeros, change in the following loop if needed
  // CaslArray2D<double> fx(nX, nY), fy(nX, nY);
  // For example 4: linear advection
  // For example 5: linear SDE a : c1
  CaslArray2D<double> c1(nX, nY);
  CaslArray2D<double> c2(nX, nY);
  for (int i = 1; i <= nX; i++)
  {
    for (int j = 1; j <= nY; j++)
    {
      c1(i, j) = 0.0;
      c2(i, j) = 0.0;
    }
  }

for (int i = 1; i <= nX; ++i)
    for (int j = 1; j <= nY; ++j)
    {
      if (pow(grid.x(i), 2) + pow(grid.y(j), 2) < 1) phi_n(i, j) = 1;
      else phi_n(i, j) = 0; 
    }
  // Export phi_0 to MATLAB
  CaslCppToMATLAB2D cppToMatlab;
  cppToMatlab.exportDataToMatlab(grid, phi_n, "../__Output/" + projectName + "/phi_0.dat");

  CaslOptionPaddingWith boundaryCondition = withConstantExtrapolation;
  CaslOptionNumericalFirstDerivative firstDerivativeScheme = WENO5;
  CaslOptionSecondOrderTermDirection secondOrderTermDirection = XY;
  CaslOptionNumericalSecondDerivative secondDerivativeScheme = BackwardTimeCentralSpacing;

  CaslHamiltonianDiffusion2D hamiltonianForDiffusion(grid, c1, c2);
  CaslHamiltonJacobi2D<double> HJSolver(grid, hamiltonianForDiffusion, dt, current_time, firstDerivativeScheme);
  CaslSecondOrderDerivative2D<double> LaplacianSolver(grid, dt, current_time, diffusionC, secondOrderTermDirection, secondDerivativeScheme, boundaryCondition, HJSolver, heatSourceFunction);

  int timestepNumber = 0;

  // Backward in time simulation:
  while (current_time < finalT)
  {
    // Compute dt using findCFL function:
    // dt = (CFL_number * LaplacianSolver.findCFLHamiltonianDiffusion(phi_n, current_time));
    // dt = HJSolver.findCFL(phi_n, current_time);
    if (current_time + dt > finalT)
    {
      dt = finalT - current_time;
      current_time = finalT;
    }

    // TVD-RK3:
    phi_n.fillPaddingPoints(boundaryCondition);
    HJSolver.eulerStep(phi_n, phi_np1);

    phi_np1.fillPaddingPoints(boundaryCondition);
    HJSolver.eulerStep(phi_np1, phi_np2);

    phi_np12 = .75 * phi_n + .25 * phi_np2;
    phi_np12.fillPaddingPoints(boundaryCondition);
    HJSolver.eulerStep(phi_np12, phi_np32);

    phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;

    LaplacianSolver.backwardTimeCentralSpacing(phi_n, phi_np1, 1e-30, 1e6);
    LaplacianSolver.backwardTimeCentralSpacing(phi_np1, phi_np2, 1e-30, 1e6);

    phi_np12 = .75 * phi_n + .25 * phi_np2;

    LaplacianSolver.backwardTimeCentralSpacing(phi_np12, phi_np32, 1e-30, 1e6);

    phi_np1 = 1. / 3 * phi_n + 2. / 3 * phi_np32;

    // Update for next time step:
    current_time += dt;
    phi_n = phi_np1;
    timestepNumber++;

    cout << to_string(current_time) << endl;
    cppToMatlab.exportDataToMatlab(grid, phi_n, "../__Output/" + projectName + "/phi_" + to_string(timestepNumber) + ".dat");
  }

  return EXIT_SUCCESS;
}

double heatSourceFunction(CaslGrid2D &grid, int i, int j, double time, int dimension)
{
  // For example 1, 3, 4 & 5: (heat1D_x), homogenous
  return 0.0;

  // For example 2: (heat1D_x_nonhomogeneous)
  // double x = (dimension == 1) ? grid.x(i) : grid.y(j);
  // return x;
}
