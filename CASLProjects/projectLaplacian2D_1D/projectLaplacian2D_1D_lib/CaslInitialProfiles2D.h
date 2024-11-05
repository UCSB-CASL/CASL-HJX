//
// Created by Faranak Rajabi on 2/17/24.
//

#ifndef CASL_INITIALPROFILES2D_H
#define CASL_INITIALPROFILES2D_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../../../CASLCommonLibrary/CaslArray2D.h"
#include "../../../CASLCommonLibrary/CaslGrid2D.h"

template<class T> class CaslInitialProfiles2D {
private:
    CaslGrid2D                          & _grid;                     // grid
    T                                   & _dt;                       // time step
    T                                   & _time;                     // time
    CaslArray2D<T>                        _initialVelocity;          // For advection type equations
    T                                     _diffusionCoefficient;     // Diffusion coefficient
    T                                     _L;                        // Rod length for heat type equations

public:
    // Constructor for CaslInitialProfiles2D class, defining
     explicit CaslInitialProfiles2D(CaslGrid2D & grid, T & dt, T & time, T diffusionCoefficient = T(), T  L = T(), CaslArray2D<T> initialVelocity = CaslArray2D<T> ());
    ~CaslInitialProfiles2D() = default;

    void heatSinIP1D                        (CaslArray2D<T> & phi_n    , int dimension          );
    void heatCosIP1D                        (CaslArray2D<T> & phi_n    , int dimension          );
    void heatCosIP1DExact                   (CaslArray2D<T> & phi_exact, int dimension, double t);
    void heatSinePolyIP1D                   (CaslArray2D<T> & phi_n    , int dimension, double t);
    void heatNonSmoothLineIP1D              (CaslArray2D<T> & phi_n    , int dimension          );
    void heatSinIP1DExact                   (CaslArray2D<T> & phi_exact, int dimension, double t);
    void heatSinePolyIP1DNonHomogeneousExact(CaslArray2D<T> & phi_exact, int dimension, double t);
    void heatNonSmoothLineIP1DExact         (CaslArray2D<T> & phi_n    , int dimension, double t);
    void advectionDiffusionCosIP1D          (CaslArray2D<T> & phi_n    , int dimension);
    void advectionDiffusionCosIP1DExact     (CaslArray2D<T> & phi_exact, int dimension, double t);
    void linearSDE                          (CaslArray2D<T> & phi_n_1  , CaslArray2D<double> & phi_n_2 , int dimension, double t);
    void linearSDEExact                     (CaslArray2D<T> & phi_exact_1, CaslArray2D<T> & phi_exact_2, int dimension, double t);
    void advectionSquaredIP1D               (CaslArray2D<T> & phi_n    , int dimension          );
    // needs modification
    void advectionSquaredIP1DExact          (CaslArray2D<T> & phi_exact, int dimension, double t);
    void advectionSinIP1D                   (CaslArray2D<T> & phi_n    , int dimension          );
    void advectionSinIP1DExact              (CaslArray2D<T> & phi_exact, int dimension, double t);
    void advectionDiffusionSinIP1D          (CaslArray2D<T> & phi_n    , int dimension          );
    void advectionDiffusionSinIP1DExact     (CaslArray2D<T> & phi_exact, int dimension, double t);
    void advectionDiffusionExpIP1D          (CaslArray2D<T> & phi_n    , int dimension);
    void advectionDiffusionExpIP1DExact     (CaslArray2D<T> & phi_exact, int dimension, double t);
    void linearSDEForward(CaslArray2D<T> & phi_n, int dimension);
    void linearSDEForwardExact(CaslArray2D<T> & phi_n, int dimension, double t);
};

#include "CaslInitialProfiles2D.cpp"

#endif //CASL_INITIALPROFILES2D_H

