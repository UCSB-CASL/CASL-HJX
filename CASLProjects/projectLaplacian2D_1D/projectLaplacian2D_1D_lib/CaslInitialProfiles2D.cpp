//
// Created by Faranak Rajabi on 2/17/24.
//

#ifndef CASL_INITIAL_PROFILES_2D_CPP
#define CASL_INITIAL_PROFILES_2D_CPP

#include "CaslInitialProfiles2D.h"

// Constructor for Advection-Diffusion type equations
template<class T>
CaslInitialProfiles2D<T>::CaslInitialProfiles2D(CaslGrid2D & grid, T & dt, T & time, T diffusionCoefficient, T L, CaslArray2D<T> initialVelocity) :
        _grid(grid),
        _dt(dt),
        _time(time),
        _diffusionCoefficient(diffusionCoefficient),
        _L(L),
        _initialVelocity(initialVelocity)
{}

template<class T>
void CaslInitialProfiles2D<T>::heatSinIP1D(CaslArray2D<T> & phi_n, int dimension) {
    // From: http://dma.ing.uniroma1.it/users/lsa_adn/MATERIALE/FDheat.pdf (page 12)
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatSinIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T phi_n_ij = (dimension == 1) ? sin(pi * _grid.x(i) / _L) : sin(pi * _grid.y(j) / _L);
            phi_n(i, j) = phi_n_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatSinIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    // From: http://dma.ing.uniroma1.it/users/lsa_adn/MATERIALE/FDheat.pdf (page 12)
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    // Compute phi_0 using heatSinIP1D
    CaslArray2D<T> phi_0(nX, nY);
    heatSinIP1D(phi_0, dimension);

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            phi_exact(i, j) = phi_0(i, j) * exp(- (_diffusionCoefficient * pi * pi * t) / (_L * _L));
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatCosIP1D(CaslArray2D<T> & phi_n, int dimension) {
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T phi_n_ij = (dimension == 1) ? cos(pi * (_grid.x(i) - 0.5)) : cos(pi * (_grid.y(j) - 0.5));
            phi_n(i, j) = phi_n_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatSinePolyIP1D(CaslArray2D<T> & phi_n, int dimension, double t) {
    // From: https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/07%3A_Green's_Functions/7.03%3A_The_Nonhomogeneous_Heat_Equation
    // Problem description:
    //    ut - uxx = x, 0 <= x <= 1.
    //    u(0, t) = 2, u(L, t) = t, t > 0
    //    u(x, 0) = 3 * sin(2 * π * x) + 2 * (1 - x),
    // Exact Solution:
    //    u(x, t) = 3 * sin(2𝜋𝑥) * exp(−4𝜋^2𝑡) + 2 + (𝑡 − 2) * 𝑥.
    // Note: This function computes u(x, 0), take L = 1.

    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            // Special consideration for x = 0, and x = L
            if (i == 1) {
                phi_n(i, j) = 2;
            }

            else if (i == nX) {
                phi_n(i, j) = t;
            }

            else {
                T sin_term_ij = (dimension == 1) ? sin(2.0 * pi * _grid.x(i)) : sin(2.0 * pi * _grid.y(j));
                T poly_term_ij = (dimension == 1) ? (1.0 - _grid.x(i)) : (1.0 - _grid.y(j));
                phi_n(i, j) = 3.0 * sin_term_ij + 2.0 * poly_term_ij;
            }
        }
    }

}

template<class T>
void CaslInitialProfiles2D<T>::heatSinePolyIP1DNonHomogeneousExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    // From: https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/07%3A_Green's_Functions/7.03%3A_The_Nonhomogeneous_Heat_Equation
    // Exact Solution for heatSinePolyIP1D(CaslArray2D<T> & phi_n, int dimension, double t)
    // u(x, t) = 3 * sin(2𝜋𝑥) * exp(−4𝜋^2𝑡) + 2 + (𝑡 − 2) * 𝑥.
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T sin_term_ij = (dimension == 1) ? sin(2.0 * pi * _grid.x(i)) : sin(2.0 * pi * _grid.y(j));
            T poly_term_ij = (dimension == 1) ? _grid.x(i) : _grid.y(j);
            phi_exact(i, j) = 3.0 * sin_term_ij * exp(-4.0 * pi * pi * t) + 2.0 + (t - 2.0) * poly_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatNonSmoothLineIP1D(CaslArray2D<T> & phi_n, int dimension) {
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T line_term_ij = (dimension == 1) ? _grid.x(i):_grid.y(j);
            if (line_term_ij < 1.0 && line_term_ij > 0.0) {
                phi_n(i, j) = line_term_ij;
            }
            else {
                phi_n(i, j) = 2.0 - line_term_ij;
            }
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatNonSmoothLineIP1DExact(CaslArray2D<T> & phi_n, int dimension, double t) {
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

        for (int i = 1; i <= nX; i++) {
            for (int j = 1; j <= nY; j++) {
                double sum = 0.0;
                for (int n = 1; n <= 100000; n++) {
                        sum = sum + (sin((2.0 * n - 1) * pi / 4.0) / (pow((2.0 * n - 1.0), 2))) * exp(- (pow((2.0 * n - 1.0), 2) * pi * pi * t / 16.0)) * sin((2.0 * n - 1 / 4.0) * pi * _grid.x(i));
                        phi_n(i, j) =  (32.0 / pow(pi, 2)) * sum;
                     }
                }
            }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionCosIP1D(CaslArray2D<T> & phi_n, int dimension) {
    // Note: for this function, domain should be in Pi scale!
    int nX = _grid.nX(), nY = _grid.nY();

    double k = 1.0;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T cos_term_ij = (dimension == 1) ? cos(k * _grid.x(i)) : cos(k * _grid.y(j));
            phi_n(i, j) = cos_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionCosIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    // Periodic BC
    // ut = -c * u_x + v * u_xx
    // v: diffusion coefficient
    // c: velocity
    // c = 1, v = 0.1

    // Test case 4:
    // ut = -c * u_x + v * u_xx
    // Periodic BC
    // v: diffusion coefficient
    // c: velocity
    // c = 1, v = 0.1
    // u(x, 0) = cos(k*pi*x);

    int nX = _grid.nX(), nY = _grid.nY();

    double k = 1.0;
    double c;
    double v = _diffusionCoefficient;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            c = _initialVelocity(i, j);
            T cos_term_ij = (dimension == 1) ? cos(k * (_grid.x(i) - c * t)) : cos(k * (_grid.y(j) - c * t));
            phi_exact(i, j) = exp(-v * k * k * t) * cos_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::linearSDE(CaslArray2D<T> & phi_n_1, CaslArray2D<double> & phi_n_2, int dimension, double t) {
    // Ian Mitchel Linear SDE:
    // phi(T) = x;
    // phi_n_1: expected value function as a terminal cost E[x(T)]
    // phi_n_2: Variance function at the final point Var[x(T)]

    // Analytic results.
    // mean_t0 = stateX;
    // var_t0 = zeros(size(mean_t0));
    // mean_tf = exp(a * tf) * mean_t0;
    // var_tf = 0.5 / a * ((b^2 + 2 * a * var_t0) * exp(2 * a * tf) - b^2);

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T poly_term_ij = (dimension == 1) ?  _grid.x(i): _grid.y(j);
            phi_n_1(i, j) = poly_term_ij;
            phi_n_2(i, j) = pow(_grid.x(i), 2);
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::linearSDEExact(CaslArray2D<T> & phi_exact_1, CaslArray2D<T> & phi_exact_2, int dimension, double t) {
    // E[x]   = exp(a * t) * E(x(T)]
    // Var[x] = 0.5/a * [(b^2 + 2a * Var[x(T)])exp(2at) - b^2]
    // Assuming E(x(T)) = x, Var[x] = 0 for T = tf;
    int nX = _grid.nX(), nY = _grid.nY();

    CaslArray2D<double> phi_0_1(phi_exact_1.nX(), phi_exact_1.nY()); // E[x(T)]
    CaslArray2D<double> phi_0_2(phi_exact_1.nX(), phi_exact_1.nY()); // Var[x(T)]
    linearSDE(phi_0_1, phi_0_2, 1, t);

    double a;
    auto b_squared = _diffusionCoefficient;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            a = _initialVelocity(i, j);
            phi_exact_1(i, j) = phi_0_1(i, j) * exp(a * t);
            phi_exact_2(i, j) = (0.5 / a) * ((b_squared + 2.0 * a * phi_0_2(i, j)) * exp(2.0 * a * t) - b_squared);
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::heatCosIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;

    CaslArray2D<double> phi_0(phi_exact.nX(), phi_exact.nY());
    heatCosIP1D(phi_0, 1);
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            phi_exact(i, j) = phi_0(i, j) * exp(-t);
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionSquaredIP1D(CaslArray2D<T> & phi_n, int dimension) {
    int nX = _grid.nX(), nY = _grid.nY();

    CaslArray2D<double> phi_0(phi_n.nX(), phi_n.nY());
    heatCosIP1D(phi_0, 1);
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T x_term_ij = (dimension == 1) ? _grid.x(i) : _grid.y(j);
            if (x_term_ij > 1 & x_term_ij < 2) {
                phi_n(i, j) = 1.0;
            }
            else {
                phi_n(i, j) = 0.0;
            }
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionSquaredIP1DExact(CaslArray2D<T> &phi_exact, int dimension, double t) {

    int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T c = _initialVelocity(i, j);
            T x;
            if (dimension == 1) {
                x = _grid.x(i);
            } else {
                x = _grid.y(j);
            }

            if (max(1.0, x - c*t) < x && x < min(2.0, x + c*t)) {
                phi_exact(i,j) = 1.0;
            } else {
                phi_exact(i,j) = 0.0;
            }

        }
    }

}

template<class T>
void CaslInitialProfiles2D<T>::advectionSinIP1D(CaslArray2D<T> & phi_n, int dimension) {
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;
    T k = 2.0;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T sin_term_ij = (dimension == 1) ? sin(k * pi * _grid.x(i)) : sin(k * pi * _grid.y(j));
            phi_n(i, j) = sin_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionSinIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;
    T k = 2.0;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T c1 = _initialVelocity(i, j);
            T sin_term_ij = (dimension == 1) ? sin(k * pi * (_grid.x(i) - c1 * t)) : sin(k * pi * (_grid.y(j) - c1 * t));
            phi_exact(i, j) = sin_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionSinIP1D(CaslArray2D<T> & phi_n, int dimension) {
    // u(x, 0) = sin(pi * x)
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }

    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;
    T k = 1.0;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T c1 = _initialVelocity(i, j);
            T sin_term_ij = (dimension == 1) ? sin(2.0 * pi * k * _grid.x(i)) : sin(2.0 * pi * k * _grid.y(j));
            phi_n(i, j) = sin_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionSinIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t){
    // u(x, t) = A * sin(2 * pi * k * (x - c1*t)) * exp(- 4 * pi^2 * k^2 * D * t)
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }
    int nX = _grid.nX(), nY = _grid.nY();
    auto pi = M_PI;
    T k = 2.0;

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T c1 = _initialVelocity(i, j);
            // needs modification
            T sin_term_ij = (dimension == 1) ? sin(2.0 * pi * k * (_grid.x(i) - c1 * t)) : sin(2.0 * pi * k * (_grid.y(j) - c1 * t));
            T exp_term_ij = exp(- 4.0 * pi * pi * k * k * _diffusionCoefficient * t);
            phi_exact(i, j) = sin_term_ij * exp_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionExpIP1D(CaslArray2D<T> & phi_n, int dimension) {
    // From: https://www.ajbasweb.com/old/ajbas/2011/june-2011/1536-1543.pdf
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }
    int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T exp_term_ij = (dimension == 1) ? exp(-pow((_grid.x(i) + 0.5), 2) / 0.00125) : exp(-pow((_grid.y(j) + 0.5), 2) / 0.00125);
            phi_n(i, j) = exp_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::advectionDiffusionExpIP1DExact(CaslArray2D<T> & phi_exact, int dimension, double t) {
    // From: https://www.ajbasweb.com/old/ajbas/2011/june-2011/1536-1543.pdf
    if (dimension < 1 || dimension > 2) {
        std::cerr << "In CaslInitialProfiles2D::heatCosIP1D, invalid dimension. EXITING." << std::endl;
        exit(1);
    }
    int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            T exp_term_ij = (dimension == 1) ? exp(-pow((_grid.x(i) + 0.5 - t), 2) / (0.00125 + 0.04 * t)) : exp(-pow((_grid.y(j) + 0.5 - t), 2) / (0.00125 + 0.04 * t)) ;
            T sqrtT_term_ij =  0.025 / sqrt(0.000625 + 0.02 * t);
            phi_exact(i, j) = exp_term_ij * sqrtT_term_ij;
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::linearSDEForward(CaslArray2D<T> & phi_n, int dimension) {
    // E[x]   = exp(a * t) * E(x(T)]
    // Var[x] = 0.5/a * [(b^2 + 2a * Var[x(T)])exp(2at) - b^2]
    // Assuming E(x(T)) = x, Var[x] = 0 for T = tf;
    int nX = _grid.nX(), nY = _grid.nY();

    T tf = 0.25;
    double a;
    auto b_squared = _diffusionCoefficient;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            a = _initialVelocity(i, j);
            phi_n(i, j) = _grid.x(i) * exp(a * tf);
        }
    }
}

template<class T>
void CaslInitialProfiles2D<T>::linearSDEForwardExact(CaslArray2D<T> & phi_n, int dimension, double t) {
    // E[x]   = exp(a * t) * E(x(T)]
    // Var[x] = 0.5/a * [(b^2 + 2a * Var[x(T)])exp(2at) - b^2]
    // Assuming E(x(T)) = x, Var[x] = 0 for T = tf;
    int nX = _grid.nX(), nY = _grid.nY();

    T tf = 0.25;
    double a;
    auto b_squared = _diffusionCoefficient;
    for (int i = 1; i <= nX; i++) {
        for (int j = 1; j <= nY; j++) {
            a = _initialVelocity(i, j);
            phi_n(i, j) = _grid.x(i) * exp(a * (tf - t));
        }
    }
}

#endif //CASL_INITIAL_PROFILES_2D_CPP
