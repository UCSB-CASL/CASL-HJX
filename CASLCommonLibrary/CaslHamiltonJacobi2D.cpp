//
// Created by Frederic Gibou on 1/5/23.
//

#ifndef CASL_HAMILTON_JACOBI2D_CPP
#define CASL_HAMILTON_JACOBI2D_CPP

#include <cmath>
#include "CaslHamiltonJacobi2D.h"

template<class T>
CaslHamiltonJacobi2D<T>::CaslHamiltonJacobi2D(CaslGrid2D & grid, CaslHamiltonian2D & hamiltonian, T & dt, T & time, CaslOptionNumericalFirstDerivative firstDerivativeScheme) :
        _grid(grid), _hamiltonian(hamiltonian), _time(time), _dt(dt), _firstDerivativeScheme(firstDerivativeScheme) {}

template<class T>
CaslHamiltonJacobi2D<T>::~CaslHamiltonJacobi2D() = default;

template<class T>
CaslArray2D<T> CaslHamiltonJacobi2D<T>::computeNumericalHamiltonian(const CaslArray2D<double>& un) {
    int nX = _grid.nX();
    int nY = _grid.nY();

    CaslArray2D<double> DxMinus(nX, nY), DxPlus(nX, nY);
    CaslArray2D<double> DyMinus(nX, nY), DyPlus(nX, nY);

    computeDxMinusDxPlus(un, DxMinus, DxPlus);
    computeDyMinusDyPlus(un, DyMinus, DyPlus);

    CaslArray2D<T> numericalHamiltonian(un.nX(), un.nY());

    double alphaX, alphaY;
    double ux_min = std::min(DxMinus(1, 1), DxPlus(1, 1));
    double ux_max = std::max(DxMinus(1, 1), DxPlus(1, 1));
    double uy_min = std::min(DyMinus(1, 1), DyPlus(1, 1));
    double uy_max = std::max(DyMinus(1, 1), DyPlus(1, 1));

    // Get interval:
    for (int i = 1; i <= nX; ++i)
        for (int j = 1; j <= nY; ++j) {
            ux_min = std::min(ux_min, std::min(DxMinus(i, j), DxPlus(i, j)));
            ux_max = std::max(ux_max, std::max(DxMinus(i, j), DxPlus(i, j)));
            uy_min = std::min(uy_min, std::min(DyMinus(i, j), DyPlus(i, j)));
            uy_max = std::max(uy_max, std::max(DyMinus(i, j), DyPlus(i, j)));
        }

    //  Update the solution:
    double ux_ave, uy_ave, ux_min_local, ux_max_local, uy_min_local, uy_max_local;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            ux_ave = 0.5 * (DxMinus(i, j) + DxPlus(i, j));
            uy_ave = 0.5 * (DyMinus(i, j) + DyPlus(i, j));

            ux_min_local = min(DxMinus(i, j), DxPlus(i, j));
            ux_max_local = max(DxMinus(i, j), DxPlus(i, j));
            uy_min_local = min(DyMinus(i, j), DyPlus(i, j));
            uy_max_local = max(DyMinus(i, j), DyPlus(i, j));

            alphaX = _hamiltonian.maxAbsH1(ux_min_local, ux_max_local, uy_min, uy_max, i, j, _time);
            alphaY = _hamiltonian.maxAbsH2(ux_min, ux_max, uy_min_local, uy_max_local, i, j, _time);
            numericalHamiltonian(i, j) = _hamiltonian.H(ux_ave, uy_ave, i, j, _time)
                                   - alphaX * 0.5 * (DxPlus(i, j) - DxMinus(i, j))
                                   - alphaY * 0.5 * (DyPlus(i, j) - DyMinus(i, j));
        }
    }
    return numericalHamiltonian;
}

template<class T>
void CaslHamiltonJacobi2D<T>::eulerStep(const CaslArray2D<double>& un, CaslArray2D<double>& unp1) {
    int nX = _grid.nX();
    int nY = _grid.nY();
    CaslArray2D<T> numericalHamiltonian = computeNumericalHamiltonian(un);

    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            unp1(i, j) = un(i, j) - _dt * numericalHamiltonian(i, j); // + dt * (cos(grid.x(i) + sin(grid.y(j) ...)
        }
    }

}

template<class T>
double CaslHamiltonJacobi2D<T>::findCFL(const CaslArray2D<double>& un, const double time)
{
    int nX = _grid.nX();
    int nY = _grid.nY();

    CaslArray2D<double> DxMinus(nX, nY), DxPlus(nX, nY);
    CaslArray2D<double> DyMinus(nX, nY), DyPlus(nX, nY);

    computeDxMinusDxPlus(un, DxMinus, DxPlus);
    computeDyMinusDyPlus(un, DyMinus, DyPlus);

    double maxAbsH1 = _hamiltonian.maxAbsH1(DxMinus(1, 1), DxPlus(1, 1), DyMinus(1, 1), DyPlus(1, 1), 1, 1, time);
    double maxAbsH2 = _hamiltonian.maxAbsH2(DxMinus(1, 1), DxPlus(1, 1), DyMinus(1, 1), DyPlus(1, 1), 1, 1, time);
    for(int i=1;i<=nX;i++) for(int j=1;j<=nY;j++){
            double DxMinLocal = std::min(DxMinus(i, j), DxPlus(i, j));
            double DxMaxLocal = std::max(DxMinus(i, j), DxPlus(i, j));
            double DyMinLocal = std::min(DyMinus(i, j), DyPlus(i, j));
            double DyMaxLocal = std::max(DyMinus(i, j), DyPlus(i, j));
            maxAbsH1=max(maxAbsH1, _hamiltonian.maxAbsH1(DxMinLocal, DxMaxLocal, DyMinLocal, DyMaxLocal, i, j, time));
            maxAbsH2=max(maxAbsH2, _hamiltonian.maxAbsH2(DxMinLocal, DxMaxLocal, DyMinLocal, DyMaxLocal, i, j, time));}

    // return  1.0 / (maxAbsH1 / _grid.dx() + maxAbsH2 / _grid.dy());
    /* changed to the inverse shape to use in CaslSecondOrderDerivative2D class */
     return (maxAbsH1 / _grid.dx() + maxAbsH2 / _grid.dy());
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeOptimalControl(const CaslArray2D<double>& un, int K, double uMax, CaslArray2D<double>& uStar) {
    int nX = _grid.nX();
    int nY = _grid.nY();

    CaslArray2D<double> DxMinus(nX, nY), DxPlus(nX, nY);
    computeDxMinusDxPlus(un, DxMinus, DxPlus);

    //  Update the solution:
    double ux_ave;
    for (int i = 1; i <= nX; ++i) {
        for (int j = 1; j <= nY; ++j) {
            ux_ave = 0.5 * (DxMinus(i, j) + DxPlus(i, j));
            if (fabs(ux_ave) <= 2.0 * K * uMax) {
                uStar(i, j) = (-0.5 * ux_ave) / K;
            }
            if (fabs(ux_ave) > 2.0 * K * uMax) {
                if (ux_ave > 0) {
                    uStar(i, j) = -uMax;
                }
                if (ux_ave <= 0) {
                    uStar(i, j) = uMax;
                }
            }
        }
    }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDxMinusDxPlus(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus) {
    if (_firstDerivativeScheme == Upwind) { computeDxMinusDxPlusUpwind(un, DxMinus, DxPlus); return;}
    if (_firstDerivativeScheme == ENO2  ) { computeDxMinusDxPlusENO2  (un, DxMinus, DxPlus); return;}
    if (_firstDerivativeScheme == ENO3  ) { computeDxMinusDxPlusENO3  (un, DxMinus, DxPlus); return;}
    if (_firstDerivativeScheme == WENO5 ) { computeDxMinusDxPlusWENO5 (un, DxMinus, DxPlus); return;}

    std::cout << "In CaslHamiltonJacobi2D::ComputeDxMinus, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDxMinusDxPlusUpwind(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus) {
    T dx = _grid.dx(); int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            DxMinus(i, j) = ( un(i,j) - un(i-1,j) ) / dx;
            DxPlus (i, j) = ( un(i+1,j) - un(i,j) ) / dx;
        }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDxMinusDxPlusENO2(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus) {
    T dx = _grid.dx(); int nX = _grid.nX(), nY = _grid.nY();

    int orderENO = 2;

    std::vector< CaslArray2D<T> > DD(orderENO + 1, CaslArray2D<T>(nX, nY, orderENO, orderENO, orderENO, orderENO) );

    // compute the divided difference table:
    int iMin = DD[0].iMin(), iMax = DD[0].iMax();
    for (int i = iMin; i <= iMax; ++i) for (int j = 1; j <= nY; ++j) DD[0](i,j) = un(i,j);
    for (int level = 1; level <= orderENO; ++level)
        for (int i = iMin + level; i <= iMax-level; ++i)
            for (int j = 1; j <= nY; ++j)
                DD[level](i,j) = ( DD[level - 1](i+1,j) - DD[level - 1](i,j)) / level / dx;

    T dQ1dx, dQ2dx, c;

    // DxMinus:
    for (int i = 1; i <= nX; ++i) {
        int k = i - 1;
        for (int j = 1; j <= nY; ++j) {
            dQ1dx = DD[1](k, j);

            if (fabs(DD[2](k - 1, j)) < fabs(DD[2](k, j))) {
                c = DD[2](k - 1, j);
            } else {
                c = DD[2](k, j);
            }
            dQ2dx = c * (2 * (i - k) - 1) * dx;

            DxMinus(i, j) = dQ1dx + dQ2dx;
        }
    }
    // DxPlus:
    for (int i = 1; i <= nX; ++i) {
        int k = i;
        for (int j = 1; j <= nY; ++j) {
            dQ1dx = DD[1](k, j);
            if (fabs(DD[2](k - 1, j)) < fabs(DD[2](k, j))) {
                c = DD[2](k - 1, j);
            } else {
                c = DD[2](k, j);
            }
            dQ2dx = c * (2 * (i - k) - 1) * dx;

            DxPlus(i, j) = dQ1dx + dQ2dx;
        }
    }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDxMinusDxPlusENO3(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus) {
    T dx = _grid.dx(); int nX = _grid.nX(), nY = _grid.nY();

    int orderENO = 3;

    std::vector< CaslArray2D<T> > DD(orderENO + 1, CaslArray2D<T>(nX, nY, orderENO, orderENO, orderENO, orderENO) );

    // compute the divided difference table:
    int iMin = DD[0].iMin(), iMax = DD[0].iMax();
    for (int i = iMin; i <= iMax; ++i) for (int j = 1; j <= nY; ++j) DD[0](i,j) = un(i,j);
    for (int level = 1; level <= orderENO; ++level)
        for (int i = iMin + level; i <= iMax-level; ++i)
            for (int j = 1; j <= nY; ++j)
                DD[level](i,j) = ( DD[level - 1](i+1,j) - DD[level - 1](i,j)) / level / dx;

    T dQ1dx, dQ2dx, dQ3dx, c, cStar;
    long kStar;

    // DxMinus:
    for (int i = 1; i <= nX; ++i) {
        int k = i - 1;
        for (int j = 1; j <= nY; ++j) {
            dQ1dx = DD[1](k, j);

            if (fabs(DD[2](k - 1, j)) < fabs(DD[2](k, j))) {
                c = DD[2](k - 1, j); kStar = k - 1;
            } else {
                c = DD[2](k, j);     kStar = k;
            }
            dQ2dx = c * (2 * (i - k) - 1) * dx;

            if (fabs(DD[3](kStar - 1, j)) < fabs(DD[3](kStar, j))) cStar = DD[3](kStar - 1, j);
            else cStar = DD[3](kStar, j);
            dQ3dx = cStar * (3 * (i - kStar) * (i - kStar) - 6 * (i - kStar) + 2) * dx * dx;

            DxMinus(i, j) = dQ1dx + dQ2dx + dQ3dx;
        }
    }
    // DxPlus:
    for (int i = 1; i <= nX; ++i) {
        int k = i;
        for (int j = 1; j <= nY; ++j) {
            dQ1dx = DD[1](k, j);
            if (fabs(DD[2](k - 1, j)) < fabs(DD[2](k, j))) {
                c = DD[2](k - 1, j);
                kStar = k - 1;
            } else {
                c = DD[2](k, j);
                kStar = k;
            }
            dQ2dx = c * (2 * (i - k) - 1) * dx;

            if (fabs(DD[3](kStar - 1, j)) < fabs(DD[3](kStar, j))) cStar = DD[3](kStar - 1, j);
            else cStar = DD[3](kStar, j);
            dQ3dx = cStar * (3 * (i - kStar) * (i - kStar) - 6 * (i - kStar) + 2) * dx * dx;

            DxPlus(i, j) = dQ1dx + dQ2dx + dQ3dx;
        }
    }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDxMinusDxPlusWENO5(const CaslArray2D<double>& un, CaslArray2D<double>& DxMinus, CaslArray2D<double>& DxPlus) {
    T dx = _grid.dx(); int nX = _grid.nX(), nY = _grid.nY();

    double ux1, ux2, ux3;
    double omega1, omega2, omega3;
    double alpha1, alpha2, alpha3, sumAlphas;
    double d1, d2, d3, d4, d5, d1Square, d2Square, d3Square, d4Square, d5Square;
    double S1, S2, S3;
    double epsilon;

    // DxMinus
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            d1 = ( un(i-2, j) - un(i-3, j) ) / dx;
            d2 = ( un(i-1, j) - un(i-2, j) ) / dx;
            d3 = ( un(i-0, j) - un(i-1, j) ) / dx;
            d4 = ( un(i+1, j) - un(i-0, j) ) / dx;
            d5 = ( un(i+2, j) - un(i+1, j) ) / dx;

            S1 = 13./12 * pow(d1 - 2*d2 + d3, 2) + 0.25 * pow(  d1 - 4*d2 + 3*d3,2);
            S2 = 13./12 * pow(d2 - 2*d3 + d4, 2) + 0.25 * pow(  d2        -   d4,2);
            S3 = 13./12 * pow(d3 - 2*d4 + d5, 2) + 0.25 * pow(3*d3 - 4*d4 +   d5,2);

            d1Square = pow(d1, 2);
            d2Square = pow(d2, 2);
            d3Square = pow(d3, 2);
            d4Square = pow(d4, 2);
            d5Square = pow(d5, 2);
            epsilon = 1e-6 * std::max(d1Square, std::max(d2Square, std::max(d3Square, std::max(d4Square, d5Square)))) + 10e-99;

            alpha1 = 0.1 / pow( S1 + epsilon, 2 );
            alpha2 = 0.6 / pow( S2 + epsilon, 2 );
            alpha3 = 0.3 / pow( S3 + epsilon, 2 );

            sumAlphas = alpha1 + alpha2 + alpha3;
            omega1 = alpha1 / sumAlphas;
            omega2 = alpha2 / sumAlphas;
            omega3 = 1. - omega1 - omega2;

            ux1 = (  2.*d1 - 7.*d2 + 11.*d3 ) / 6.0;
            ux2 = ( -1.*d2 + 5.*d3 + 2. *d4 ) / 6.0;
            ux3 = (  2.*d3 + 5.*d4 - 1. *d5 ) / 6.0;

            DxMinus(i, j) = omega1 * ux1 + omega2 * ux2 + omega3 * ux3;
        }

    // DxPlus:
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            d1 = ( un(i+3, j) - un(i+2, j) ) / dx;
            d2 = ( un(i+2, j) - un(i+1, j) ) / dx;
            d3 = ( un(i+1, j) - un(i+0, j) ) / dx;
            d4 = ( un(i+0, j) - un(i-1, j) ) / dx;
            d5 = ( un(i-1, j) - un(i-2, j) ) / dx;

            S1 = 13./12 * pow(d1 - 2*d2 + d3, 2) + 0.25 * pow(  d1 - 4*d2 + 3*d3,2);
            S2 = 13./12 * pow(d2 - 2*d3 + d4, 2) + 0.25 * pow(  d2        -   d4,2);
            S3 = 13./12 * pow(d3 - 2*d4 + d5, 2) + 0.25 * pow(3*d3 - 4*d4 +   d5,2);

            d1Square = d1 * d1;
            d2Square = d2 * d2;
            d3Square = d3 * d3;
            d4Square = d4 * d4;
            d5Square = d5 * d5;
            epsilon = 1e-6 * std::max(d1Square, std::max(d2Square, std::max(d3Square, std::max(d4Square, d5Square)))) + 10e-99;

            alpha1 = 0.1 / pow( S1 + epsilon, 2 );
            alpha2 = 0.6 / pow( S2 + epsilon, 2 );
            alpha3 = 0.3 / pow( S3 + epsilon, 2 );

            sumAlphas = alpha1 + alpha2 + alpha3;
            omega1 = alpha1 / sumAlphas;
            omega2 = alpha2 / sumAlphas;
            omega3 = 1. - omega1 - omega2;

            ux1 = (  2.*d1 - 7.*d2 + 11.*d3 ) / 6.0;
            ux2 = ( -1.*d2 + 5.*d3 + 2. *d4 ) / 6.0;
            ux3 = (  2.*d3 + 5.*d4 - 1. *d5 ) / 6.0;

            DxPlus(i, j) = omega1 * ux1 + omega2 * ux2 + omega3 * ux3;
        }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDyMinusDyPlus(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus) {

    if (_firstDerivativeScheme == Upwind) { computeDyMinusDyPlusUpwind(un, DyMinus, DyPlus); return;}
    if (_firstDerivativeScheme == ENO2  ) { computeDyMinusDyPlusENO2  (un, DyMinus, DyPlus); return;}
    if (_firstDerivativeScheme == ENO3  ) { computeDyMinusDyPlusENO3  (un, DyMinus, DyPlus); return;}
    if (_firstDerivativeScheme == WENO5 ) { computeDyMinusDyPlusWENO5 (un, DyMinus, DyPlus); return;}

    std::cout << "In CaslHamiltonJacobi2D::computeDyMinusDyPlus, the solver does not exists. EXITING." << std::endl; exit(1);
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDyMinusDyPlusUpwind(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus) {
    T dy = _grid.dy(); int nX = _grid.nX(), nY = _grid.nY();

    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            DyMinus(i, j) = ( un(i,j) - un(i,j-1) ) / dy;
            DyPlus (i, j) = ( un(i,j+1) - un(i,j) ) / dy;
        }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDyMinusDyPlusENO2(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus) {
    T dy = _grid.dy(); int nX = _grid.nX(), nY = _grid.nY();

    int orderENO = 2;

    std::vector<CaslArray2D<T> > DD(orderENO + 1, CaslArray2D<T>(nX, nY, orderENO, orderENO, orderENO, orderENO));

    // compute divided difference tables:
    int jMin = DD[0].jMin(), jMax = DD[0].jMax();
    for (int i = 1; i <= nX; ++i) for (int j = jMin; j <= jMax; ++j) DD[0](i, j) = un(i, j);
    for (int level = 1; level <= orderENO; ++level)
        for (int i = 1; i <= nX; ++i)
            for (int j = jMin + level; j <= jMax - level; ++j)
                DD[level](i, j) = (DD[level - 1](i, j + 1) - DD[level - 1](i, j)) / level / dy;

    T dQ1dy, dQ2dy, c;

    // DyMinus:
    for (int j = 1; j <= nY; ++j) {
        int k = j - 1;
        for (int i = 1; i <= nX; ++i) {
            dQ1dy = DD[1](i, k);

            if (fabs(DD[2](i, k - 1)) < fabs(DD[2](i, k))) {
                c = DD[2](i, k - 1);
            } else {
                c = DD[2](i, k);
            }
            dQ2dy = c * (2 * (j - k) - 1) * dy;

            DyMinus(i, j) = dQ1dy + dQ2dy;
        }
    }

    // DyPlus:
    for (int j = 1; j <= nY; ++j) {
        int k = j;
        for (int i = 1; i <= nX; ++i) {
            dQ1dy = DD[1](i, k);

            if (fabs(DD[2](i, k - 1)) < fabs(DD[2](i, k))) {
                c = DD[2](i, k - 1);
            } else {
                c = DD[2](i, k);
            }
            dQ2dy = c * (2 * (j - k) - 1) * dy;

            DyPlus(i, j) = dQ1dy + dQ2dy;
        }
    }
}

template<class T>
void CaslHamiltonJacobi2D<T>::computeDyMinusDyPlusENO3(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus) {
    T dy = _grid.dy(); int nX = _grid.nX(), nY = _grid.nY();

    int order_ENO = 3;

    std::vector<CaslArray2D<T> > DD(order_ENO + 1, CaslArray2D<T>(nX, nY, order_ENO, order_ENO, order_ENO, order_ENO));

    // compute divided difference tables:
    int jMin = DD[0].jMin(), jMax = DD[0].jMax();
    for (int i = 1; i <= nX; ++i) for (int j = jMin; j <= jMax; ++j) DD[0](i, j) = un(i, j);
    for (int level = 1; level <= order_ENO; ++level)
        for (int i = 1; i <= nX; ++i)
            for (int j = jMin + level; j <= jMax - level; ++j)
                DD[level](i, j) = (DD[level - 1](i, j + 1) - DD[level - 1](i, j)) / level / dy;

    T dQ1dy, dQ2dy, dQ3dy, c, cStar;
    int kStar;

    // DyMinus:
    for (int j = 1; j <= nY; ++j) {
        int k = j - 1;
        for (int i = 1; i <= nX; ++i) {
            dQ1dy = DD[1](i, k);

            if (fabs(DD[2](i, k - 1)) < fabs(DD[2](i, k))) {
                c = DD[2](i, k - 1);
                kStar = k - 1;
            } else {
                c = DD[2](i, k);
                kStar = k;
            }
            dQ2dy = c * (2 * (j - k) - 1) * dy;

            if (fabs(DD[3](i, kStar - 1)) < fabs(DD[3](i, kStar))) cStar = DD[3](i, kStar - 1);
            else cStar = DD[3](i, kStar);
            dQ3dy = cStar * (3 * (j - kStar) * (j - kStar) - 6 * (j - kStar) + 2) * dy * dy;

            DyMinus(i, j) = dQ1dy + dQ2dy + dQ3dy;
        }
    }

    // DyPlus:
    for (int j = 1; j <= nY; ++j) {
        int k = j;
        for (int i = 1; i <= nX; ++i) {
            dQ1dy = DD[1](i, k);

            if (fabs(DD[2](i, k - 1)) < fabs(DD[2](i, k))) {
                c = DD[2](i, k - 1);
                kStar = k - 1;
            } else {
                c = DD[2](i, k);
                kStar = k;
            }
            dQ2dy = c * (2 * (j - k) - 1) * dy;

            if (fabs(DD[3](i, kStar - 1)) < fabs(DD[3](i, kStar))) cStar = DD[3](i, kStar - 1);
            else cStar = DD[3](i, kStar);
            dQ3dy = cStar * (3 * (j - kStar) * (j - kStar) - 6 * (j - kStar) + 2) * dy * dy;

            DyPlus(i, j) = dQ1dy + dQ2dy + dQ3dy;
        }
    }
}


template<class T>
void CaslHamiltonJacobi2D<T>::computeDyMinusDyPlusWENO5(const CaslArray2D<double>& un, CaslArray2D<double>& DyMinus, CaslArray2D<double>& DyPlus) {
    T dy = _grid.dy(); int nX = _grid.nX(), nY = _grid.nY();

    double uy1, uy2, uy3;
    double omega1, omega2, omega3;
    double alpha1, alpha2, alpha3, sumAlphas;
    double d1, d2, d3, d4, d5, d1Square, d2Square, d3Square, d4Square, d5Square;
    double S1, S2, S3;
    double epsilon;

    // DyMinus:
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            d1 = ( un(i, j-2) - un(i, j-3) ) / dy;
            d2 = ( un(i, j-1) - un(i, j-2) ) / dy;
            d3 = ( un(i, j-0) - un(i, j-1) ) / dy;
            d4 = ( un(i, j+1) - un(i, j-0) ) / dy;
            d5 = ( un(i, j+2) - un(i, j+1) ) / dy;

            S1 = 13./12 * pow(d1 - 2*d2 + d3, 2) + 0.25 * pow(  d1 - 4*d2 + 3*d3,2);
            S2 = 13./12 * pow(d2 - 2*d3 + d4, 2) + 0.25 * pow(  d2        -   d4,2);
            S3 = 13./12 * pow(d3 - 2*d4 + d5, 2) + 0.25 * pow(3*d3 - 4*d4 +   d5,2);

            d1Square = pow(d1, 2);
            d2Square = pow(d2, 2);
            d3Square = pow(d3, 2);
            d4Square = pow(d4, 2);
            d5Square = pow(d5, 2);
            epsilon = 1e-6 * std::max(d1Square, std::max(d2Square, std::max(d3Square, std::max(d4Square, d5Square)))) + 10e-99;

            alpha1 = 0.1 / pow( S1 + epsilon, 2 );
            alpha2 = 0.6 / pow( S2 + epsilon, 2 );
            alpha3 = 0.3 / pow( S3 + epsilon, 2 );

            sumAlphas = alpha1 + alpha2 + alpha3;
            omega1 = alpha1 / sumAlphas;
            omega2 = alpha2 / sumAlphas;
            omega3 = 1. - omega1 - omega2;

            uy1 = ( 2. * d1 - 7. * d2 + 11.* d3 ) / 6.0;
            uy2 = (-1. * d2 + 5. * d3 + 2. * d4 ) / 6.0;
            uy3 = ( 2. * d3 + 5. * d4 - 1. * d5 ) / 6.0;

            DyMinus(i, j) = omega1 * uy1 + omega2 * uy2 + omega3 * uy3;
        }

    // DyPlus:
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            d1 = ( un(i, j+3) - un(i, j+2) ) / dy;
            d2 = ( un(i, j+2) - un(i, j+1) ) / dy;
            d3 = ( un(i, j+1) - un(i, j+0) ) / dy;
            d4 = ( un(i, j+0) - un(i, j-1) ) / dy;
            d5 = ( un(i, j-1) - un(i, j-2) ) / dy;

            S1 = 13./12 * pow(d1 - 2*d2 + d3, 2) + 0.25 * pow(  d1 - 4*d2 + 3*d3,2);
            S2 = 13./12 * pow(d2 - 2*d3 + d4, 2) + 0.25 * pow(  d2        -   d4,2);
            S3 = 13./12 * pow(d3 - 2*d4 + d5, 2) + 0.25 * pow(3*d3 - 4*d4 +   d5,2);

            d1Square = d1 * d1;
            d2Square = d2 * d2;
            d3Square = d3 * d3;
            d4Square = d4 * d4;
            d5Square = d5 * d5;
            epsilon = 1e-6 * std::max(d1Square, std::max(d2Square, std::max(d3Square, std::max(d4Square, d5Square)))) + 10e-99;

            alpha1 = 0.1 / pow( S1 + epsilon, 2 );
            alpha2 = 0.6 / pow( S2 + epsilon, 2 );
            alpha3 = 0.3 / pow( S3 + epsilon, 2 );

            sumAlphas = alpha1 + alpha2 + alpha3;
            omega1 = alpha1 / sumAlphas;
            omega2 = alpha2 / sumAlphas;
            omega3 = 1. - omega1 - omega2;

            uy1 = ( 2. * d1 - 7. * d2 + 11.* d3 ) / 6.0;
            uy2 = (-1. * d2 + 5. * d3 + 2. * d4 ) / 6.0;
            uy3 = ( 2. * d3 + 5. * d4 - 1. * d5 ) / 6.0;

            DyPlus(i, j) = omega1 * uy1 + omega2 * uy2 + omega3 * uy3;
        }
}

#endif // CASL_HAMILTON_JACOBI2D_CPP