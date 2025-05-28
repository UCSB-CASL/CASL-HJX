//
// Created by Faranak Rajabi on 5/17/25.
//

// Provides a high-precision solver for the LQR Riccati equation
// to use for validation and hybrid methods

#ifndef RICCATI_SOLVER_H
#define RICCATI_SOLVER_H

#include <vector>
#include <algorithm>
#include <cmath>
#include "CaslGrid2D.h"
#include "CaslArray2D.h"

class RiccatiSolver {
private:
    // System matrices
    double _A[2][2];
    double _B[2];
    double _Q[2][2];
    double _R;

    // Riccati solution storage
    std::vector<double> _times;
    std::vector<double> _P11;
    std::vector<double> _P12;
    std::vector<double> _P22;

public:
    // Constructor sets up system matrices for standard double integrator
    RiccatiSolver() {
        // A = [0 1; 0 0]
        _A[0][0] = 0.0; _A[0][1] = 1.0;
        _A[1][0] = 0.0; _A[1][1] = 0.0;

        // B = [0; 1]
        _B[0] = 0.0; _B[1] = 1.0;

        // Q = I
        _Q[0][0] = 1.0; _Q[0][1] = 0.0;
        _Q[1][0] = 0.0; _Q[1][1] = 1.0;

        // R = 1
        _R = 1.0;
    }

    // Constructor with custom system matrices
    RiccatiSolver(double A[2][2], double B[2], double Q[2][2], double R) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                _A[i][j] = A[i][j];
                _Q[i][j] = Q[i][j];
            }
            _B[i] = B[i];
        }
        _R = R;
    }

    // Solve Riccati equation backwards from terminal time
    void solve(double tFinal, const std::vector<double>& exportTimes) {
        // Clear previous solutions
        _times.clear();
        _P11.clear();
        _P12.clear();
        _P22.clear();

        // Sort times for Riccati integration (from tFinal down to 0)
        std::vector<double> tauVec;
        for (double t : exportTimes) {
            tauVec.push_back(tFinal - t);
        }
        std::sort(tauVec.begin(), tauVec.end());

        // Terminal condition P(tFinal) = I
        double P11 = 1.0, P12 = 0.0, P22 = 1.0;

        // Add terminal condition to solution
        _times.push_back(tFinal);
        _P11.push_back(P11);
        _P12.push_back(P12);
        _P22.push_back(P22);

        // Use a small time step for accuracy
        double tau = 0.0;
        double dt_riccati = 0.0001;

        // Backward integration using 4th order Runge-Kutta
        int saveIndex = 0;
        while (tau < tFinal) {
            // RK4 integration step
            // K1
            double k1_11 = -(_Q[0][0] - P12*P12/_R + 2*P11*_A[0][1]);
            double k1_12 = -(_Q[0][1] + P11*_A[1][1] + P12*_A[0][0] - P12*P22/_R);
            double k1_22 = -(_Q[1][1] - P22*P22/_R + 2*P22*_A[1][0]);

            // K2
            double P11_mid = P11 + 0.5*dt_riccati*k1_11;
            double P12_mid = P12 + 0.5*dt_riccati*k1_12;
            double P22_mid = P22 + 0.5*dt_riccati*k1_22;

            double k2_11 = -(_Q[0][0] - P12_mid*P12_mid/_R + 2*P11_mid*_A[0][1]);
            double k2_12 = -(_Q[0][1] + P11_mid*_A[1][1] + P12_mid*_A[0][0] - P12_mid*P22_mid/_R);
            double k2_22 = -(_Q[1][1] - P22_mid*P22_mid/_R + 2*P22_mid*_A[1][0]);

            // K3
            P11_mid = P11 + 0.5*dt_riccati*k2_11;
            P12_mid = P12 + 0.5*dt_riccati*k2_12;
            P22_mid = P22 + 0.5*dt_riccati*k2_22;

            double k3_11 = -(_Q[0][0] - P12_mid*P12_mid/_R + 2*P11_mid*_A[0][1]);
            double k3_12 = -(_Q[0][1] + P11_mid*_A[1][1] + P12_mid*_A[0][0] - P12_mid*P22_mid/_R);
            double k3_22 = -(_Q[1][1] - P22_mid*P22_mid/_R + 2*P22_mid*_A[1][0]);

            // K4
            P11_mid = P11 + dt_riccati*k3_11;
            P12_mid = P12 + dt_riccati*k3_12;
            P22_mid = P22 + dt_riccati*k3_22;

            double k4_11 = -(_Q[0][0] - P12_mid*P12_mid/_R + 2*P11_mid*_A[0][1]);
            double k4_12 = -(_Q[0][1] + P11_mid*_A[1][1] + P12_mid*_A[0][0] - P12_mid*P22_mid/_R);
            double k4_22 = -(_Q[1][1] - P22_mid*P22_mid/_R + 2*P22_mid*_A[1][0]);

            // Update with 4th order accuracy
            P11 += (dt_riccati/6) * (k1_11 + 2*k2_11 + 2*k3_11 + k4_11);
            P12 += (dt_riccati/6) * (k1_12 + 2*k2_12 + 2*k3_12 + k4_12);
            P22 += (dt_riccati/6) * (k1_22 + 2*k2_22 + 2*k3_22 + k4_22);

            tau += dt_riccati;

            // Save solution at requested times
            while (saveIndex < tauVec.size() && std::abs(tau - tauVec[saveIndex]) < dt_riccati/2) {
                _times.push_back(tFinal - tau);
                _P11.push_back(P11);
                _P12.push_back(P12);
                _P22.push_back(P22);
                saveIndex++;
            }
        }
    }

    // Interpolate P matrix at time t
    void getP(double t, double& P11, double& P12, double& P22) const {
        // Find time indices for interpolation
        size_t idx = 0;
        while (idx < _times.size() && _times[idx] < t) idx++;

        if (idx == 0) {
            // Before the first time point, use the first value
            P11 = _P11.front();
            P12 = _P12.front();
            P22 = _P22.front();
        } else if (idx >= _times.size()) {
            // After the last time point, use the last value
            P11 = _P11.back();
            P12 = _P12.back();
            P22 = _P22.back();
        } else {
            // Interpolate between two time points
            double t0 = _times[idx-1];
            double t1 = _times[idx];
            double alpha = (t - t0) / (t1 - t0);

            P11 = (1-alpha) * _P11[idx-1] + alpha * _P11[idx];
            P12 = (1-alpha) * _P12[idx-1] + alpha * _P12[idx];
            P22 = (1-alpha) * _P22[idx-1] + alpha * _P22[idx];
        }
    }

    // Compute value function from Riccati solution at time t
    void getValueFunction(const CaslGrid2D& grid, double t, CaslArray2D<double>& phi) const {
        // Interpolate P(t) from stored values
        double P11, P12, P22;
        getP(t, P11, P12, P22);

        // Set value function V(x,t) = 0.5 * x^T P(t) x
        int nX = grid.nX();
        int nY = grid.nY();

        for (int i = 1; i <= nX; i++) {
            for (int j = 1; j <= nY; j++) {
                double x = grid.x(i);
                double y = grid.y(j);

                phi(i, j) = 0.5 * (P11*x*x + 2*P12*x*y + P22*y*y);
            }
        }
    }

    // Get control gain at time t
    void getControlGain(double t, double& K1, double& K2) const {
        double P11, P12, P22;
        getP(t, P11, P12, P22);

        // For system with B = [0; 1], R = 1
        // K = R^(-1) * B^T * P
        K1 = _B[0] * P11 + _B[1] * P12;
        K2 = _B[0] * P12 + _B[1] * P22;
    }

    // Get size of stored solution
    size_t size() const {
        return _times.size();
    }
};

#endif // RICCATI_SOLVER_H
