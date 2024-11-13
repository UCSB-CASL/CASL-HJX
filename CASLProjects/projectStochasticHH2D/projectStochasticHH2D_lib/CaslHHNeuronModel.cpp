// Created by Faranak Rajabi on 7/22/23.

#include "CaslHHNeuronModel.h"
#include <cmath>

CaslHHNeuronModel::CaslHHNeuronModel() {
}

CaslHHNeuronModel::~CaslHHNeuronModel() {
}

double CaslHHNeuronModel::calculateMInf(double v) {
    return calculateAM(v) / (calculateAM(v) + calculateBM(v));
}

double CaslHHNeuronModel::calculateAM(double v) {
    double nom = 0.1 * (v + 40);
    double denom = 1.0 - exp(-(v + 40) / 10.0);
    return nom / denom;
}

double CaslHHNeuronModel::calculateBM(double v) {
    return 4.0 * exp(-(v + 65) / 18.0);
}

double CaslHHNeuronModel::calculateAN(double v) {
    double nom = 0.01 * (v + 55);
    double denom = 1.0 - exp(-(v + 55) / 10.0);
    return nom / denom;
}

double CaslHHNeuronModel::calculateBN(double v) {
    return 0.125 * exp(-(v + 65) / 80.0);
}

double CaslHHNeuronModel::calculateFV(double v, double n) {
    double nom = Ib - gna * pow(calculateMInf(v), 3) * (0.8 - n) * (v - Vna) - gk * pow(n, 4) * (v - Vk) - gl * (v - Vl);
    double den = c;
    return nom / den;
}

double CaslHHNeuronModel::calculateFN(double v, double n) {
    return calculateAN(v) * (1 - n) - calculateBN(v) * n;
}

double CaslHHNeuronModel::calculateVdotSN(double v, double n) {
    /*
     * dv/dt as a state dynamics for a single neuron without noise, coupling and strength
     * The input is a scalar v (NOT a vector for the population of neurons, ith neuron of N neurons)
     */
    return calculateFV(v, n);
}

double CaslHHNeuronModel::calculateNdotSN(double v, double n) {
    /*
     * dn/dt as a state dynamics for a single neuron without noise, coupling and control
     * The input is a scalar n (NOT a vector for the population of neurons, ith neuron of N neurons)
     */
    return calculateFN(v, n);
}

// HH model as a system of two ODEs(importable for a RK4 solverr)
void CaslHHNeuronModel::calculateHHS(double t, double v, double n, std::vector<double>& dydt) {
    dydt[0] = calculateVdotSN(v, n);
    dydt[1] = calculateNdotSN(v, n);
}

void CaslHHNeuronModel::calculateZDyn(double t, double x, double y, double u, std::vector<double>& dydt) {
    /*
     * State dynamics for a single neuron with control but without noise and coupling
     * The inputs are scalar(NOT a vector for the population of neurons, ith neuron of N neurons)
     * u should be the output of Bi-linear interpolation function when U as a CaslArray2D passed to it on the grid
     */
    double v = Ks * x, n = y;
    dydt[0] = (1.0 / Ks) * calculateFV(v, n) + (1.0 / Ks) * u;
    dydt[1] = calculateFN(n, v);
}

// HH model as a system of two ODEs for two neurons without coupling and control(importable to a RK4 solver)
void CaslHHNeuronModel::calculateHHSTwoNs(double t, double v1, double n1, double v2, double n2, std::vector<double>& dydt) {
    dydt[0] = calculateVdotSN(v1, n1);
    dydt[1] = calculateNdotSN(v1, n1);
    dydt[2] = calculateVdotSN(v2, n2);
    dydt[3] = calculateNdotSN(v2, n2);
}

// Fv for two coupled neurons without control and coupling(and without time as an input)
void CaslHHNeuronModel::calculateFvTwoNs(double v1, double n1, double v2, double n2, double &fv1, double &fv2) {
    fv1 = calculateVdotSN(v1, n1);
    fv2 = calculateVdotSN(v2, n2);
}

// Fv for two coupled neurons without control(and without time as an input)
void CaslHHNeuronModel::calculateFvTwoNsCoupled(double v1, double n1, double v2, double n2, double alpha, double &fv1, double &fv2) {
    fv1 = calculateVdotSN(v1, n1) + alpha * (v2 - v1);
    fv2 = calculateVdotSN(v2, n2) + alpha * (v1 - v2);
}

// Fn for two coupled neurons without control(and without time as an input)
void CaslHHNeuronModel::calculateFnTwoNsCoupled(double v1, double n1, double v2, double n2, double &fn1, double &fn2) {
    fn1 = calculateNdotSN(v1, n1);
    fn2 = calculateNdotSN(v2, n2);
}

// HH model as a system of two ODEs for two coupled neurons without control(importable to a RK4 solver)
void CaslHHNeuronModel::calculateHHSTwoNsCoupled(double t, double v1, double n1, double v2, double n2, double alpha, std::vector<double>& dydt) {
    dydt[0] = calculateVdotSN(v1, n1) + alpha * (v2 - v1);
    dydt[1] = calculateNdotSN(v1, n1);
    dydt[2] = calculateVdotSN(v2, n2) + alpha * (v1 - v2);
    dydt[3] = calculateNdotSN(v2, n2);
}

// HH model as a system of two ODEs for two coupled neurons without control(importable to a RK4 solver)
void CaslHHNeuronModel::calculateHHSTwoNsCoupledControl(double t, double v1, double n1, double v2, double n2, double alpha,
                                                        double u1, double u2, std::vector<double>& dydt) {
    dydt[0] = calculateVdotSN(v1, n1) + alpha * (v2 - v1) + u1;
    dydt[1] = calculateNdotSN(v1, n1);
    dydt[2] = calculateVdotSN(v2, n2) + alpha * (v1 - v2) + u2;
    dydt[3] = calculateNdotSN(v2, n2);
}

