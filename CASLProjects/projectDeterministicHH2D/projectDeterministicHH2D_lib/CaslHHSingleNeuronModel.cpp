//
// Created by Faranak on 7/22/23.
//

/*
 * Implementation of the header file for single neuron model
 * in the absence of noise, coupling and control
 */

// Defined Constants
#include "CaslHHSingleNeuronModel.h"
#include <cmath>

double m_inf(double v)
{
    return am(v) / (am(v) + bm(v));
}

double am(double v)
{
    double nom = 0.1 * (v + 40.0);
    double denom = 1.0 - exp(-(v + 40.0) / 10.0);
    return nom / denom;
}

double bm(double v)
{
    return 4.0 * exp(-(v + 65.0) / 18.0);
}

double an(double v)
{
    double nom = 0.01 * (v + 55.0);
    double denom = 1.0 - exp(-(v + 55.0) / 10.0);
    return nom / denom;
}

double bn(double v)
{
    return 0.125 * exp(-(v + 65.0) / 80.0);
}

double fv(double v, double n)
{
    /*
     * dv/dt as a state dynamics for a single neuron without noise, coupling and strength
     * The input is a scalar v (NOT a vector for the population of neurons, ith neuron of N neurons)
     */
    double nom = Ib - gna * pow(m_inf(v), 3) * (0.8 - n) * (v - Vna) - gk * pow(n, 4) * (v - Vk) - gl * (v - Vl);
    double den = c;
    return nom / den;
}

double fn(double v, double n) {
    /*
     * dn/dt as a state dynamics for a single neuron without noise, coupling and control
     * The input is a scalar n (NOT a vector for the population of neurons, ith neuron of N neurons)
     */
    return an(v) * (1.0 - n) - bn(v) * n;
}

double vdotsn(double v, double n)
{
    return fv(v, n);
}

double ndotsn(double v, double n)
{
    return fn(v, n);
}

// HH model as a system of two odes
void hhs(double t, double v, double n, std::vector<double>& dydt)
{
    dydt[0] = vdotsn(v, n);
    dydt[1] = ndotsn(v, n);
}

void zdyn(double x, double y, double u, double& fvControl, double& fnControl){
    /*
     * State dynamics for a single neuron with control but without noise and coupling
     * The inputs are scalar(NOT a vector for the population of neurons, ith neuron of N neurons)
     * u should be the output of Bi-linear interpolation function when U as a CaslArray2D passed to it on the grid
     */
    double v = Ks * x;
    double n = y;
    fvControl = (1.0/Ks) * fv(v, n) + (1.0 / Ks) * u;
    fnControl = fn(n, v);
}