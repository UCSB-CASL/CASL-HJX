//
// Created by Faranak on 7/22/23.
//

/*
 * Implementation of the header file for single neuron model
 * in the absence of noise, coupling and control
 */

#include "CaslHHSingleNeuronModel.h"
#include <cmath>
#include <limits>

namespace {
    constexpr double EPS_HH = 1e-12;
}

double am(double v)
{
    // Removable singularity at v = -40:
    // lim_{v->-40} 0.1(v+40)/(1-exp(-(v+40)/10)) = 1
    if (std::fabs(v + 40.0) < EPS_HH) return 1.0;

    double nom = 0.1 * (v + 40.0);
    double denom = 1.0 - std::exp(-(v + 40.0) / 10.0);
    return nom / denom;
}

double bm(double v)
{
    return 4.0 * std::exp(-(v + 65.0) / 18.0);
}

double an(double v)
{
    // Removable singularity at v = -55:
    // lim_{v->-55} 0.01(v+55)/(1-exp(-(v+55)/10)) = 0.1
    if (std::fabs(v + 55.0) < EPS_HH) return 0.1;

    double nom = 0.01 * (v + 55.0);
    double denom = 1.0 - std::exp(-(v + 55.0) / 10.0);
    return nom / denom;
}

double bn(double v)
{
    return 0.125 * std::exp(-(v + 65.0) / 80.0);
}

double m_inf(double v)
{
    const double a = am(v);
    const double b = bm(v);
    const double denom = a + b;

    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(denom) || std::fabs(denom) < EPS_HH) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return a / denom;
}

double fv(double v, double n)
{
    /*
     * dv/dt as a state dynamics for a single neuron without noise, coupling and control
     */
    const double minf = m_inf(v);
    const double nom =
        Ib
        - gna * std::pow(minf, 3) * (0.8 - n) * (v - Vna)
        - gk  * std::pow(n, 4)    * (v - Vk)
        - gl  * (v - Vl);

    return nom / c;
}

double fn(double v, double n)
{
    /*
     * dn/dt as a state dynamics for a single neuron without noise, coupling and control
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

// HH model as a system of two ODEs
void hhs(double t, double v, double n, std::vector<double>& dydt)
{
    (void)t;
    dydt[0] = vdotsn(v, n);
    dydt[1] = ndotsn(v, n);
}

void zdyn(double x, double y, double u, double& fvControl, double& fnControl)
{
    /*
     * State dynamics for a single neuron with control but without noise and coupling
     */
    double v = Ks * x;
    double n = y;

    fvControl = (1.0 / Ks) * fv(v, n) + (1.0 / Ks) * u;
    fnControl = fn(v, n);   // FIXED: was fn(n, v)
}