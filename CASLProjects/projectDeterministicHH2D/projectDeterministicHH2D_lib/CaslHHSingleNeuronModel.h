//
// Created by Faranak on 7/22/23.
//

#ifndef CASLUNIFORM_CASLHHSINGLENEURONMODEL_H
#define CASLUNIFORM_CASLHHSINGLENEURONMODEL_H

#include <vector>

// Defined Constants
const double Vna = 50.0;               // Na(Sodium) Voltage (mV)
const double Vk = -77.0;               // K(Potassium) Voltage (mV)
const double Vl = -54.4;             // L(Leakage Channels) Voltage (mV)

const double gna = 120.0;              // Conductance of Na mS/cm^2
const double gk = 36.0;                // Conductance of K mS/cm^2
const double gl = 0.3;               // Conductance of Leakage mS/cm^2

const double c = 1.0;                  // microF/cm^2
const double Ib = 10.0;                // Neuron's baseline current (microA/cm^2), Bifurcation Param
const double Ts = 11.58;             // Period of spiking (ms)
const double TCondition = 7.0;         // Time that control should take the neuron to target point
const double Ks = 100.0;                   // The scale factor, scaling the voltage values for numerical stability

double m_inf(double v);
double am(double v);
double bm(double v);
double an(double v);
double bn(double v);
double fv(double v, double n);
double fn(double v, double n);
double vdotsn(double v, double n);
double ndotsn(double v, double n);
void hhs(double t, double v, double n, std::vector<double>& dydt);
void zdyn(double x, double y, double u, double& fvControl, double& fnControl);

#endif //CASLUNIFORM_CASLHHSINGLENEURONMODEL_H
