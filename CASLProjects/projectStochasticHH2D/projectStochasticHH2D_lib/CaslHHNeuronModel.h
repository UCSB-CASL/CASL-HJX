// Created by Faranak Rajabi on 7/22/23.

#ifndef CASL4D_CASL_HH_NEURON_MODEL_H
#define CASL4D_CASL_HH_NEURON_MODEL_H

#include <vector>

class CaslHHNeuronModel {
public:
    // Constructor
    CaslHHNeuronModel();

    // Destructor
    ~CaslHHNeuronModel();

    // Member functions
    double calculateMInf                (double v);
    double calculateAM                  (double v);
    double calculateBM                  (double v);
    double calculateAN                  (double v);
    double calculateBN                  (double v);
    double calculateFV                  (double v , double n);
    double calculateFN                  (double v , double n);
    double calculateVdotSN              (double v , double n);
    double calculateNdotSN              (double v , double n);
    void calculateHHS                   (double t , double v , double n, std::vector<double>& dydt);
    void calculateHHSTwoNs              (double t , double v1, double n1, double v2, double n2, std::vector<double>& dydt);
    void calculateFvTwoNs               (double v1, double n1, double v2, double n2, double &fv1, double &fv2);
    void calculateFvTwoNsCoupled        (double v1, double n1, double v2, double n2, double alpha, double &fv1, double &fv2);
    void calculateFnTwoNsCoupled        (double v1, double n1, double v2, double n2, double &fn1, double &fn2);
    void calculateHHSTwoNsCoupled       (double t , double v1, double n1, double v2, double n2, double alpha, std::vector<double>& dydt);
    void calculateHHSTwoNsCoupledControl(double t , double v1, double n1, double v2, double n2, double alpha, double u1, double u2, std::vector<double>& dydt);
    void calculateZDyn                  (double t , double x , double y , double u , std::vector<double>& dydt);

private:
    // Constants
    const double Vna{50.0}  ;       // Na(Sodium) Voltage (mV)
    const double Vk{-77.0}  ;       // K(Potassium) Voltage (mV)
    const double Vl{-54.4};       // L(Leakage Channels) Voltage (mV)

    const double gna{120.0};       // Conductance of Na mS/cm^2
    const double gk{36.0}  ;       // Conductance of K mS/cm^2
    const double gl{0.3} ;       // Conductance of Leakage mS/cm^2

    const double c{1.0}         ;  // microF/cm^2
    const double Ib{10}       ;  // Neuron's baseline current (microA/cm^2), Bifurcation Param
    const double Ts{11.58}    ;  // Period of spiking (ms)
    const double TCondition{7.0};  // Time that control should take the neuron to target point
    const double Ks{100.0}      ;  // The scale factor, scaling the voltage values for numerical stability
};

#endif //CASL4D_CASL_HH_NEURON_MODEL_H
