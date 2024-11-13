//
// Created by Faranak on 7/22/23.
//

#ifndef CASLUNIFORM_CASLHHWHOLESYSTEMMODEL_H
#define CASLUNIFORM_CASLHHWHOLESYSTEMMODEL_H

#include "CaslHHSingleNeuronModel.h"
#include <random>
#include <chrono>

//Defined Constants
const int N = 100;                  // Number of Neurons
const double alpha_mu = 0.1;        // Mean value for random normal distribution for Coupling strength(alpha)
const double alpha_sigma = 0.02;    // Sigma value for random normal distribution for Coupling strength parameter sigma
const double eta_mu = 0.0;            // Mean value for random normal distribution for neuron's noise(eta)
const double eta_sigma = 1.0;         // D value(for sigma) for random normal distribution for neuron's noise(eta)

// Generator class for generating random numbers with a normal distribution.
class Generator {
    std::default_random_engine generator; // Random number generator engine.
    std::normal_distribution<double> distribution; // Normal distribution.
    double min; // Minimum value of the distribution range.
    double max; // Maximum value of the distribution range.

public:
    // Constructor to initialize the Generator with mean, stddev, min, and max.
    Generator(double mean, double stddev, double min, double max) :
            distribution(mean, stddev), min(min), max(max) {

        auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        generator.seed(seed); // Seed the random number generator.
    }

    // Operator function to generate a random number using the distribution.
    double operator()() {
        double number = distribution(generator);
        return number;
    }
};

void vdot(double v[], double n[], double u, std::vector<double>& vDotsVec);
void ndot(double v[], double n[], std::vector<double>& nDotsVec);
//void HH_sys(double v[], double n[], double u, vector<double> vdot_out, vector<double> ndot_out);

#endif //CASLUNIFORM_CASLHHWHOLESYSTEMMODEL_H
