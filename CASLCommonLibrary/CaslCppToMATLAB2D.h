//
// Created by Faranak Rajabi on 2/17/24.
//

#ifndef CASL_CPPTOMATLAB2D_H
#define CASL_CPPTOMATLAB2D_H

#include <string>
#include <iostream>
#include <fstream>

#include "CaslGrid2D.h"
#include "CaslArray2D.h"


class CaslCppToMATLAB2D {
private:

public:
     CaslCppToMATLAB2D();
    ~CaslCppToMATLAB2D() = default;

    void exportDataToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName);
    void exportSimInfoToMatlab(double currentTime, double dt, int iteration,
                               double phiValue, const std::string& fileName);
};


#endif //CASL_CPPTOMATLAB2D_H
