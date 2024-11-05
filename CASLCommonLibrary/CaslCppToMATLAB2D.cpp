//
// Created by Faranak Rajabi on 2/17/24.
//

#include "CaslCppToMATLAB2D.h"

CaslCppToMATLAB2D::CaslCppToMATLAB2D() {}

void CaslCppToMATLAB2D::exportDataToMatlab(const CaslGrid2D& grid, const CaslArray2D<double>& un, const std::string& fileName) {
   int nX = grid.nX(), nY = grid.nY();
   std::ofstream ofStream; ofStream.open(fileName);

   if (ofStream.is_open()) {
        for (int i = 1; i <= nX; i++) {
            for (int j = 1; j <= nY; j++) {
                ofStream << un(i, j) << " ";
            }
            ofStream << std::endl;
        }
        ofStream.close();
        return;
   }

   // If file could not be open:
   std::cout << "File " << fileName << " could not be opened. EXITING." << std::endl; exit(1);
}

void CaslCppToMATLAB2D::exportSimInfoToMatlab(double currentTime, double dt, int iteration,
                           double phiValue, const std::string& fileName) {
    std::ofstream ofStream;
    ofStream.open(fileName, std::ios::app);  // append mode

    if (ofStream.is_open()) {
        ofStream << std::fixed;  // for fixed-point notation
        ofStream << "Export to Matlab at currentTime = " << currentTime << std::endl;
        ofStream << "dt: " << std::setprecision(16) << dt << std::endl;
        ofStream << "phi_" << iteration << ": " << std::setprecision(16) << phiValue << std::endl;
        ofStream << "----------------------------------------" << std::endl;
        ofStream.close();
        return;
    }

    // If file couldn't be opened
    std::cout << "File " << fileName << " could not be opened. EXITING." << std::endl;
    exit(1);
}


