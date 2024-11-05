#ifndef CASL_LQR_SYSTEM_DYNAMICS_H
#define CASL_LQR_SYSTEM_DYNAMICS_H

#include "CaslGrid2D.h"
#include "CaslArray2D.h"

class CaslLQRSystemDynamics {
private:
    void solveRiccati(double t, double& P11, double& P12, double& P22) const;

public:
    CaslLQRSystemDynamics();
    ~CaslLQRSystemDynamics();

    void LQRDynamics(double x1, double x2, double &fx1, double &fx2);
    double exactSolution(double x1, double x2, double t) const;
    double getOptimalControl(double x1, double x2, double t) const;
    double computeMaxError(const CaslGrid2D& grid,
                           const CaslArray2D<double>& numericalSolution,
                           double t);
    double computeL2Error(const CaslGrid2D& grid,
                          const CaslArray2D<double>& numericalSolution,
                          double t);
    double verifyHJB(double x1, double x2, double t) const;
};

#endif //CASL_LQR_SYSTEM_DYNAMICS_H