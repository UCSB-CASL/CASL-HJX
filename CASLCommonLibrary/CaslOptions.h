//
// Created by Frederic Gibou on 6/26/22.
//

#ifndef CASL_OPTIONS_H
#define CASL_OPTIONS_H

enum CaslOptionPaddingDirection         {directionPositiveX, directionNegativeX, directionX, directionPositiveY, directionNegativeY, directionY, directionXY};
enum CaslOptionPaddingWith              {withConstantExtrapolation, withLinearExtrapolation, withQuadraticExtrapolation, withPeriodicCondition};
enum CaslOptionNumericalFirstDerivative {Upwind, ENO2, ENO3, WENO5};
enum CaslOptionNumericalSecondDerivative{ForwardTimeCentralSpacing, BackwardTimeCentralSpacing, CrankNicolson};
enum CaslOptionSecondOrderTermDirection {X, Y, XY};

#endif //CASL_OPTIONS_H