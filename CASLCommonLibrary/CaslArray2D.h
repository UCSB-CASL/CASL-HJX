//
// Created by Frederic Gibou on 12/28/22.
//
// Implements a CASLArray2D with the following structure in the x-direction (the y-direction is similar):
//
// In the schematic, we call L the number of pads in the negative x-direction (named _nPadsNegX in code).
// In the schematic, we call R the number of pads in the positive x-direction (named _nPadsPosX in code).
// In the schematic, we call nX the number of interior points in the x-direction (named _nX in code).
//
//       Pads in                                  Interior                               Pads in
//      NegativeX                                  Points                               PositiveX
// |------------------|   |-----------------------------------------------------|   |----------------|
//   o      o  ...   o      *      *      *      *   ...  *      *      *      *     o      o  ...  o
//   |               |      |                                                  |     |              |
// -L+1              0      1                                                  nX   nX+1           nX+R
//  =                                                                                               =
// iMin()                                                                                         iMax()

#ifndef CASL_ARRAY2D_H
#define CASL_ARRAY2D_H

#include <iostream>
#include <vector>

// todo: assert not out of bound.
// todo: add comments
// todo: setup unit tests.


#include "CaslOptions.h"

template<class T> class CaslArray2D {
private:
    std::vector<std::vector<T>> _data;         // Vector storing the data.
    int                         _nX;           // Number of interior points going from 1 to nX in x-direction.
    int                         _nY;           // Number of interior points going from 1 to nY in y-direction.
    int                         _nPadsNegX;    // Number of pads in the negative x-direction.
    int                         _nPadsPosX;    // Number of pads in the positive x-direction.
    int                         _nPadsNegY;    // Number of pads in the negative y-direction.
    int                         _nPadsPosY;    // Number of pads in the positive y-direction.

public:

    explicit CaslArray2D(int nX = 1, int nY = 1, int nPadsNegX = 0, int nPadsPosX = 0, int nPadsNegY = 0, int nPadsPosY = 0);
    CaslArray2D(const CaslArray2D<T>& rhs);
    virtual ~CaslArray2D();

    // Implements the access operator.
          T& operator()(const int& i, const int& j);
    const T& operator()(const int& i, const int& j) const;

    // Implements the operator for "this = that".
    // Vectors must have the same number of interior points and same number of positive and negative pads.
    CaslArray2D<T>& operator=(const CaslArray2D<T>& that);

    // Implements the operator for "this += that".
    // Vectors must have the same number of interior points and same number of positive and negative pads.
    CaslArray2D<T>& operator+=(const CaslArray2D<T>& that);

    // Implements the operator for "this + that".
    // Vectors must have the same number of interior points and same number of positive and negative pads.
    CaslArray2D<T> operator+ (const CaslArray2D<T>& that);

    // Implements the operator for "this -= that".
    // Vectors must have the same number of interior points and same number of positive and negative pads.
    CaslArray2D<T>& operator-=(const CaslArray2D<T>& that);

    // Implements the operator for "this - that".
    // Vectors must have the same number of interior points and same number of positive and negative pads.
    CaslArray2D<T> operator- (const CaslArray2D<T>& that);

    // Implements the operator for "this *= scalar".
    CaslArray2D<T>& operator*=(const T& scalar);

    // Implements the operator for "this * scalar".
    CaslArray2D<T> operator*(const T& scalar);

    // Implements the operator for "this /= scalar".
    CaslArray2D<T>& operator/=(const T& scalar);

    // Implements the operator for "this / scalar".
    CaslArray2D<T> operator/(const T& scalar);

    // Overload the << operator
    template <class U> friend std::ostream& operator<<(std::ostream& os, const CaslArray2D<U>& caslArray2D);

    /*
     * fillPaddingPoints: fill the padding points with values defined by the chosen method
     * @param [in]: with - the method used to define the values.
     * @param [in]: direction - the direction of the padding, i.e. the pads to the left, right, bottom, top
     *              By default we pad the outside of the entire domain.
     */
    void fillPaddingPoints(const CaslOptionPaddingWith& with, const CaslOptionPaddingDirection& direction = directionXY);

    /*
     * PadWithPeriodicX(): fill the padding points with values defined by periodicity in the x-direction.
     */
    void PadWithPeriodicX();

    /*
     * PadWithPeriodicY(): fill the padding points with values defined by periodicity in the y-direction.
     */
    void PadWithPeriodicY();

    /*
     * padWithConstantExtrapolationPosX(): fill the padding points with values defined by constant extrapolation in the positive x-direction.
     */
    void padWithConstantExtrapolationPosX();

    /*
     * padWithConstantExtrapolationNegX(): fill the padding points with values defined by constant extrapolation in the negative x-direction.
     */
    void padWithConstantExtrapolationNegX();

    /*
     * padWithConstantExtrapolationPosY(): fill the padding points with values defined by constant extrapolation in the positive y-direction.
     */
    void padWithConstantExtrapolationPosY();

    /*
     * padWithConstantExtrapolationNegY(): fill the padding points with values defined by constant extrapolation in the negative y-direction.
     */
    void padWithConstantExtrapolationNegY();

    /*
     * padWithLinearExtrapolationPosX(): fill the padding points with values defined by linear extrapolation in the positive x-direction.
     */
    void padWithLinearExtrapolationPosX();

    /*
     * padWithLinearExtrapolationNegX(): fill the padding points with values defined by linear extrapolation in the negative x-direction.
     */
    void padWithLinearExtrapolationNegX();

    /*
     * padWithLinearExtrapolationPosY(): fill the padding points with values defined by linear extrapolation in the positive y-direction.
     */
    void padWithLinearExtrapolationPosY();

    /*
     * padWithLinearExtrapolationNegY(): fill the padding points with values defined by linear extrapolation in the negative y-direction.
     */
    void padWithLinearExtrapolationNegY();

    /*
     * padWithQuadraticExtrapolationPosX(): fill the padding points with values defined by quadratic extrapolation in the positive x-direction.
     */
    void padWithQuadraticExtrapolationPosX();

    /*
     * padWithQuadraticExtrapolationNegX(): fill the padding points with values defined by quadratic extrapolation in the negative x-direction.
     */
    void padWithQuadraticExtrapolationNegX();

    /*
     * padWithQuadraticExtrapolationPosY(): fill the padding points with values defined by quadratic extrapolation in the positive y-direction.
     */
    void padWithQuadraticExtrapolationPosY();

    /*
     * padWithQuadraticExtrapolationNegY(): fill the padding points with values defined by quadratic extrapolation in the negative y-direction.
     */
    void padWithQuadraticExtrapolationNegY();

    /*
     * assertSameDimensionsAs(): asserts that this = that in terms of their dimensions and number of pads
     */
    void assertSameDimensionsAs(const CaslArray2D<T>& that);

    /*
     * sizeDirectionX(): returns the size (number of interior points + pads) of this in the x-direction
     */
    int sizeDirectionX() const;

    /*
     * sizeDirectionY(): returns the size (number of interior points + pads) of this in the y-direction
     */
    int sizeDirectionY() const;

    /*
     * iMin(): returns the smallest index (can be negative) in the x-direction
     */
    int           iMin() const;

    /*
     * iMax(): returns the largest index (can be negative) in the x-direction
     */
    int           iMax() const;

    /*
     * jMin(): returns the smallest index (can be negative) in the y-direction
     */
    int           jMin() const;

    /*
     * jMax(): returns the largest index (can be negative) in the y-direction
     */
    int           jMax() const;

    /*
     * nX(): returns the number of interior points in the x-direction
     */
    int             nX() const;

    /*
     * nY(): returns the number of interior points in the y-direction
     */
    int             nY() const;

    /*
     * nPadsNegX(): returns the number of pads in the negative x-direction
     */
    int       nPadsNegX() const;

    /*
     * nPadsPosX(): returns the number of pads in the positive x-direction
     */
    int       nPadsPosX() const;

    /*
     * nPadsNegY(): returns the number of pads in the negative y-direction
     */
    int       nPadsNegY() const;

    /*
     * nPadsPosY(): returns the number of pads in the positive y-direction
     */
    int       nPadsPosY() const;

    T maxAbs() const;
};

#include "CaslArray2D.cpp"

#endif // CASL_ARRAY2D_H
