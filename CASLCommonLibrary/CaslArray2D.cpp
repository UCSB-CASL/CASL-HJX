//
// Created by Frederic Gibou on 12/28/22.
//

#ifndef CASL_ARRAY2D_CPP
#define CASL_ARRAY2D_CPP

#include "CaslArray2D.h"
#include "CaslOptions.h"

/*
 * Parameter Constructor
 */
template<class T> CaslArray2D<T>::CaslArray2D(int nX, int nY, const int nPadsNegX, const int nPadsPosX, const int nPadsNegY, const int nPadsPosY) {
    _nX        = nX;
    _nY        = nY;
    _nPadsNegX = nPadsNegX;
    _nPadsPosX = nPadsPosX;
    _nPadsNegY = nPadsNegY;
    _nPadsPosY = nPadsPosY;
    _data.resize(_nX + _nPadsNegX + _nPadsPosX);   // Number of columns of the grid.
    for (int col = 0; col < _data.size(); ++col) _data[col].resize(nY + nPadsNegY + nPadsPosY, 0);
}

/*
 * Copy Constructor
 */
template<class T> CaslArray2D<T>::CaslArray2D(const CaslArray2D<T>& rhs) {
    _data       = rhs._data;
    _nX         = rhs._nX;
    _nY         = rhs._nY;
    _nPadsNegX  = rhs._nPadsNegX;
    _nPadsPosX  = rhs._nPadsPosX;
    _nPadsNegY  = rhs._nPadsNegY;
    _nPadsPosY  = rhs._nPadsPosY;
}

/*
 * (Virtual) Destructor
 */
template<class T> CaslArray2D<T>::~CaslArray2D() = default;

/*
 * Access the individual elements
 */
template<class T> T& CaslArray2D<T>::operator()(const int& i, const int& j) {
    return _data[i - 1 + _nPadsNegX][j - 1 + _nPadsNegY];
}

/*
 * Access the individual elements (const)
 */
template<class T> const T& CaslArray2D<T>::operator()(const int& i, const int& j) const {
    return _data[i - 1 + _nPadsNegX][j - 1 + _nPadsNegY];
}

/*
 * CaslArray2D this = that
 */
template<class T>
CaslArray2D<T>& CaslArray2D<T>::operator=(const CaslArray2D<T>& that)
{
    if (this != &that) { // Check that it is not a self-copy.
        _nX        = that._nX;
        _nY        = that._nY;
        _nPadsNegX = that._nPadsNegX;
        _nPadsPosX = that._nPadsPosX;
        _nPadsNegY = that._nPadsNegY;
        _nPadsPosY = that._nPadsPosY;
        _data.resize(_nX + _nPadsNegX + _nPadsPosX);   // Number of columns of the grid.
        for (int col = 0; col < _data.size(); ++col) _data[col].resize(_nY + _nPadsNegY + _nPadsPosY, 0);
        for (int i = 0; i < this->sizeDirectionX(); ++i)
            for (int j = 0; j < this->sizeDirectionY(); ++j)
                this->_data[i][j] = that._data[i][j];
    }
    return *this;
}

/*
 * CaslArray2D this += that
 */
template<class T> CaslArray2D<T>& CaslArray2D<T>::operator+=(const CaslArray2D<T>& that) {
    this->assertSameDimensionsAs(that);
    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            this->_data[i][j] += that._data[i][j];
    return *this;
}

/*
 * CaslArray2D this -= that
 */
template<class T> CaslArray2D<T>& CaslArray2D<T>::operator-=(const CaslArray2D<T>& that) {
    this->assertSameDimensionsAs(that);
    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            this->_data[i][j] -= that._data[i][j];
    return *this;
}

/*
 * CaslArray2D this *= that
 */
template<class T> CaslArray2D<T>& CaslArray2D<T>::operator*=(const T& scalar) {
    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            this->_data[i][j] *= scalar;
    return *this;
}

/*
 * CaslArray2D this /= that
 */
template<class T> CaslArray2D<T>& CaslArray2D<T>::operator/=(const T& scalar) {
    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            this->_data[i][j] /= scalar;
    return *this;
}

/*
 * CaslArray2D this + that.
 */
template<class T>
CaslArray2D<T> CaslArray2D<T>::operator+ (const CaslArray2D<T>& that){
    this->assertSameDimensionsAs(that);
    CaslArray2D<T> result(_nX, _nY, _nPadsNegX, _nPadsPosX, _nPadsNegY, _nPadsPosY);

    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            result._data[i][j] = this->_data[i][j] + that._data[i][j];

    return result;
}

/*
 * CaslArray2D this - that.
 */
template<class T>
CaslArray2D<T> CaslArray2D<T>::operator- (const CaslArray2D<T>& that){
    this->assertSameDimensionsAs(that);
    CaslArray2D<T> result(_nX, _nY, _nPadsNegX, _nPadsPosX, _nPadsNegY, _nPadsPosY);

    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            result._data[i][j] = this->_data[i][j] - that._data[i][j];

    return result;
}

/*
 * CaslArray2D result = this / scalar
 */
template<class T>
CaslArray2D<T> CaslArray2D<T>::operator/ (const T& scalar){
    CaslArray2D<T> result(_nX, _nY, _nPadsNegX, _nPadsPosX, _nPadsNegY, _nPadsPosY);

    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            result._data[i][j] = scalar / this->_data[i][j];
    return result;
}

/*
 * CaslArray2D result = this * scalar
 */
template<class T>
CaslArray2D<T> CaslArray2D<T>::operator* (const T& scalar){
    CaslArray2D<T> result(_nX, _nY, _nPadsNegX, _nPadsPosX, _nPadsNegY, _nPadsPosY);

    for (int i = 0; i < this->sizeDirectionX(); ++i)
        for (int j = 0; j < this->sizeDirectionY(); ++j)
            result._data[i][j] = scalar * this->_data[i][j];
    return result;
}

/*
 * Non-member multiplication by a scalar from left
 */
template<class T>
CaslArray2D<T> operator* ( T scalar, CaslArray2D<T>& array ) {
    return array * scalar;
}

/*
 * Non-member multiplication by a scalar from right
 */
template<class T>
CaslArray2D<T> operator* (CaslArray2D<T>& array, T scalar) {
    return array * scalar;
}

/*
 * Non-member division by a scalar from left
 */
template<class T>
CaslArray2D<T> operator/ ( T scalar, CaslArray2D<T>& array ) {
    return array / scalar;
}

/*
 * Non-member division by a scalar from right
 */
template<class T>
CaslArray2D<T> operator/ (CaslArray2D<T>& array, T scalar) {
    return array / scalar;
}

/*
 * Overload <<
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const CaslArray2D<T>& caslArray2D) {
    int nPadsNegX = caslArray2D._nPadsNegX, nPadsPosX = caslArray2D._nPadsPosX;
    int nPadsNegY = caslArray2D._nPadsNegY, nPadsPosY = caslArray2D._nPadsPosY;
    int nX = caslArray2D._nX, nY = caslArray2D._nY;

    // Bottom pad
    for (int j = 1; j <= nPadsNegY; ++j) {
        for (int i = 1 - nPadsNegX; i <= nX + nPadsPosX; ++i) {
            os << caslArray2D(i, j) << "\t";
        }
        os << std::endl;
    }

    // Interior points
    for (int j = 1; j <= nY; ++j) {
        for (int i = 1 - nPadsNegX; i <= nX + nPadsPosX; ++i) {
            os << caslArray2D(i, j) << "\t";
        }
        os << std::endl;
    }

    // Top pad
    for (int j = nY + 1; j <= nY + nPadsPosY; ++j) {
        for (int i = 1 - nPadsNegX; i <= nX + nPadsPosX; ++i) {
            os << caslArray2D(i, j) << "\t";
        }
        os << std::endl;
    }

    return os;
}

/*
 * fillPaddingPoints: fill the padding points with values defined by the chosen method
 * @param [in]: with - the method used to define the values.
 * @param [in]: direction - the direction of the padding, i.e. the pads to the left, right, bottom, top
 *              By default we pad the outside of the entire domain.
 */
template<class T>
void CaslArray2D<T>::fillPaddingPoints(const CaslOptionPaddingWith& with, const CaslOptionPaddingDirection& direction) {
    if ((direction == directionNegativeX || direction == directionPositiveX) && with == withPeriodicCondition) {
        std::cout << "CASL ERROR in CASLArray2D::fillPaddingPoints():" << std::endl <<
                  "\t periodic can only be used with directionX, not directionNegativeX or directionPositiveX." <<
                  "\n Exiting." << std::endl;
        exit(1);
    }

    if ((direction == directionNegativeY || direction == directionPositiveY) && with == withPeriodicCondition) {
        std::cout << "CASL ERROR in CASLArray2D::fillPaddingPoints():" << std::endl <<
                  "\t periodic can only be used with directionY, not directionNegativeY or directionPositiveY." <<
                  "\n Exiting." << std::endl;
        exit(1);
    }

    // Constant extrapolation case:
    if (with == withConstantExtrapolation) {
        if (direction == directionPositiveX || direction == directionX || direction == directionXY)
            padWithConstantExtrapolationPosX();
        if (direction == directionPositiveY || direction == directionY || direction == directionXY)
            padWithConstantExtrapolationPosY();

        if (direction == directionNegativeX || direction == directionX || direction == directionXY)
            padWithConstantExtrapolationNegX();
        if (direction == directionNegativeY || direction == directionY || direction == directionXY)
            padWithConstantExtrapolationNegY();

        return;
    }

    // Linear extrapolation case:
    if (with == withLinearExtrapolation) {
        if (direction == directionPositiveX || direction == directionX || direction == directionXY)
            padWithLinearExtrapolationPosX();
        if (direction == directionPositiveY || direction == directionY || direction == directionXY)
            padWithLinearExtrapolationPosY();

        if (direction == directionNegativeX || direction == directionX || direction == directionXY)
            padWithLinearExtrapolationNegX();
        if (direction == directionNegativeY || direction == directionY || direction == directionXY)
            padWithLinearExtrapolationNegY();

        return;
    }

    // Quadratic extrapolation case:
    if (with == withQuadraticExtrapolation) {
        if (direction == directionPositiveX || direction == directionX || direction == directionXY)
            padWithQuadraticExtrapolationPosX();
        if (direction == directionPositiveY || direction == directionY || direction == directionXY)
            padWithQuadraticExtrapolationPosY();

        if (direction == directionNegativeX || direction == directionX || direction == directionXY)
            padWithQuadraticExtrapolationNegX();
        if (direction == directionNegativeY || direction == directionY || direction == directionXY)
            padWithQuadraticExtrapolationNegY();

        return;
    }

    // Periodic case:
    if (with == withPeriodicCondition) {
        if (direction == directionX )   PadWithPeriodicX();
        if (direction == directionY )   PadWithPeriodicY();
        if (direction == directionXY) { PadWithPeriodicX(); PadWithPeriodicY();}
        return;
    }

    // If none of the above, return an error message and exit:
    std::cout << "CASL ERROR in CASLArray2D::fillPaddingPoints():" << std::endl <<
              "\t the combination of " << direction << " and " << with <<
              " is not implemented." <<
              "\n Exiting." << std::endl;
    exit(1);
}

/*
 * PadWithPeriodicX(): fill the padding points with values defined by periodicity in the x-direction.
 */
template<class T>
void CaslArray2D<T>::PadWithPeriodicX() {
    int offsetRight = this->sizeDirectionX() - _nPadsPosX;
    int offsetLeft  = _nPadsNegX + 1;

    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        // Set the value(s) of the xPositive padding:
        for (int i = 0; i < _nPadsPosX; ++i) _data[offsetRight + i][j] = _data[offsetLeft + i][j];

        // Set the value(s) of the xNegative padding:
        int offset = this->sizeDirectionX() - 1 - _nPadsNegX - _nPadsPosX;
        for (int i = 0; i < _nPadsNegX; ++i) _data[i][j] = _data[offset + i][j];
    }
}

/*
 * PadWithPeriodicY(): fill the padding points with values defined by periodicity in the y-direction.
 */
template<class T>
void CaslArray2D<T>::PadWithPeriodicY() {
    int offsetTop = this->sizeDirectionY() - _nPadsPosY;
    int offsetBottom  = _nPadsNegY + 1;

    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        // Set the value(s) of the yPositive padding:
        for (int j = 0; j < _nPadsPosY; ++j) _data[i][offsetTop + j] = _data[i][offsetBottom + j];

        // Set the value(s) of the xNegative padding:
        int offset = this->sizeDirectionY() - 1 - _nPadsNegY - _nPadsPosY;
        for (int j = 0; j < _nPadsNegY; ++j) _data[i][j] = _data[i][offset + j];
    }
}

/*
 * padWithConstantExtrapolationPosX(): fill the padding points with values defined by constant extrapolation in the positive x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithConstantExtrapolationPosX() {
    // Set the value(s) of the xPositive padding:
    int offset = this->sizeDirectionX() - _nPadsPosX;
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        T valueToCopy = _data[offset - 1][j];
        for (int i = 0; i < _nPadsPosX; ++i) _data[offset + i][j] = valueToCopy;
    }
}

/*
 * padWithConstantExtrapolationPosY(): fill the padding points with values defined by constant extrapolation in the positive y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithConstantExtrapolationPosY() {
    // Set the value(s) of the yPositive padding:
    int offset = this->sizeDirectionY() - _nPadsPosY;
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        T valueToCopy = _data[i][offset - 1];
        for (int j = 0; j < _nPadsPosY; ++j) _data[i][offset + j] = valueToCopy;
    }
}

/*
 * padWithConstantExtrapolationNegX(): fill the padding points with values defined by constant extrapolation in the negative x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithConstantExtrapolationNegX() {
    // Set the value(s) of the xNegative padding:
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        T valueToCopy = _data[_nPadsNegX][j];
        for (int i = _nPadsNegX - 1; i >= 0; --i) _data[i][j] = valueToCopy;
    }
}

/*
 * padWithConstantExtrapolationNegY(): fill the padding points with values defined by constant extrapolation in the negative y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithConstantExtrapolationNegY() {
    // Set the value(s) of the yNegative padding:
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        T valueToCopy = _data[i][_nPadsNegY];
        for (int j = _nPadsNegY - 1; j >= 0; --j) _data[i][j] = valueToCopy;
    }
}

/*
 * padWithLinearExtrapolationPosX(): fill the padding points with values defined by linear extrapolation in the positive x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithLinearExtrapolationPosX() {
    // Set the value(s) of the xPositive padding:
    int offset = this->sizeDirectionX() - _nPadsPosX;
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        for (int i = 0; i < _nPadsPosX; ++i) _data[offset + i][j] = 2 * _data[offset + i - 1][j] - _data[offset + i - 2][j];
    }
}


/*
 * padWithLinearExtrapolationPosY(): fill the padding points with values defined by linear extrapolation in the positive y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithLinearExtrapolationPosY() {
    // Set the value(s) of the yPositive padding:
    int offset = this->sizeDirectionY() - _nPadsPosY;
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        for (int j = 0; j < _nPadsPosY; ++j) _data[i][offset + j] = 2 * _data[i][offset + j - 1] - _data[i][offset + j - 2];
    }
}

/*
 * padWithLinearExtrapolationNegX(): fill the padding points with values defined by linear extrapolation in the negative x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithLinearExtrapolationNegX() {
    // Set the value(s) of the xNegative padding:
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        for (int i = _nPadsNegX - 1; i >= 0; --i) _data[i][j] = 2 * _data[i + 1][j] - _data[i + 2][j];
    }
}


/*
 * padWithLinearExtrapolationNegY(): fill the padding points with values defined by linear extrapolation in the negative y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithLinearExtrapolationNegY() {
    // Set the value(s) of the yNegative padding:
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        for (int j = _nPadsNegY - 1; j >= 0; --j) _data[i][j] = 2 * _data[i][j + 1] - _data[i][j + 2];
    }
}

/*
 * padWithQuadraticExtrapolationPosX(): fill the padding points with values defined by quadratic extrapolation in the positive x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithQuadraticExtrapolationPosX() {
    // Set the value(s) of the xPositive padding:
    int offset = this->sizeDirectionX() - _nPadsPosX;
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        for (int i = 0; i < _nPadsPosX; ++i) _data[offset + i][j] = 3 * _data[offset + i - 1][j] - 3 * _data[offset + i - 2][j] + _data[offset + i - 3][j];
    }
}

/*
 * padWithQuadraticExtrapolationPosY(): fill the padding points with values defined by quadratic extrapolation in the positive y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithQuadraticExtrapolationPosY() {
    // Set the value(s) of the yPositive padding:
    int offset = this->sizeDirectionY() - _nPadsPosY;
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        for (int j = 0; j < _nPadsPosY; ++j) _data[i][offset + j] = 3 * _data[i][offset + j - 1] - 3 * _data[i][offset + j - 2] + _data[i][offset + j - 3];
    }
}

/*
 * padWithQuadraticExtrapolationNegX(): fill the padding points with values defined by quadratic extrapolation in the negative x-direction.
 */
template<class T>
void CaslArray2D<T>::padWithQuadraticExtrapolationNegX() {
    // Set the value(s) of the xNegative padding:
    for (int j = _nPadsNegY; j < _nPadsNegY + _nY; ++j) {
        for (int i = _nPadsNegX - 1; i >= 0; --i) _data[i][j] = 3 * _data[i + 1][j] - 3 * _data[i + 2][j] + _data[i + 3][j];
    }
}

/*
 * padWithQuadraticExtrapolationNegY(): fill the padding points with values defined by quadratic extrapolation in the negative y-direction.
 */
template<class T>
void CaslArray2D<T>::padWithQuadraticExtrapolationNegY() {
    // Set the value(s) of the yNegative padding:
    for (int i = _nPadsNegX; i < _nPadsNegX + _nX; ++i) {
        for (int j = _nPadsNegY - 1; j >= 0; --j) _data[i][j] = 3 * _data[i][j + 1] - 3 * _data[i][j + 2] + _data[i][j + 3];
    }
}

/*
 * assertSameDimensionsAs(): asserts that this = that in terms of their dimensions and number of pads
 */
template<class T>
void CaslArray2D<T>::assertSameDimensionsAs(const CaslArray2D<T>& that) {
    if ( _nX        != that._nX        ||
         _nPadsNegX != that._nPadsNegX ||
         _nPadsPosX != that._nPadsPosX ||
         _nY        != that._nY        ||
         _nPadsNegY != that._nPadsNegY ||
         _nPadsPosY != that._nPadsPosY) {
        std::cout
                << "CASL ERROR in CaslArray2D assertSameDimensionsAs - DIMENSIONS DO NOT AGREE BETWEEN THIS AND THAT:"
                << std::endl <<
                "\t THIS number of interior points in x-direction = " << _nX              << std::endl <<
                "\t THAT number of interior points in x-direction = " << that.nX()        << std::endl <<
                "\t THIS number of positive pads   in x-direction = " << _nPadsPosX       << std::endl <<
                "\t THAT number of positive pads   in x-direction = " << that.nPadsPosX() << std::endl <<
                "\t THIS number of negative pads   in x-direction = " << _nPadsNegX       << std::endl <<
                "\t THAT number of negative pads   in x-direction = " << that.nPadsNegX() << std::endl <<
                "\t THIS number of interior points in y-direction = " << _nY              << std::endl <<
                "\t THAT number of interior points in y-direction = " << that.nY()        << std::endl <<
                "\t THIS number of positive pads   in y-direction = " << _nPadsPosY       << std::endl <<
                "\t THAT number of positive pads   in y-direction = " << that.nPadsPosY() << std::endl <<
                "\t THIS number of negative pads   in y-direction = " << _nPadsNegY       << std::endl <<
                "\t THAT number of negative pads   in y-direction = " << that.nPadsNegY() << std::endl <<
                "EXITING." << std::endl;
        exit(1);
    }
}

/*
 * sizeDirectionX(): returns the size (number of interior points + pads) of this in the x-direction
 */
template<class T> int CaslArray2D<T>::sizeDirectionX() const { return _nPadsNegX + _nX + _nPadsPosX; }

/*
 * sizeDirectionY(): returns the size (number of interior points + pads) of this in the y-direction
 */
template<class T> int CaslArray2D<T>::sizeDirectionY() const { return _nPadsNegY + _nY + _nPadsPosY; }

/*
 * nX(): returns the number of interior points in the x-direction
 */
template<class T> int CaslArray2D<T>::nX() const { return _nX; }

/*
 * nY(): returns the number of interior points in the y-direction
 */
template<class T> int CaslArray2D<T>::nY() const { return _nY; }

/*
 * iMin(): returns the smallest index (can be negative) in the x-direction
 */
template<class T> int CaslArray2D<T>::iMin() const { return   1 - _nPadsNegX; }

/*
 * iMax(): returns the largest index (can be negative) in the x-direction
 */
template<class T> int CaslArray2D<T>::iMax() const { return _nX + _nPadsPosX; }

/*
 * jMin(): returns the smallest index (can be negative) in the y-direction
 */
template<class T> int CaslArray2D<T>::jMin() const { return   1 - _nPadsNegY; }

/*
 * jMax(): returns the largest index (can be negative) in the y-direction
 */
template<class T> int CaslArray2D<T>::jMax() const { return _nY + _nPadsPosY; }

/*
 * nPadsNegX(): returns the number of pads in the negative x-direction
 */
template<class T> int  CaslArray2D<T>::nPadsNegX() const { return _nPadsNegX; }

/*
 * nPadsPosX(): returns the number of pads in the positive x-direction
 */
template<class T> int  CaslArray2D<T>::nPadsPosX() const { return _nPadsPosX; }

/*
 * nPadsNegY(): returns the number of pads in the negative y-direction
 */
template<class T> int  CaslArray2D<T>::nPadsNegY() const { return _nPadsNegY; }

/*
 * nPadsPosY(): returns the number of pads in the positive y-direction
 */
template<class T> int  CaslArray2D<T>::nPadsPosY() const { return _nPadsPosY; }

/*
 * maxAbs(): returns the maximum magnitude of the elements in a CaslArray2D
 */
template<class T> T CaslArray2D<T>::maxAbs() const {
    T maxAbs = 0;
    for (int i = 0; i < _nX; ++i) for (int j = 0; j < _nY; ++j) if ( abs(_data[i][j]) > maxAbs ) maxAbs = abs(_data[i][j]);
    return maxAbs;
}


#endif // CASL_ARRAY2D_CPP