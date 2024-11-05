//
// Created by Frederic Gibou on 12/28/22.
//

#ifndef CASL_GRID2D_H
#define CASL_GRID2D_H

class CaslGrid2D {
private:
    int    _nX, _nY;                    // number of grid points in x- and y- directions
    double _xMin, _xMax, _yMin, _yMax;  // domain dimensions
    double _dx, _dy;                    // grid spacing in x- and y- directions

public:
    CaslGrid2D(double xMin, double xMax, double yMin, double yMax, int nX, int nY) {
        _xMin = xMin; _xMax = xMax; _nX = nX;
        _yMin = yMin; _yMax = yMax; _nY = nY;
        _dx = ( _xMax - _xMin ) / (double) ( _nX - 1 );
        _dy = ( _yMax - _yMin ) / (double) ( _nY - 1 );
    }

    CaslGrid2D() {
        _xMin = 0; _xMax = 0; _nX = 0;
        _yMin = 0; _yMax = 0; _nY = 0;
        _dx = 0;
        _dy = 0;
    }

    ~CaslGrid2D() = default;

    // nX(): returns the number of grid points in the x-direction
    int nX() const { return _nX; }

    // nY(): returns the number of grid points in the y-direction
    int nY() const { return _nY; }

    // dx(): returns the grid spacing in the x-direction
    double dx() const { return _dx; }

    // dy(): returns the grid spacing in the y-direction
    double dy() const { return _dy; }

    // xMin(): returns the minimum value of grid in the x-direction
    double xMin() const { return _xMin; }

    // xMax(): returns the maximum value of grid in the x-direction
    double xMax() const { return _xMax; }

    // yMin(): returns the minimum value of grid in the x-direction
    double yMin() const { return _yMin; }

    // xMax(): returns the maximum value of grid in the x-direction
    double yMax() const { return _yMax; }

    // x(i): returns the x-coordinate associated with grid point index (i, .)
    double x(int i) const { return _xMin + (double) (i-1) * _dx; }  // start indices at 1 in simulations.

    // y(j): returns the y-coordinate associated with grid point index (., j)
    double y(int j) const { return _yMin + (double) (j-1) * _dy; }  // start indices at 1 in simulations.
};

#endif // CASL_GRID2D_H