#ifndef CASLUNIFORM_SESOLVER_H
#define CASLUNIFORM_SESOLVER_H

#include <vector>
#include "CaslArray2D.h"
#include <cmath>
#include <complex>


using namespace std;



/// Implementation of efficient representation of discrete Poisson equation matrix (2d)
template<class T> class DPMatrix2D {
protected:
    vector<vector<T>> _data;    // Represents grid as [-y, -x, 0, +x, +y]
    int nX{};                   // Number of values on x axis of grid
    int nY{};                   // Number of values on y axis of grid
    int _n{};                   // Length of matrix (nX * nY)
    T _scale;                   // Matrix scale factor when calculating operator/preprocessing

    /// Returns dot product between vectors u and v
    T dot(const CaslArray2D<T>& u, const CaslArray2D<T>& v);
    // Returns dfft
private:
public:
    vector<complex<T>> fftPadded(const vector<complex<T>>& x, bool inverse=false);
    vector<vector<complex<T>>> fftPadded(const vector<vector<complex<T>>>& x, bool inverse=false);
    vector<complex<T>> fft(const vector<complex<T>>& x, bool inverse=false);
    vector<vector<complex<T>>> fft(const vector<vector<complex<T>>>& x, bool inverse=false);

    vector<vector<complex<T>>> fft2D(const CaslArray2D<T>& x);
    CaslArray2D<T> ifft2D(const vector<vector<complex<T>>>& x);

    /// Instantiate DPMatrix of zeros with associated grid dimensions nX * nY
    explicit DPMatrix2D(int nX, int nY);
    /// Instantiate DPMatrix of zeros with associated square grid dimensions m * m
    explicit DPMatrix2D(int m) : DPMatrix2D(m, m) {}



    /// Return length of matrix
    int n() {return _n;}

    /// Return value of matrix at (I, J)
    T operator()(int I, int J) const;

    /// Get value of data representation at [i][j]
    const T& data(int i, int j) const {return _data[i][j];}

    /// Set value of data representation at [i][j]
    T& data(int i, int j) {return _data[i][j];}

    /// Set individual element of matrix (primarily for debug)
    void init(int I, int J, T val);



    /// Fill matrix with a dirichlet boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    /// @param scale data representation divides arguments by scale for better convergence
    void dirichletBoundary(T diag, T offX, T offY, T scale=1);

    /// Fill matrix with a neumann boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    /// @param scale data representation divides arguments by scale for better convergence
    void neumannBoundary(T diag, T offX, T offY, T scale=1);

    /// Fill matrix with a neumann boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param off value of off diagonals associated with adjacent cells
    void neumannBoundary(T diag, T off) { neumannBoundary(diag, off, off); }

    /// Fill matrix with a constant interpolation boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    void constantBoundary(T diag, T offX, T offY, T scale=1);

    /// Fill matrix with a linear interpolation boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    void linearBoundary(T diag, T offX, T offY, T scale=1);


    /// Apply scale factor to inner grid cells
    void preprocessInner(CaslArray2D<T>& b);

    /// Apply scale factor to all grid cells
    void preprocessFull(CaslArray2D<T>& b);

    /// Return the result of applying the associated operator of the matrix on the vector
    CaslArray2D<T> applyOperator(CaslArray2D<T>& x);


    /// Left-hand side matrix multiplication
    /// @param x applies matrix to x vector
    /// @param Ax where to place the resulting grid (overwrites old values)
    virtual void matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax);



    /// ILU(0) decomposition of matrix
    /// @param L where to store resulting lower triangular matrix (overwrites old values)
    /// @param U where to store resulting lower triangular matrix (overwrites old values)
    virtual void ILU(DPMatrix2D<T>& L, DPMatrix2D<T>& U);

    /// Solve L x0 = b for lower triangular L
    virtual void lSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b);

    /// Solve U x0 = b for upper triangular U
    virtual void uSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b);

    /// Solve L U x0 = b for lower triangular L and upper triangular U
    /// @param x0 LHS vector (resulting solution to system)
    /// @param b RHS vector
    /// @param L Lower triangular matrix
    /// @param U upper triangular matrix
    /// @param temp temporary vector used for solving system (overwrites old values)
    void luSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                 DPMatrix2D<T>& L, DPMatrix2D<T>& U, CaslArray2D<T>& temp);



    /// Solves for x in Ax=b using Gauss-Seidel method
    /// @param x0 initial guess for x (changed in place)
    /// @param b target value for Ax
    /// @param niter maximum number of iterations before end
    void gaussSeidel(CaslArray2D<T>& x0, const CaslArray2D<T> &b, int niter);


    /// Solves for x in Ax=b using Conjugate-Gradient method
    /// @param x0 initial guess for x (changed in place)
    /// @param b target value for Ax
    /// @param tol value of residual to end iteration
    /// @param niter maximum number of iterations before end
    void conjGrad(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter= 1000);

    /// Solves for x in Ax=b using Conjugate-Gradient method with ILU(0) preconditioner
    /// @param x0 initial guess for x (changed in place)
    /// @param b target value for Ax
    /// @param L lower triangular matrix part of preconditioner
    /// @param U upper triangular matrix part of preconditioner
    /// @param tol value of residual to end iteration
    /// @param niter maximum number of iterations before end
    void conjGradILU(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                     DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter=1000);

    /// Solves for x in Ax=b using Conjugate-Gradient method.
    /// NULL version solves systems with a ones vector null space by adjusting average value to 0
    /// @param x0 initial guess for x (changed in place)
    /// @param b target value for Ax
    /// @param tol value of residual to end iteration
    /// @param niter maximum number of iterations before end
    void conjGrad_NULL(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter=1000);

    /// Solves for x in Ax=b using Conjugate-Gradient method with ILU(0) preconditioner.
    /// NULL version solves systems with a ones vector null space by adjusting average value to 0
    /// @param x0 initial guess for x (changed in place)
    /// @param b target value for Ax
    /// @param L lower triangular matrix part of preconditioner
    /// @param U upper triangular matrix part of preconditioner
    /// @param tol value of residual to end iteration
    /// @param niter maximum number of iterations before end
    void conjGradILU_NULL(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                          DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter=1000);
};

template <class T> class PeriodicSolver {
private:
    T center, ax, ay, tol{1e-10};
    int nX, nY;
public:
    PeriodicSolver(T center, T ax, T ay, int nX, int nY):
        center(center), ax(ax), ay(ay), nX(nX), nY(nY) {};
    T getTolerance() const {return tol;}
    void setTolerance(T tolerance) {tol = tolerance;}

    vector<complex<T>> fftPadded(const vector<complex<T>>& x, bool inverse=false);
    vector<vector<complex<T>>> fftPadded(const vector<vector<complex<T>>>& x, bool inverse=false);
    vector<complex<T>> fft(const vector<complex<T>>& x, bool inverse=false);
    vector<vector<complex<T>>> fft(const vector<vector<complex<T>>>& x, bool inverse=false);
    vector<vector<complex<T>>> fft2D(const CaslArray2D<T>& x);
    CaslArray2D<T> ifft2D(const vector<vector<complex<T>>>& x);
    CaslArray2D<T> solve(const CaslArray2D<T>& x);
};

template<class T> class DPMatrixPeriodic2D: public DPMatrix2D<T>{
protected:
    using DPMatrix2D<T>::_data, DPMatrix2D<T>::nX, DPMatrix2D<T>::nY, DPMatrix2D<T>::_n, DPMatrix2D<T>::_scale;
public:
    using DPMatrix2D<T>::DPMatrix2D;

    /// Fill matrix with a periodic boundary and symmetrical adjacent grid contributions
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    void periodicBoundary(T diag, T offX, T offY, T scale=1);

    /// Left-hand side matrix multiplication with periodic
    /// @param x applies matrix to x vector
    /// @param Ax where to place the resulting grid (overwrites old values)
    void matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax);

    /// ILU(0) decomposition of matrix
    /// @param L where to store resulting lower triangular matrix (overwrites old values)
    /// @param U where to store resulting lower triangular matrix (overwrites old values)
    void ILU(DPMatrixPeriodic2D<T>& L, DPMatrixPeriodic2D<T>& U);

    /// Solve L x0 = b for lower triangular L
    void lSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b);

    /// Solve U x0 = b for upper triangular U
    void uSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b);
};

template<typename T, class Matrix> class DPMatrixExtended2D: public Matrix {
protected:
    using Matrix::_data, Matrix::nX, Matrix::nY, Matrix::_n, Matrix::dot, Matrix::luSolve;
    vector<vector<T>> _extend;          // Extends grid as [-2x +2x -2y +2y]
    CaslArray2D<T>* preconditioner;     // Row norms of original matrix
public:
    DPMatrixExtended2D(int nX, int nY);
    explicit DPMatrixExtended2D(int m) : DPMatrixExtended2D(m, m) {};
    ~DPMatrixExtended2D() {delete preconditioner;}

    /// Fill matrix with a quadratic boundary condition
    /// @param diag value of diagonals of matrix
    /// @param offX value of off diagonals associated with x-adjacent cells
    /// @param offY value of off diagonals associated with y-adjacent cells
    void quadraticBoundary(T diag, T offX, T offY);


    // Left-hand side matrix multiplication extended a grid square
    /// @param x applies matrix to x vector
    /// @param Ax where to place the resulting grid (overwrites old values)
    void matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax);

    void conjGradQUAD(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter);
    void biconjGradQUAD(CaslArray2D<T>& x0, const CaslArray2D<T> &b, DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter);
};

#endif //CASLUNIFORM_SESOLVER_H
