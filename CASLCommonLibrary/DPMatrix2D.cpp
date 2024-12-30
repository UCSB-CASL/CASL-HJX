#include "DPMatrix2D.h"

template class DPMatrix2D<double>;
template class DPMatrix2D<float>;
template class PeriodicSolver<double>;
template class PeriodicSolver<float>;
template class DPMatrixPeriodic2D<double>;
template class DPMatrixPeriodic2D<float>;
template class DPMatrixExtended2D<double, DPMatrix2D<double>>;

int mod(int x, int m) {return (x % m + m) % m;}

/// Checks to make sure CASL grid has pads
template<class T>
void assertPads(CaslArray2D<T>& x){
    if (   x.nPadsPosX() == 0
           or x.nPadsNegX() == 0
           or x.nPadsPosY() == 0
           or x.nPadsNegY() == 0)
    {cout << "ERROR: No padding in grid"; exit(1);}
}

template <class T>
T DPMatrix2D<T>::dot(const CaslArray2D<T>& u, const CaslArray2D<T>& v){
    T output = 0;
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
        output += u(i, j) * v(i, j);
    return output;
}

//template<class T>
//vector<complex<T>> DPMatrix2D<T>::fft(const vector<complex<T>>& x, bool inverse){
//    _n = x.size();
//
//    vector<complex<T>> a(_n);
//    vector<complex<T>> b(_n);
//
//    for(int n = 0; n < _n; ++n) {
//        complex<T> temp = exp(complex<T>(0, (n*n) * -M_PI / _n));
//        if (inverse) temp = conj(temp);
//
//        a[n] = x[n] * temp;
//        b[n] = conj(temp);
//    }
//
//    vector<complex<T>> ainv = fftPadded(a);
//    vector<complex<T>> output = fftPadded(b);
//
//    for (int n = 0; n < ainv.size(); ++n)
//        output[n] = ainv[n] * output[n];
//
//    output = fftPadded(output, true);
//    output.resize(_n);
//
//    for (int n = 0; n < _n; ++n)
//        output[n] = conj(b[n]) * output[n];
//
//    if (inverse)
//        for (int n = 0; n < _n; ++n) output[n] /= complex<T>(_n);
//
//    return output;
//}

/*
 * Matrix Implementation
 */

template <class T> DPMatrix2D<T>::DPMatrix2D(int nX, int nY) {
    _scale = 1;
    this->_n = nX * nY;
    this->nX = nX;
    this->nY = nY;

    _data = vector<vector<T>>(_n);
    for (int i = 0; i < _n; ++i) _data[i] = vector<T>(5);
}

template<class T>
T DPMatrix2D<T>::operator()(int I, int J) const{
    // Not especially performant but may be useful for debug

    if      (J == mod(I - nX, _n) )
        return _data[I][0];

    else if (J == (mod(I - 1, nX) + nX * (I / nX)))
        return _data[I][1];

    else if (J == I )
        return _data[I][2];

    else if (J == (mod(I + 1, nX) + nX * (I / nX)))
        return _data[I][3];

    else if (J == ((I + nX) % _n))
        return _data[I][4];

    else return 0;
}

template<class T>
void DPMatrix2D<T>::init(int I, int J, T val){
    // Only works if editing cell in matrix
    if      (J == mod(I - nY, _n))
        _data[I][0] = val;

    else if (J == (mod(I - 1, nX) + nY * (I / nY)))
        _data[I][1] = val;

    else if (J == I)
        _data[I][2] = val;

    else if (J == (mod(I + 1, nX) + nY * (I / nY)))
        _data[I][3] = val;

    else if (J == ((I + nY) % _n))
        _data[I][4] = val;
}

template<class T>
void DPMatrix2D<T>::dirichletBoundary(T diag, T offX, T offY, T scale){
    // Dirichlet boundary condition
    _scale = scale;

    diag /= scale;
    offX /= scale;
    offY /= scale;

    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i){
        // Boundary handling
        if(    i == 1 or i == nX
            or j == 1 or j == nY)
        {
            _data[id][2] = 1;
            id += 1;
            continue;
        }

        // Inner cells
        _data[id][0] = offY;
        _data[id][1] = offX;
        _data[id][2] = diag;
        _data[id][3] = offX;
        _data[id][4] = offY;

        id += 1;
    }
}

template<class T>
void DPMatrix2D<T>::neumannBoundary(T diag, T offX, T offY, T scale){
    // Neumann boundary condition
    _scale = scale;

    diag /= scale;
    offX /= scale;
    offY /= scale;

    const double offX2 = offX + offX;
    const double offY2 = offY + offY;

    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
        _data[id][2] = diag;

        // Determine whether X boundary
        if      (i == 1 ) {_data[id][1] = 0;     _data[id][3] = offX2;} // Left X
        else if (i == nX) {_data[id][1] = offX2; _data[id][3] = 0;    } // Right X
        else              {_data[id][1] = offX;  _data[id][3] = offX; } // Inner

        // Determine whether Y boundary
        if      (j == 1 ) {_data[id][0] = 0;     _data[id][4] = offY2;} // Left Y
        else if (j == nY) {_data[id][0] = offY2; _data[id][4] = 0;    } // Right Y
        else              {_data[id][0] = offY;  _data[id][4] = offY; } // Inner

        id += 1;
    }
}

template<class T>
void DPMatrix2D<T>::constantBoundary(T diag, T offX, T offY, T scale){
    // Constant boundary condition
    _scale = scale;

    diag /= scale;
    offX /= scale;
    offY /= scale;

    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
            _data[id][2] = diag;

            // Determine whether X boundary
            if      (i == 1 ) {_data[id][3] = offX; _data[id][2] += offX; } // Left X
            else if (i == nX) {_data[id][1] = offX; _data[id][2] += offX; } // Right X
            else              {_data[id][1] = offX;  _data[id][3] = offX; } // Inner

            // Determine whether Y boundary
            if      (j == 1 ) {_data[id][4] = offY; _data[id][2] += offY; } // Left Y
            else if (j == nY) {_data[id][0] = offY; _data[id][2] += offY; } // Right Y
            else              {_data[id][0] = offY; _data[id][4] = offY; } // Inner

            id += 1;
        }
}

template<class T>
void DPMatrix2D<T>::linearBoundary(T diag, T offX, T offY, T scale){
    // Linear boundary condition
    _scale = scale;

    diag /= scale;
    offX /= scale;
    offY /= scale;

    const double offX2 = offX + offX;
    const double offY2 = offY + offY;

    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
        _data[id][2] = diag;

        // Determine whether X boundary
        if (i == 1 || i == nX) {
            _data[id][1] = 0; _data[id][3] = 0;
            _data[id][2] += offX2;
        }
        else {_data[id][1] = offX;  _data[id][3] = offX; } // Inner


        // Determine whether Y boundary
        if (j == 1 || j == nY) {
            _data[id][0] = 0; _data[id][4] = 0;
            _data[id][2] += offY2;
        }
        else { _data[id][0] = offY;  _data[id][4] = offY; } // Inner

        id += 1;
    }
}

template<class T>
void DPMatrix2D<T>::preprocessInner(CaslArray2D<T>& b){
    for (int i = 2; i < b.nX(); ++i) for (int j = 2; j < b.nY(); ++j)
        b(i, j) /= _scale;
}

template<class T>
void DPMatrix2D<T>::preprocessFull(CaslArray2D<T>& b){
    for (int i = 1; i <= b.nX(); ++i) for (int j = 1; j <= b.nY(); ++j)
        b(i, j) /= _scale;
}

template<class T>
CaslArray2D<T> DPMatrix2D<T>::applyOperator(CaslArray2D<T>& x){
    // Apply operator on vector
    CaslArray2D<T> output(nX, nY);
    matMult(x, output);

    // Scale output
    if (_scale == 1) return output;

    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            output(i, j) *= _scale;

    return output;
}

template<class T>
void DPMatrix2D<T>::matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax) {
    assertPads(x);

    int id = 0;

    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
        Ax(i, j) =_data[id][0] * x(i, j - 1)
                + _data[id][1] * x(i - 1, j)
                + _data[id][2] * x(i, j)
                + _data[id][3] * x(i + 1, j)
                + _data[id][4] * x(i, j + 1);
        id += 1;
    }
}

template <class T>
void DPMatrix2D<T>::ILU(DPMatrix2D<T>& L, DPMatrix2D<T>& U) {
    int kx, ky;
    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
        L.data(id, 2) = 1;                                  // L_ii = 1
        U.data(id, 2) = _data[id][2];                       // U_ii = a_ii

        if (i > 1) {
            kx = id - 1;                                    // -x adjacent index
            U.data(kx, 3) = _data[kx][3];                   // U_ij = a_ij
            L.data(id, 1) = _data[id][1] / (U.data(kx, 2) ? U.data(kx, 2) : 1);     // L_ik = a_ik / u_kk
            U.data(id, 2) -= L.data(id, 1) * _data[kx][3];    // U_ii = a_ii - l_ik*u_ki - ...
        }

        if (j > 1) {
            ky = id - nX;                                   // -y adjacent index
            U.data(ky, 4) = _data[ky][4];                   // U_ij = a_ij
            L.data(id, 0) = _data[id][0] / (U.data(ky, 2) ? U.data(ky, 2) : 1);     // L_ik = a_ik / u_kk
            U.data(id, 2) -= L.data(id, 0) * _data[ky][4];    // U_ii = a_ii - l_ik*u_ki - ...
        }

        id += 1;
    }
}


template <class T>
void DPMatrix2D<T>::lSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b){
    // Solve Lx = b for lower triangular L
    int id = 0;

    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
        x0(i, j) = (b(i,j) - x0(i, j-1) * _data[id][0]
                           - x0(i-1, j) * _data[id][1]) / _data[id][2];
        id += 1;
    }
}

template <class T>
void DPMatrix2D<T>::uSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b){
    // Solve Ux = b for upper triangular U
    int id = _n;

    // Center rows
    for (int j = nY; j >= 1; --j) for (int i = nX; i >= 1; --i) {
        id -= 1;
        x0(i, j) = (b(i,j) - x0(i, j+1) * _data[id][4]
                           - x0(i+1, j) * _data[id][3]) / _data[id][2];
    }
}

template <class T>
void DPMatrix2D<T>::luSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                            DPMatrix2D<T>& L, DPMatrix2D<T>& U, CaslArray2D<T>& temp){
    assertPads(temp); assertPads(x0);
    L.lSolve(temp, b);      // temp = L^-1 b
    U.uSolve(x0, temp);     // x0 = U^-1 L^-1 b
    // LU x0 = b
}


template<class T>
void DPMatrix2D<T>::gaussSeidel(CaslArray2D<T>& x0, const CaslArray2D<T> &b, int niter){
    int I, J;
    for(int n = 0; n < niter; ++n){
        for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i)  {
                x0(i, j) = b(i, j);
                I = i-1 + (j-1)*nX;
                for (int l = 1; l <= nY; ++l) for (int k = 1; k <= nX; ++k) {
                        J = k - 1 + (l - 1) * nX;
                        x0(i, j) -= (*this)(I, J) * x0(k, l) * (I != J);
                    }

                x0(i, j) /= _data[i][2];
            }

    }
}


template<class T>
void DPMatrix2D<T>::conjGrad(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter){
    // Variable initialization
    tol = tol*tol;

    CaslArray2D<T> r(nX, nY);
    matMult(x0, r);
    CaslArray2D<T> p(nX, nY, 1);
    CaslArray2D<T> Ap(nX, nY);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j){
        r(i, j) = b(i, j) - r(i, j);
        p(i, j) = r(i, j);
    }

    T a;
    T B;
    T rTr = dot(r,r);
    // Main loop
    for(int k = 0; k < niter; ++k){
        matMult(p, Ap);

        a = rTr / dot(p, Ap);

        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)  {
            x0(i, j) += a * p(i, j);
            r(i, j)  -= a * Ap(i, j);
        }

        rTr = dot(r,r);
        if (rTr < tol) break;

        B = -(dot(r, Ap) / dot(p, Ap));
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            p(i, j) = r(i, j) + B * p(i, j);
    }
}

template<class T>
void DPMatrix2D<T>::conjGradILU(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                                DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter){
    // Variable initialization
    tol = tol*tol;
    CaslArray2D<T> z(nX, nY, 1, 1), p(nX, nY, 1, 1), Ap(nX, nY, 1, 1);
    T rTz_k, rTz_k1, a, B;

    CaslArray2D<T> r(x0.nX(), x0.nY());
    matMult(x0, r);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
        r(i, j) = b(i, j) - r(i, j);     // r = b - Ax0

    luSolve(z, r, L, U, Ap);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
        p(i, j) = z(i, j);
    rTz_k1 = dot(r, z);
    // Main loop
    for(int k = 0; k < niter; ++k){
        matMult(p, Ap);
        a = rTz_k1 / dot(p, Ap);
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j){
            x0(i, j) += a * p(i, j);
            r(i, j)  -= a * Ap(i, j);
        }
        if (dot(r, r) < tol) break;
        luSolve(z, r, L, U, Ap);
        rTz_k = rTz_k1;
        rTz_k1 = dot(r, z);
        B = rTz_k1 / rTz_k;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            p(i, j) = B * p(i, j) + z(i, j);
    }
}

template<class T>
void DPMatrix2D<T>::conjGrad_NULL(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter){
    // Variable initialization
    tol = tol*tol;
    T a, B;
    double sum, average;
    CaslArray2D<T> p(nX, nY, 1, 1), Ap(nX, nY);

    CaslArray2D<T> r(nX, nY);
    matMult(x0, r);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j){
        r(i, j) = b(i, j) - r(i, j);
        p(i, j) = r(i, j);
    }

    T rTr = dot(r,r);
    // Main loop
    for(int k = 0; k < niter; ++k){
        matMult(p, Ap);

        a = rTr / dot(p, Ap);
        sum = 0;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            x0(i, j) += a * p(i, j);
            sum += x0(i, j);
            r(i, j)  -= a * Ap(i, j);
        }
        average = sum / _n;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            x0(i, j) -= average;
        rTr = dot(r,r);

//        cout << rTr << ", ";
        if (rTr < tol) break;

        B = -(dot(r, Ap) / dot(p, Ap));
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            p(i, j) = r(i, j) + B * p(i, j);
    }
}

template<class T>
void DPMatrix2D<T>::conjGradILU_NULL(CaslArray2D<T>& x0, const CaslArray2D<T> &b,
                                     DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter){
    // Variable initialization
    tol = tol*tol;
    CaslArray2D<T> z(nX, nY, 1, 1), p(nX, nY, 1, 1), Ap(nX, nY, 1, 1);
    T rTz_k, rTz_k1, a, B;

    CaslArray2D<T> r(nX, nY);
    matMult(x0, r);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
        r(i, j) = b(i, j) - r(i, j);     // r = b - Ax0

    luSolve(z, r, L, U, Ap);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
        p(i, j) = z(i, j);
    rTz_k1 = dot(r, z);
    double sum, average;
    // Main loop
    for(int k = 0; k < niter; ++k) {
        matMult(p, Ap);
        a = rTz_k1 / dot(p, Ap);
        sum = 0;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
                x0(i, j) += a * p(i, j);
                sum += x0(i, j);
                r(i, j) -= a * Ap(i, j);
            }
        average = sum / _n;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
                x0(i, j) -= average;
        if (dot(r, r) < tol) break;
        luSolve(z, r, L, U, Ap);
        rTz_k = rTz_k1;
        rTz_k1 = dot(r, z);
        B = rTz_k1 / rTz_k;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
                p(i, j) = B * p(i, j) + z(i, j);
    }
}


/// Periodic

int bit_reverse(int I, int N) {
    int rev = 0;

    for (int i = 0; i < N-1; ++i) {
        rev |= I & 1;
        rev <<= 1;
        I >>= 1;
    }

    return rev | I & 1;
}

template <class T>
vector<complex<T>> PeriodicSolver<T>::fftPadded(const vector<complex<T>>& x, bool inverse){
    int N = 2;
    int N_bits = 1;
    int N0 = x.size() - 1;
    while ((N0 = N0 >> 1)) {
        N = N << 1;
        N_bits++;
    }
    vector<complex<T>> x_hat(N);
    for(int i = 0; i < x.size(); ++i)
        x_hat[bit_reverse(i, N_bits)] = x[i];

    int n = 2;
    for (int i = 1; i <= N_bits; ++i) {
        complex<T> wi = exp(complex<T>(0, 2*M_PI / n));
        if (inverse) wi = conj(wi);

        for (int k = 0; k < N; k += n){
//            cout << i <<" " << k << ":\n";
//            for(auto it: inv) cout << it << ", "; cout << endl;

            complex<T> w = 1;
            for (int j = 0; j < n/2; ++j) {
                complex<T> second = w * x_hat[k + j + n / 2];
                x_hat[k + j + n / 2] = x_hat[k + j] - second;
                x_hat[k + j] += second;

                w = w * wi;
            }
        }
        n <<= 1;
    }

    if (inverse)
        for (int i = 0; i < x_hat.size(); ++i)
            x_hat[i] /= complex<T>(N);

    return x_hat;
}

template <class T>
vector<vector<complex<T>>> PeriodicSolver<T>::fftPadded(const vector<vector<complex<T>>>& x, bool inverse){
    int N = 2;
    int N_bits = 1;
    int N0 = x.size() - 1;
    while ((N0 = N0 >> 1)) {
        N = N << 1;
        N_bits++;
    }

    int M = x[0].size();
    vector<vector<complex<T>>> x_hat(N, vector<complex<T>>(M));
    for(int i = 0; i < x.size(); ++i)
        x_hat[bit_reverse(i, N_bits)] = x[i];

    int n = 2;
    for (int i = 1; i <= N_bits; ++i) {
        complex<T> wi = exp(complex<T>(0, 2*M_PI / n));
        if (inverse) wi = conj(wi);

        for (int k = 0; k < N; k += n){
            complex<T> w = 1;
            for (int j = 0; j < n/2; ++j) {
                for (int p = 0; p < M; ++p) {
                    complex<T> second = w * x_hat[k + j + n / 2][p];
                    x_hat[k + j + n / 2][p] = x_hat[k + j][p] - second;
                    x_hat[k + j][p] += second;
                }
                w = w * wi;
            }
        }
        n <<= 1;
    }

    if (inverse)
        for (int i = 0; i < x_hat.size(); ++i) for (int j = 0; j < M; ++j)
                x_hat[i][j] /= complex<T>(N);

    return x_hat;
}

template <class T>
vector<complex<T>> PeriodicSolver<T>::fft(const vector<complex<T>>& x, bool inverse){
    // Calculate zero padding
    int M = 2;
    int N_bits = 1;
    int N0 = 2*x.size() - 1;
    while ((N0 = N0 >> 1)) {
        M = M << 1;
        N_bits++;
    }

    int N = x.size();
    // Preprocess convolution
    vector<complex<T>> a(M), b(M);
    for(int n = 0; n < N; ++n) {
        complex<T> factor = exp(complex<T>(0,-(n*n) * M_PI / N));
        if (inverse) factor = conj(factor);

        a[n] = x[n] * factor;
        b[n] = conj(factor);
        if (n > 0) b[M - n] = b[n];
    }

    // Transform
    vector<complex<T>> ainv = fftPadded(a);
    vector<complex<T>> output = fftPadded(b);

    // Multiply
    for(int i = 0; i < M; ++i)
        output[i] *= ainv[i];

    // Invert
    output = fftPadded(output, true);

    // Multiply factor
    for (int i = 0; i < N; ++i)
        output[i] *= conj(b[i]);

    if (inverse) for (int i = 0; i < N; ++i)
            output[i] /= x.size();

    // Postprocess
    output.resize(N); // Fix length

    return output;
}

template <class T>
vector<vector<complex<T>>> PeriodicSolver<T>::fft(const vector<vector<complex<T>>>& x, bool inverse){
    // Calculate zero padding
    int M = 2;
    int M_bits = 1;
    int N0 = 2*x.size() - 1;
    while ((N0 = N0 >> 1)) {
        M = M << 1;
        M_bits++;
    }
    int N = x.size();
    int K = x[0].size();
    // Preprocess convolution
    vector<vector<complex<T>>> a(M, vector<complex<T>>(K)), b(M, vector<complex<T>>(K));
    for(int n = 0; n < N; ++n) {

        complex<T> factor = exp(complex<T>(0,-(n*n) * M_PI / N));
        if (inverse) factor = conj(factor);

        for (int k = 0; k < K; ++k) {
            a[n][k] = x[n][k] * factor;
            b[n][k] = conj(factor);
            if (n > 0) b[M - n][k] = b[n][k];
        }
    }
    // Transform
    vector<vector<complex<T>>> ainv = fftPadded(a);
    vector<vector<complex<T>>> output = fftPadded(b);

    // Multiply
    for(int i = 0; i < M; ++i) for(int k = 0; k < K; ++k)
            output[i][k] *= ainv[i][k];

    // Invert
    output = fftPadded(output, true);

    // Multiply factor
    for(int i = 0; i < M; ++i) for(int k = 0; k < K; ++k)
            output[i][k] *= conj(b[i][k]);

    if (inverse) for(int i = 0; i < M; ++i) for(int k = 0; k < K; ++k)
                output[i][k] /= N;

    // Postprocess
    output.resize(N); // Fix length

    return output;
}


template <class T>
vector<vector<complex<T>>> PeriodicSolver<T>::fft2D(const CaslArray2D<T>& x) {
    vector<vector<complex<T>>> x_hat(nX, vector<complex<T>>(nY));
    for (int j = 1; j <= x.nY(); ++j)
        for(int i = 1; i <= x.nX(); ++i)
            x_hat[j-1][i-1] = x(i, j);


    x_hat = fft(x_hat);

    for (int i = 0; i < x_hat.size(); ++i)
        x_hat[i] = fft(x_hat[i]);

    return x_hat;
}

template <class T>
CaslArray2D<T> PeriodicSolver<T>::ifft2D(const vector<vector<complex<T>>>& x){
    vector<vector<complex<T>>> x_hat(nX, vector<complex<T>>(nY));
    for (int j = 0; j < x.size(); ++j)
        for(int i = 0; i < x[0].size(); ++i)
            x_hat[j][i] = x[j][i];

    for (int i = 0; i < x_hat.size(); ++i)
        x_hat[i] = fft(x_hat[i], true);

    x_hat = fft(x_hat, true);

    CaslArray2D<T> output(nX, nY);

    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i)
            output(i, j) = real(x_hat[j-1][i-1]);

    return output;
}

template <class T>
CaslArray2D<T> PeriodicSolver<T>::solve(const CaslArray2D<T>& x){
    vector<vector<complex<T>>> x_hat = fft2D(x);
    int N = x_hat.size();
    int M = x_hat[0].size();

    for(int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j){
            double p = center
                    + 2 * ay*cos(2 * M_PI * i / N)
                    + 2 * ax*cos(2 * M_PI * j / M);

            if (abs(p) < tol) x_hat[i][j] = 0;
            else x_hat[i][j] /= p;
        };


    return ifft2D(x_hat);
}

/// Periodic Matrix

template<class T>
void DPMatrixPeriodic2D<T>::periodicBoundary(T diag, T offX, T offY, T scale){
    _scale = scale;

    diag /= scale;
    offX /= scale;
    offY /= scale;

    for (int I = 0; I < _n; ++I){
        _data[I][0] = offY; _data[I][4] = offY;
        _data[I][1] = offX; _data[I][3] = offX;
        _data[I][2] = diag;
    }
}

template<class T>
void DPMatrixPeriodic2D<T>::matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax) {
    const int a = _n-1; const int b = nX-1;
    // Top corners
    Ax(1, 1) = _data[0][0] * x(1, nY) + _data[0][4] * x(1, 2)
             + _data[0][1] * x(nX, 1) + _data[0][3] * x(2, 1) + _data[0][2] * x(1, 1);


    Ax(nX, 1) = _data[b][0] * x(nX, nY) + _data[b][4] * x(nX, 2)
              + _data[b][1] * x(b , 1 ) + _data[b][3] * x(1 , 1) + _data[b][2] * x(nX, 1);

    int id = 1; int id2 = _n - nX + 1;
    // First and Last Rows
    for (int i = 2; i < nX; ++i) {
        Ax(i, 1) = _data[id][0] * x(i, nY) + _data[id][4] * x(i, 2)
                 + _data[id][1] * x(i-1, 1) + _data[id][3] * x(i+1, 1) + _data[id][2] * x(i, 1);

        Ax(i, nY) = _data[id2][0] * x(i, nY-1) + _data[id2][4] * x(i, 1)
                  + _data[id2][1] * x(i-1, nY) + _data[id2][3] * x(i+1, nY) + _data[id2][2] * x(i, nY);

        id += 1;
        id2 += 1;
    }
    id += 1;

    // Middle Rows
    for (int j = 2; j < nY; ++j){
        // Left side
        Ax(1, j) = _data[id][0] * x(1 , j-1) + _data[id][4] * x(1, j+1)
                 + _data[id][1] * x(nX, j  ) + _data[id][3] * x(2, j  ) + _data[id][2] * x(1, j);

        id += 1;

        // Inside
        for (int i = 2; i < nX; ++i) {
            Ax(i, j) = _data[id][0] * x(i, j-1) + _data[id][4] * x(i, j+1)
                     + _data[id][1] * x(i-1, j) + _data[id][3] * x(i+1, j) + _data[id][2] * x(i, j);
            id += 1;
        }

        // Right side
        Ax(nX, j) = _data[id][0] * x(nX, j-1) + _data[id][4] * x(nX, j+1)
                  + _data[id][1] * x(b , j  ) + _data[id][3] * x(1 , j  ) + _data[id][2] * x(nX, j);
        id += 1;
    }

    // Bottom Corners
    Ax(1, nY) = _data[id][0] * x(1, nY-1) + _data[id][4] * x(1, 1)
              + _data[id][1] * x(nX, nY ) + _data[id][3] * x(2, nY) + _data[id][2] * x(1, nY);

    Ax(nX, nY) = _data[a][0] * x(nX, nY-1) + _data[a][4] * x(nX, 1)
               + _data[a][1] * x(b , nY  ) + _data[a][3] * x(1, nY) + _data[a][2] * x(nX, nY);
}


template <class T>
void DPMatrixPeriodic2D<T>::ILU(DPMatrixPeriodic2D<T>& L, DPMatrixPeriodic2D<T>& U) {
    int kx, ky;
    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
            L.data(id, 2) = 1;                                  // L_ii = 1
            U.data(id, 2) = _data[id][2];                       // U_ii = a_ii

            if (i > 1) {
                kx = id - 1;                                    // -x adjacent index (kx < id)
                U.data(kx, 3) = _data[kx][3];                   // U_ij = a_ij
                L.data(id, 1) = _data[id][1] / U.data(kx, 2);     // L_ik = a_ik / u_kk
                U.data(id, 2) -= L.data(id, 1) * _data[kx][3];    // U_ii = a_ii - l_ik*u_ki - ...
            }
            if (i == nX) {
                kx = id + 1 - nX;                               // +x adjacent index (kx < id)
                U.data(kx, 1) = _data[kx][1];                   // U_ij = a_ij
                L.data(id, 3) = _data[id][3] / U.data(kx, 2);     // L_ik = a_ik / u_kk
                U.data(id, 2) -= L.data(id, 3) * _data[kx][1];    // U_ii = a_ii - l_ik*u_ki - ...
            }

            if (j > 1) {
                ky = id - nX;                                   // -y adjacent index (ky < id)
                U.data(ky, 4) = _data[ky][4];                   // U_ij = a_ij
                L.data(id, 0) = _data[id][0] / U.data(ky, 2);     // L_ik = a_ik / u_kk
                U.data(id, 2) -= L.data(id, 0) * _data[ky][4];    // U_ii = a_ii - l_ik*u_ki - ...
            }
            if (j == nY) {
                ky = id + nY - _n;                              // +y adjacent index (ky < id)
                U.data(ky, 0) = _data[ky][0];                   // U_ij = a_ij
                L.data(id, 4) = _data[id][4] / U.data(ky, 2);     // L_ik = a_ik / u_kk
                U.data(id, 2) -= L.data(id, 4) * _data[ky][0];    // U_ii = a_ii - l_ik*u_ki - ...
            }

            id += 1;
        }
}


template <class T>
void DPMatrixPeriodic2D<T>::lSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b){
    // Solve Lx = b for lower triangular L
    int id = 0;

    for (int j = 1; j < nY; ++j) {
        for (int i = 1; i < nX; ++i) {
            x0(i, j) = (b(i, j) - x0(i, j-1) * _data[id][0]
                                - x0(i-1, j) * _data[id][1]) / _data[id][2];
            id += 1;
        }
        // X = nX
        x0(nX, j) = (b(nX, j) - x0(nX, j-1) * _data[id][0]
                   - x0(1, j) * _data[id][3] - x0(nX-1, j) * _data[id][1]) / _data[id][2];
        id += 1;
    }

    // nY boundary
    for (int i = 1; i < nX; ++i) {
        x0(i, nY) = (b(i, nY)
                    - x0(i, nY-1) * _data[id][0] - x0(i, 1) * _data[id][4]
                    - x0(i-1, nY) * _data[id][1]) / _data[id][2];
        id += 1;
    }

    // Far corner
    x0(nX, nY) = (b(nX, nY)
                 - x0(nX, nY-1) * _data[id][0] - x0(nX, 1) * _data[id][4]
                 - x0(nX-1, nY) * _data[id][1] - x0(1, nY) * _data[id][3]) / _data[id][2];

}

template <class T>
void DPMatrixPeriodic2D<T>::uSolve(CaslArray2D<T>& x0, const CaslArray2D<T> &b){
    // Solve Ux = b for upper triangular U
    int id = _n;

    // Center rows
    for (int j = nY; j > 1; --j) {
        for (int i = nX; i > 1; --i) {
            id -= 1;
            x0(i, j) = (b(i, j) - x0(i, j + 1) * _data[id][4]
                                - x0(i + 1, j) * _data[id][3]) / _data[id][2];
        }
        // X = 1
        id -= 1;
        x0(1, j) = (b(1, j) - x0(1, j+1) * _data[id][4]
                  - x0(2, j) * _data[id][3] - x0(nX, j) * _data[id][0]) / _data[id][2];
    }

    // Y = 1 boundary
    for (int i = nX; i > 1; --i) {
        id -= 1;
        x0(i, 1) = (b(i, 1)
                      - x0(i, nY) * _data[id][0] - x0(i, 2) * _data[id][4]
                      - x0(i+1, 1) * _data[id][3]) / _data[id][2];
    }

    // Far corner
    id -= 1;
    x0(1, 1) = (b(1, 1)
                   - x0(1, nY) * _data[id][0] - x0(1, 2) * _data[id][4]
                   - x0(nX, 1) * _data[id][1] - x0(2, 1) * _data[id][3]) / _data[id][2];
}


/// Extended

template<typename T, class Matrix>
DPMatrixExtended2D<T, Matrix>::DPMatrixExtended2D(int nX, int nY): Matrix(nX, nY) {
    _extend = vector<vector<T>>(_n);
    for (int i = 0; i < _n; ++i) _extend[i] = vector<T>(4);
    preconditioner = new CaslArray2D<T>(nX, nY);
}

template<typename T, class Matrix>
void DPMatrixExtended2D<T, Matrix>::quadraticBoundary(T diag, T offX, T offY) {
    // Quadratic boundary condition

    const T offX2 = offX + offX; const T offX3 = offX2 + offX;
    const T offY2 = offY + offY; const T offY3 = offY2 + offY;

    int id = 0;
    for (int j = 1; j <= nY; ++j) for (int i = 1; i <= nX; ++i) {
            bool edgeX = true, edgeY = true;

            _data[id][2] = diag;

            // Determine whether X boundary
            // Left X
            if      (i == 1 ) {_data[id][1] = 0;      _data[id][3] = -offX2;
                               _data[id][2] += offX3; _extend[id][1] = offX;}
            // Right X
            else if (i == nX) {_data[id][1] = -offX2; _data[id][3] = 0;
                               _data[id][2] += offX3; _extend[id][0] = offX;}
            // Inner
            else              {_data[id][1] = offX;  _data[id][3] = offX; edgeX = false;}


            // Determine whether Y boundary
            // Left Y
            if      (j == 1 ) {_data[id][0] = 0;      _data[id][4] = -offY2;
                               _data[id][2] += offY3; _extend[id][3] = offY;}
            // Right Y
            else if (j == nY) {_data[id][0] = -offY2; _data[id][4] = 0;
                               _data[id][2] += offY3; _extend[id][2] = offY;}
            // Inner
            else              {_data[id][0] = offY;  _data[id][4] = offY; edgeY = false;}


            // Corner handle
            if (edgeX && edgeY) {
                for (int I = 0; I < 5; ++I) _data[id][I] = 0;
                for (int I = 0; I < 4; ++I) _extend[id][I] = 0;
                _data[id][2] = 1;
            }

            id += 1;
        }
}


template<typename T, class Matrix>
void DPMatrixExtended2D<T, Matrix>::matMult(CaslArray2D<T>& x, CaslArray2D<T>& Ax) {
    Matrix::matMult(x, Ax);

    int id2 = 0;
    int id = nX;

    // Top Left Corners
    Ax(1, 1) += _extend[id2][0] * x(nX-1, 1   ) + _extend[id2][1] * x(3, 1)
              + _extend[id2][2] * x(1   , nY-1) + _extend[id2][3] * x(1, 3);
    id2++;
    Ax(2, 1) += _extend[id2][0] * x(nX, 1   ) + _extend[id2][1] * x(4, 1)
              + _extend[id2][2] * x(2 , nY-1) + _extend[id2][3] * x(2, 3);
    id2++;

    Ax(1, 2) += _extend[id][0] * x(nX-1, 2 ) + _extend[id][1] * x(3, 2)
              + _extend[id][2] * x(1   , nY) + _extend[id][3] * x(1, 4);
    id++;
    Ax(2, 2) += _extend[id][0] * x(nX, 2 ) + _extend[id][1] * x(4, 2)
              + _extend[id][2] * x(2 , nY) + _extend[id][3] * x(2, 4);
    id++;


    // Top Rows
    for (int i = 3; i <= nX-2; ++i){
        Ax(i, 1) += _extend[id2][0] * x(i-2, 1   ) + _extend[id2][1] * x(i+2, 1)
                  + _extend[id2][2] * x(i  , nY-1) + _extend[id2][3] * x(i  , 3);
        id2++;
        Ax(i, 2) += _extend[id][0] * x(i-2, 2 ) + _extend[id][1] * x(i+2, 2)
                  + _extend[id][2] * x(i  , nY) + _extend[id][3] * x(i  , 4);
        id++;
    }

    // Top Right Corners
    Ax(nX-1, 1) += _extend[id2][0] * x(nX-3, 1   ) + _extend[id2][1] * x(1   , 1)
                 + _extend[id2][2] * x(nX-1, nY-1) + _extend[id2][3] * x(nX-1, 3);
    id2++;
    Ax(nX, 1) += _extend[id2][0] * x(nX-2, 1   ) + _extend[id2][1] * x(2 , 1)
               + _extend[id2][2] * x(nX  , nY-1) + _extend[id2][3] * x(nX, 3);

    Ax(nX-1, 2) += _extend[id][0] * x(nX-3, 2 ) + _extend[id][1] * x(1   , 2)
                 + _extend[id][2] * x(nX-1, nY) + _extend[id][3] * x(nX-1, 4);
    id++;
    Ax(nX  , 2) += _extend[id][0] * x(nX-2, 2 ) + _extend[id][1] * x(2 , 2)
                 + _extend[id][2] * x(nX  , nY) + _extend[id][3] * x(nX, 4);
    id++;

    for (int j = 3; j <= nY - 2; ++j) {
        // Left Boundary
        Ax(1, j) += _extend[id][0] * x(nX-1, j  ) + _extend[id][1] * x(3, j  )
                  + _extend[id][2] * x(1   , j-2) + _extend[id][3] * x(1, j+2);
        id++;
        Ax(2, j) += _extend[id][0] * x(nX  , j  ) + _extend[id][1] * x(4, j  )
                  + _extend[id][2] * x(2   , j-2) + _extend[id][3] * x(2, j+2);
        id++;

        // Inner
        for (int i = 3; i <= nX - 2; ++i) {
            Ax(i, j) += _extend[id][0] * x(i-2, j  ) + _extend[id][1] * x(i+2, j  )
                      + _extend[id][2] * x(i  , j-2) + _extend[id][3] * x(i  , j+2);
            id++;
        }

        // Right Boundary
        Ax(nX-1, j) += _extend[id][0] * x(nX-3, j  ) + _extend[id][1] * x(1   , j  )
                     + _extend[id][2] * x(nX-1, j-2) + _extend[id][3] * x(nX-1, j+2);
        id++;
        Ax(nX  , j) += _extend[id][0] * x(nX-2, j  ) + _extend[id][1] * x(2 , j  )
                     + _extend[id][2] * x(nX  , j-2) + _extend[id][3] * x(nX, j+2);
        id++;

    }

    // Bottom Left Corners
    id2 = id + nX;

    Ax(1, nY-1) += _extend[id][0] * x(nX-1, nY-1) + _extend[id][1] * x(3, nY-1)
                 + _extend[id][2] * x(1   , nY-3) + _extend[id][3] * x(1, 1    );
    id++;
    Ax(2, nY-1) += _extend[id][0] * x(nX, nY-1) + _extend[id][1] * x(4, nY-1)
                 + _extend[id][2] * x(1 , nY-3) + _extend[id][3] * x(2, 1    );
    id++;

    Ax(1, nY  ) += _extend[id2][0] * x(nX-1, nY  ) + _extend[id2][1] * x(3, nY)
                 + _extend[id2][2] * x(1   , nY-2) + _extend[id2][3] * x(1, 2 );
    id2++;
    Ax(2, nY  ) += _extend[id2][0] * x(nX, nY  ) + _extend[id2][1] * x(4, nY)
                 + _extend[id2][2] * x(2 , nY-2) + _extend[id2][3] * x(2, 2 );
    id2++;

    // Bottom Rows
    for (int i = 3; i <= nX-2; ++i){
        Ax(i, nY-1) += _extend[id][0] * x(i-2, nY-1) + _extend[id][1] * x(i+2, nY-1)
                     + _extend[id][2] * x(i  , nY-3) + _extend[id][3] * x(i  , 1    );
        id++;
        Ax(i, nY  ) += _extend[id2][0] * x(i-2, nY  ) + _extend[id2][1] * x(i+2, nY)
                     + _extend[id2][2] * x(i  , nY-2) + _extend[id2][3] * x(i  , 2 );
        id2++;
    }


    // Bottom Right Corners
    Ax(nX-1, nY-1) += _extend[id][0] * x(nX-3, nY-1) + _extend[id][1] * x(1   , nY-1)
                    + _extend[id][2] * x(nX-1, nY-3) + _extend[id][3] * x(nX-1, 1    );
    id++;
    Ax(nX  , nY-1) += _extend[id][0] * x(nX-2, nY-1) + _extend[id][1] * x(2 , nY-1)
                    + _extend[id][2] * x(nX  , nY-3) + _extend[id][3] * x(nX, 1    );

    Ax(nX-1, nY  ) += _extend[id2][0] * x(nX-3, nY  ) + _extend[id2][1] * x(1   , nY)
                    + _extend[id2][2] * x(nX-1, nY-2) + _extend[id2][3] * x(nX-1, 2 );
    id2++;
    Ax(nX  , nY  ) += _extend[id2][0] * x(nX-2, nY  ) + _extend[id2][1] * x(2 , nY)
                    + _extend[id2][2] * x(nX  , nY-2) + _extend[id2][3] * x(nX, 2 );


}

template<typename T, class Matrix>
void DPMatrixExtended2D<T, Matrix>::conjGradQUAD(CaslArray2D<T>& x0, const CaslArray2D<T> &b, T tol, int niter){
    // Null space vectors preprocess
    CaslArray2D<T> xNull(nX, nY), yNull(nX, nY), xyNull(nX, nY), x2y2Null(nX, nY);
    T xNorm = nX * (nX*nX - 1) / 12.0;
    T yNorm = nY * (nY*nY - 1) / 12.0;
    T xyNorm = xNorm * yNorm;
    xNorm *= nY; yNorm *= nX;
    T x2y2Norm = (0.15*nX*nX - 0.35)*(xNorm + yNorm) - 2*xyNorm;

    T xAvg = 0.5*(nX + 1), yAvg = 0.5*(nY + 1);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
        T I = i - xAvg, J = j - yAvg;
        xNull(i, j)    = I           / sqrt(xNorm);
        yNull(i, j)    = J           / sqrt(yNorm);
        xyNull(i, j)   = I*J         / sqrt(xyNorm);
        x2y2Null(i, j) = (I+J)*(I-J) / sqrt(x2y2Norm);

    }

    // Variable initialization
    tol = tol*tol;
    T a, B;
    CaslArray2D<T> p(nX, nY, 1, 1), Ap(nX, nY);

    CaslArray2D<T> r(nX, nY);
    matMult(x0, r);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j){
        r(i, j) = b(i, j) - r(i, j);
        p(i, j) = r(i, j);
    }

    T rTr = dot(r,r);
    // Main loop
    for(int k = 0; k < niter; ++k){
        matMult(p, Ap);

        a = rTr / dot(p, Ap);
        T sum = 0;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            x0(i, j) += a * p(i, j);
            r(i, j)  -= a * Ap(i, j);

            sum += x0(i, j);
        }
        T average = sum / _n;
        T xProj = dot(x0, xNull);
        T yProj = dot(x0, yNull);
        T xyProj = dot(x0, xyNull);
        T x2y2Proj = dot(x0, x2y2Null);
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            x0(i, j) -= average;
            x0(i, j) -= xProj * xNull(i, j);
            x0(i, j) -= yProj * yNull(i, j);
            x0(i, j) -= xyProj * xyNull(i, j);
            x0(i, j) -= x2y2Proj * x2y2Null(i, j);
        }
        rTr = dot(r,r);

        cout << rTr << ", ";
        if (rTr < tol) break;

        B = -(dot(r, Ap) / dot(p, Ap));
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j)
            p(i, j) = r(i, j) + B * p(i, j);
    }
}

template<typename T, class Matrix>
void DPMatrixExtended2D<T, Matrix>::biconjGradQUAD(CaslArray2D<T>& x0, const CaslArray2D<T>& b, DPMatrix2D<T>& L, DPMatrix2D<T>& U, T tol, int niter) {
    // Variable initialization
    tol = tol * tol;

    CaslArray2D<T> p(nX, nY, 1, 1), r(nX, nY, 1, 1), y(nX, nY, 1, 1), Ks(nX, nY, 1, 1), Kt(nX, nY, 1, 1);
    CaslArray2D<T> rHat(nX, nY), v(nX, nY, 1, 1), t(nX, nY, 1, 1);

    matMult(x0, r);
    for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
        r(i, j) = b(i, j) - r(i, j);
        rHat(i, j) = r(i, j); // Initialize rHat as r
        p(i, j) = r(i, j);
    }

    T rho = dot(rHat, r);

    // Main loop
    for (int k = 0; k < niter; ++k) {
        luSolve(y, p, L, U, v);
        matMult(y, v);

        T alpha = rho / dot(rHat, v);

        // Update x0
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            x0(i, j) += alpha * p(i, j);
            r(i, j) -= alpha * v(i, j);
        }

        T rTr = dot(r, r);

        cout << rTr << "(1), ";
        // Check for convergence
        if (rTr < tol) break;

        // Compute omega
        luSolve(y, r, L, U, t);
        matMult(y, t);

        L.lSolve(Ks, r);
        L.lSolve(Kt, t);
        T omega = dot(Kt, Ks) / dot(Kt, Kt);

        // Update r and x
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            x0(i, j) += omega * y(i, j);
            r(i, j) -= omega * t(i, j);
        }

        rTr = dot(r, r);
        cout << rTr << "(2), ";
        if (rTr < tol) break;

        // Compute beta
        T rho2 = dot(rHat, r);
        T beta = (rho2 / rho) * (alpha / omega) ;
        rho = rho2;
        for (int i = 1; i <= nX; ++i) for (int j = 1; j <= nY; ++j) {
            p(i, j) = r(i, j) + beta*(p(i, j) - omega*v(i, j));
        }
    }
}
