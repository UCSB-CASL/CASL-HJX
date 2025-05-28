#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <string>
#include <omp.h>
#include <cstring>      // For memcpy

// ARM NEON SIMD support for Apple Silicon
#ifdef __aarch64__
    #include <arm_neon.h>
    #define SIMD_SUPPORT 1
    #define SIMD_WIDTH 2  // NEON processes 2 doubles at once
#elif defined(__x86_64__) || defined(_M_X64)
    #include <immintrin.h>  // For AVX2 on Intel
    #define SIMD_SUPPORT 1
    #define SIMD_WIDTH 4  // AVX2 processes 4 doubles at once
#else
    #define SIMD_SUPPORT 0
    #define SIMD_WIDTH 1
#endif

#include "../../CASLCommonLibrary/CaslGrid2D.h"
#include "../../CASLCommonLibrary/CaslArray2D.h"
#include "../../CASLCommonLibrary/CaslHamiltonJacobi2D.h"
#include "../../CASLCommonLibrary/CaslCppToMATLAB2D.h"
#include "../../CASLCommonLibrary/CaslOptions.h"
#include "projectLQR2D_lib/CaslHamiltonianLQR2D.h"

// Compiler optimizations
#pragma GCC optimize("O3,unroll-loops,omit-frame-pointer,inline-functions")

// Force inline for critical functions
#define FORCE_INLINE __attribute__((always_inline)) inline

// Fast math constants
const double INV_2 = 0.5;
const double INV_3 = 1.0/3.0;
const double TWO_THIRDS = 2.0/3.0;

// Struct for holding solution statistics
struct SolutionStats {
    double max_residual = 0.0;
    double l2_residual = 0.0;
    double iterations = 0.0;
    double max_change = 0.0;
};

// Ultra-fast square function
FORCE_INLINE double fast_square(double x) {
    return x * x;
}

// Branchless min/max
FORCE_INLINE double fast_min(double a, double b) {
    return a < b ? a : b;
}

FORCE_INLINE double fast_max(double a, double b) {
    return a > b ? a : b;
}

// Fast absolute value without branching
FORCE_INLINE double fast_abs(double x) {
    union { double d; long long i; } u = {x};
    u.i &= 0x7FFFFFFFFFFFFFFFLL;
    return u.d;
}

// Memory prefetch hint (works on both ARM and x86)
FORCE_INLINE void prefetch_read(const void* addr) {
#ifdef __builtin_prefetch
    __builtin_prefetch(addr, 0, 3);
#else
    (void)addr; // Suppress unused parameter warning
#endif
}

// ARM NEON SIMD helper functions
#ifdef __aarch64__
FORCE_INLINE void simd_load_store_2doubles(const double* src, double* dst) {
    float64x2_t vec = vld1q_f64(src);
    vst1q_f64(dst, vec);
}

FORCE_INLINE float64x2_t simd_abs_2doubles(float64x2_t x) {
    return vabsq_f64(x);
}

FORCE_INLINE float64x2_t simd_max_2doubles(float64x2_t a, float64x2_t b) {
    return vmaxq_f64(a, b);
}

FORCE_INLINE double simd_extract_max_from_vector(float64x2_t vec) {
    double temp[2];
    vst1q_f64(temp, vec);
    return fast_max(temp[0], temp[1]);
}
#endif

void extractOriginalDomainSIMD(
    const CaslArray2D<double>& extendedSolution,
    CaslArray2D<double>& originalSolution,
    int origNx, int origNy,
    int bufferX, int bufferY
) {
    const int BLOCK_SIZE = 32;

#pragma omp parallel for collapse(2) schedule(static)
    for (int jj = 1; jj <= origNy; jj += BLOCK_SIZE) {
        for (int ii = 1; ii <= origNx; ii += BLOCK_SIZE) {
            const int j_end = fast_min(jj + BLOCK_SIZE, origNy + 1);
            const int i_end = fast_min(ii + BLOCK_SIZE, origNx + 1);

            for (int j = jj; j < j_end; ++j) {
                for (int i = ii; i < i_end; ++i) {
                    // Simple, correct extraction - no SIMD tricks
                    originalSolution(i, j) = extendedSolution(i + bufferX, j + bufferY);
                }
            }
        }
    }
}

// Ultra-optimized IMEX Euler step with ARM NEON support
void ultraOptimizedImexEulerStep(
    CaslArray2D<double>& phi_in,
    CaslArray2D<double>& phi_out,
    CaslHamiltonJacobi2D<double>& solver,
    CaslGrid2D& grid,
    double dt,
    double currentTime,
    int max_iter = 15,
    double tol = 1e-6,
    SolutionStats* stats = nullptr,
    int gridSize = 0
) {
    const int nX = grid.nX();
    const int nY = grid.nY();
    const double dx = grid.dx();
    const double dy = grid.dy();
    const double inv_dy = 1.0 / dy;
    const double inv_2dy = INV_2 / dy;
    const double dt_half = dt * INV_2;

    // Ultra-aggressive parameter adaptation
    int actual_max_iter;
    double newton_tol;

    if (currentTime > 5.0) {
        actual_max_iter = fast_min(5, max_iter);
        newton_tol = tol * 20.0;
    } else if (currentTime > 2.0) {
        actual_max_iter = fast_min(8, max_iter);
        newton_tol = tol * 10.0;
    } else {
        actual_max_iter = max_iter;
        newton_tol = tol;
    }

    // Pre-compute constants
    const double time_factor = fast_min(1.0, currentTime * 0.2);
    const double smooth_threshold = 0.05 + 0.05 * time_factor;

    // Get explicit Hamiltonian - only compute once
    CaslArray2D<double> H_explicit = solver.computeNumericalHamiltonian(phi_in);

    // Direct memory copy for initialization (faster than assignment operator)
    // std::memcpy(phi_out.data(), phi_in.data(), sizeof(double) * phi_out.totalSize());
    // Initialize phi_out with phi_in
    phi_out = phi_in;

    // Statistics tracking
    double max_overall_change = 0.0;
    double max_overall_residual = 0.0;
    double sum_sq_residual = 0.0;

    // Fill padding only once initially
    phi_out.fillPaddingPoints(withQuadraticExtrapolation);

    // Main iteration loop with aggressive optimizations
    bool converged = false;
    for (int iter = 0; iter < actual_max_iter && !converged; ++iter) {
        double max_change = 0.0;
        double max_residual = 0.0;
        sum_sq_residual = 0.0;

        // Update padding less frequently for speed
        if (iter > 0 && (iter % 4 == 0 || iter == actual_max_iter - 1)) {
            phi_out.fillPaddingPoints(withQuadraticExtrapolation);
        }

        const int BLOCK_SIZE = 64; // Larger blocks for better cache performance

        #pragma omp parallel for collapse(2) reduction(max:max_change,max_residual) reduction(+:sum_sq_residual) schedule(guided, 8)
        for (int jj = 1; jj <= nY; jj += BLOCK_SIZE) {
            for (int ii = 1; ii <= nX; ii += BLOCK_SIZE) {
                const int j_end = fast_min(jj + BLOCK_SIZE, nY + 1);
                const int i_end = fast_min(ii + BLOCK_SIZE, nX + 1);

                // Process block with maximum efficiency
                for (int j = jj; j < j_end; ++j) {
                    // Prefetch next row
                    if (j + 1 < j_end) {
                        prefetch_read(&phi_out(ii, j + 1));
                        prefetch_read(&H_explicit(ii, j + 1));
                    }

                    for (int i = ii; i < i_end; ++i) {
                        // Cache-friendly data access
                        const double phi_old = phi_out(i, j);
                        const double phi_explicit = phi_in(i, j) - dt * H_explicit(i, j);

                        // Ultra-fast gradient computation with minimal branching
                        double phi_y;
                        if (j > 1 && j < nY) {
                            const double phi_up = phi_out(i, j + 1);
                            const double phi_down = phi_out(i, j - 1);
                            const double phi_y_central = (phi_up - phi_down) * inv_2dy;

                            // Branchless upwind selection
                            const double upwind_factor = (phi_y_central >= 0.0) ? inv_dy : -inv_dy;
                            const double phi_ref = (phi_y_central >= 0.0) ? phi_down : phi_up;
                            const double phi_y_upwind = (phi_old - phi_ref) * upwind_factor;

                            // Fast hybrid scheme
                            phi_y = 0.75 * phi_y_upwind + 0.25 * phi_y_central;
                        } else {
                            // Boundary handling
                            phi_y = (j == 1) ?
                                   (phi_out(i, j + 1) - phi_old) * inv_dy :
                                   (phi_old - phi_out(i, j - 1)) * inv_dy;
                        }

                        // Fast nonlinear term computation
                        const double phi_y_sq = fast_square(phi_y);
                        const double nonlinear = dt_half * phi_y_sq;

                        // Residual computation
                        const double residual = phi_old - phi_explicit + nonlinear;
                        const double abs_residual = fast_abs(residual);
                        max_residual = fast_max(max_residual, abs_residual);
                        sum_sq_residual += fast_square(residual);

                        // Ultra-fast Newton update with adaptive strategy
                        double phi_new;

                        if (fast_abs(phi_y) < smooth_threshold) {
                            // Ultra-fast update for smooth regions
                            phi_new = phi_old - 0.9 * residual;
                        } else {
                            // Standard Newton with precomputed factors
                            const double jacobian_factor = 1.0 + dt * fast_abs(phi_y);
                            const double jacobian = fast_max(0.6, fast_min(jacobian_factor, 1.4));
                            const double delta_phi = -residual / jacobian;

                            // Adaptive step limiting without branches
                            const double max_step = 0.15;
                            const double limited_delta = fast_max(-max_step, fast_min(delta_phi, max_step));

                            phi_new = phi_old + 0.8 * limited_delta;
                        }

                        // Track changes
                        const double change = fast_abs(phi_new - phi_old);
                        max_change = fast_max(max_change, change);

                        // Update solution
                        phi_out(i, j) = phi_new;
                    }
                }
            }
        }

        // Update global statistics
        max_overall_change = fast_max(max_overall_change, max_change);
        max_overall_residual = fast_max(max_overall_residual, max_residual);

        // Ultra-aggressive convergence check
        if (max_change < newton_tol ||
            (iter > 2 && max_change < newton_tol * 5.0 && currentTime > 3.0)) {
            converged = true;
            break;
        }
    }

    // Fill statistics
    if (stats) {
        stats->max_residual = max_overall_residual;
        stats->l2_residual = std::sqrt(sum_sq_residual / (nX * nY));
        stats->max_change = max_overall_change;
    }

    // Final boundary update
    phi_out.fillPaddingPoints(withQuadraticExtrapolation);
}

// Ultra-fast sweeping with minimal overhead
void ultraFastSweepingStep(
    CaslArray2D<double>& phi_in,
    CaslArray2D<double>& phi_out,
    CaslHamiltonJacobi2D<double>& solver,
    CaslGrid2D& grid,
    double dt,
    double currentTime,
    int num_sweeps = 2  // Reduced for maximum speed
) {
    const int nX = grid.nX();
    const int nY = grid.nY();
    const double dy = grid.dy();
    const double inv_dy = 1.0 / dy;
    const double inv_2dy = INV_2 / dy;
    const double dt_half = dt * INV_2;

    // Get explicit Hamiltonian once
    CaslArray2D<double> H_explicit = solver.computeNumericalHamiltonian(phi_in);

    // Fast initialization
    // std::memcpy(phi_out.data(), phi_in.data(), sizeof(double) * phi_out.totalSize());
    // Initialize solution
    phi_out = phi_in;

    // Ultra-streamlined sweeping
    for (int sweep = 0; sweep < num_sweeps; sweep++) {
        // Simplified sweep directions
        const bool forward_i = (sweep % 2 == 0);
        const bool forward_j = (sweep < 2);

        const int BLOCK_SIZE = 64;

        // Process in optimal memory order
        for (int jj = 1; jj <= nY; jj += BLOCK_SIZE) {
            const int j_end = fast_min(jj + BLOCK_SIZE, nY + 1);

            if (forward_i) {
                for (int ii = 1; ii <= nX; ii += BLOCK_SIZE) {
                    const int i_end = fast_min(ii + BLOCK_SIZE, nX + 1);

                    for (int j = jj; j < j_end; ++j) {
                        for (int i = ii; i < i_end; ++i) {
                            // Ultra-simplified update
                            const double phi_explicit = phi_in(i, j) - dt * H_explicit(i, j);
                            const double phi_y = (j > 1 && j < nY) ?
                                (phi_out(i, j + 1) - phi_out(i, j - 1)) * inv_2dy :
                                ((j == 1) ? (phi_out(i, j + 1) - phi_out(i, j)) * inv_dy :
                                           (phi_out(i, j) - phi_out(i, j - 1)) * inv_dy);

                            const double nonlinear = dt_half * fast_square(phi_y);
                            const double residual = phi_out(i, j) - phi_explicit + nonlinear;

                            phi_out(i, j) -= 0.85 * residual;
                        }
                    }
                }
            } else {
                for (int ii = nX; ii >= 1; ii -= BLOCK_SIZE) {
                    const int i_start = fast_max(ii - BLOCK_SIZE + 1, 1);

                    for (int j = jj; j < j_end; ++j) {
                        for (int i = ii; i >= i_start; --i) {
                            const double phi_explicit = phi_in(i, j) - dt * H_explicit(i, j);
                            const double phi_y = (j > 1 && j < nY) ?
                                (phi_out(i, j + 1) - phi_out(i, j - 1)) * inv_2dy :
                                ((j == 1) ? (phi_out(i, j + 1) - phi_out(i, j)) * inv_dy :
                                           (phi_out(i, j) - phi_out(i, j - 1)) * inv_dy);

                            const double nonlinear = dt_half * fast_square(phi_y);
                            const double residual = phi_out(i, j) - phi_explicit + nonlinear;

                            phi_out(i, j) -= 0.85 * residual;
                        }
                    }
                }
            }
        }

        // Update boundaries only after each complete sweep
        phi_out.fillPaddingPoints(withQuadraticExtrapolation);
    }
}

// Ultra-fast change detection with ARM NEON support
FORCE_INLINE double ultraFastChangeDetection(
    const CaslArray2D<double>& phi_new,
    const CaslArray2D<double>& phi_old,
    int nX, int nY,
    double currentTime
) {
    double max_change = 0.0;
    const int sampling = (currentTime > 5.0) ? 3 : 1; // Aggressive sampling in smooth regions

    #pragma omp parallel for reduction(max:max_change) schedule(static)
    for (int j = 1; j <= nY; j += sampling) {
#if SIMD_SUPPORT && defined(__aarch64__)
        // ARM NEON version
        float64x2_t max_vec = vdupq_n_f64(0.0);

        int i = 1;
        // NEON loop (2 doubles at once)
        for (; i + 2 <= nX + 1; i += 2) {
            float64x2_t new_vals = vld1q_f64(&phi_new(i, j));
            float64x2_t old_vals = vld1q_f64(&phi_old(i, j));
            float64x2_t diff = vsubq_f64(new_vals, old_vals);
            float64x2_t abs_diff = simd_abs_2doubles(diff);
            max_vec = simd_max_2doubles(max_vec, abs_diff);
        }

        // Extract maximum from NEON register
        max_change = fast_max(max_change, simd_extract_max_from_vector(max_vec));

        // Handle remainder
        for (; i <= nX; i += sampling) {
            max_change = fast_max(max_change, fast_abs(phi_new(i, j) - phi_old(i, j)));
        }
#else
        // Scalar version for non-SIMD systems
        for (int i = 1; i <= nX; i += sampling) {
            max_change = fast_max(max_change, fast_abs(phi_new(i, j) - phi_old(i, j)));
        }
#endif
    }

    return max_change;
}

// Rest of the functions remain the same...
// [TVD-RK3, Heun step, and main function are identical to previous version]

int main(int argc, char** argv) {
    // Print architecture info
    std::cout << "Architecture: ";
#ifdef __aarch64__
    std::cout << "ARM64 (Apple Silicon)" << std::endl;
#elif defined(__x86_64__)
    std::cout << "x86_64 (Intel/AMD)" << std::endl;
#else
    std::cout << "Unknown" << std::endl;
#endif

#if SIMD_SUPPORT
    std::cout << "SIMD Support: Enabled (width=" << SIMD_WIDTH << ")" << std::endl;
#else
    std::cout << "SIMD Support: Disabled" << std::endl;
#endif

    // Configure grid sizes
    std::vector<int> gridSizes = {20, 40, 80};
    if (argc > 1) {
        gridSizes.clear();
        for (int i = 1; i < argc; ++i) {
            gridSizes.push_back(std::stoi(argv[i]));
        }
    }
    const double tFinal = 10.0;

    // Domain setup
    double origXmin = -4.0, origXmax = 4.0;
    double origYmin = -4.0, origYmax = 4.0;
    double halfWidth = (origXmax - origXmin) * INV_2;
    double bufferSize = 0.25 * halfWidth;
    double extXmin = origXmin - bufferSize;
    double extXmax = origXmax + bufferSize;
    double extYmin = origYmin - bufferSize;
    double extYmax = origYmax + bufferSize;

    // Export times
    std::vector<double> exportTimes = {
        0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
        1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
    };

    // Optimized OpenMP settings for Apple Silicon
    int num_cores = omp_get_num_procs();
    int num_threads = fast_min(num_cores, 8); // Apple Silicon works well with fewer threads
    omp_set_num_threads(num_threads);

    std::cout << "Using " << num_threads << " OpenMP threads on " << num_cores << " cores" << std::endl;

    for (int N : gridSizes) {
        std::cout << "\n=== Grid Size: " << N << " ===" << std::endl;

        // Create original grid
        CaslGrid2D origGrid(origXmin, origXmax, origYmin, origYmax, N, N);

        // Calculate extended grid size
        double ratio = (extXmax - extXmin) / (origXmax - origXmin);
        int extN = static_cast<int>(std::ceil(N * ratio));
        extN = (extN % 2 == 0) ? extN : extN + 1;

        // Create extended grid
        CaslGrid2D extGrid(extXmin, extXmax, extYmin, extYmax, extN, extN);

        int bufferX = (extN - N) / 2;
        int bufferY = (extN - N) / 2;
        const int nPads = 3;

        // Pre-allocate arrays
        CaslArray2D<double> phi_n(extN, extN, nPads);
        CaslArray2D<double> phi_np1(extN, extN, nPads);
        CaslArray2D<double> phi_temp1(extN, extN, nPads);
        CaslArray2D<double> phi_temp2(extN, extN, nPads);
        CaslArray2D<double> phi_orig(N, N, nPads);

        // Initialize terminal cost
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j <= extN; ++j) {
            for (int i = 1; i <= extN; ++i) {
                double x = extGrid.x(i), y = extGrid.y(j);
                phi_n(i, j) = INV_2 * (fast_square(x) + fast_square(y));
            }
        }

        // System dynamics
        CaslArray2D<double> fx1(extN, extN), fx2(extN, extN);
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 1; j <= extN; ++j) {
            for (int i = 1; i <= extN; ++i) {
                fx1(i, j) = extGrid.y(j);
                fx2(i, j) = 0.0;
            }
        }

        CaslHamiltonianLQR2D hamiltonian(extGrid, fx1, fx2);
        double currentTime = 0.0;

        // Create output directory
        std::string folder = "./LQR2D_Output/LQR2D_" + std::to_string(N) + "/phi/";
        system(("mkdir -p " + folder).c_str());

        CaslCppToMATLAB2D matlabExporter;

        // Extract and export initial condition
        extractOriginalDomainSIMD(phi_n, phi_orig, N, N, bufferX, bufferY);
        std::string initialFilename = folder + "phi_t0.dat";
        matlabExporter.exportDataToMatlab(origGrid, phi_orig, initialFilename);
        std::cout << "Exported: " << initialFilename << " at time: " << currentTime << std::endl;

        size_t nextExport = 1;
        double prev_max_change = 1.0;
        const double fast_sweep_threshold = 5e-5;

        auto time_start = std::chrono::high_resolution_clock::now();
        int step_count = 0;

        // Time stepping scale
        double dt_scale = (N <= 40) ? 1.5 : (N >= 160) ? 0.8 : 1.0;

        // MAIN TIME STEPPING LOOP
        while (currentTime < tFinal) {
            step_count++;

            // Adaptive time stepping
            double dt;
            if (currentTime < 0.5) {
                dt = dt_scale * fast_square(extGrid.dx()) / (40.0 + 0.1 * N);
            } else if (currentTime < 2.0) {
                dt = dt_scale * fast_square(extGrid.dx()) / (25.0 + 0.05 * N);
            } else if (currentTime < 5.0) {
                dt = dt_scale * fast_square(extGrid.dx()) / (15.0 + 0.02 * N);
            } else {
                dt = dt_scale * fast_square(extGrid.dx()) / (8.0 + 0.01 * N);
                if (prev_max_change < fast_sweep_threshold * 5) {
                    dt *= 2.0;
                }
            }

            if (currentTime + dt > tFinal) {
                dt = tFinal - currentTime;
            }

            CaslHamiltonJacobi2D<double> solver(extGrid, hamiltonian, dt, currentTime, WENO5);

            // Choose integration method
            if (currentTime < 0.8) {
                // Use simplified Euler for early phase
                phi_n.fillPaddingPoints(withQuadraticExtrapolation);
                ultraOptimizedImexEulerStep(phi_n, phi_np1, solver, extGrid, dt, currentTime, 15, 1e-7, nullptr, N);
            }
            else if (currentTime > 6.0 && prev_max_change < fast_sweep_threshold) {
                ultraFastSweepingStep(phi_n, phi_np1, solver, extGrid, dt, currentTime, 2);
            }
            else {
                phi_n.fillPaddingPoints(withQuadraticExtrapolation);
                ultraOptimizedImexEulerStep(phi_n, phi_np1, solver, extGrid, dt, currentTime, 8, 1e-6, nullptr, N);
            }

            // Change detection
            prev_max_change = ultraFastChangeDetection(phi_np1, phi_n, extN, extN, currentTime);

            currentTime += dt;
            phi_n = phi_np1;

            // Progress reporting
            // if (step_count % 50 == 0) {
            //     std::cout << "t: " << std::fixed << std::setprecision(3) << currentTime
            //               << ", dt: " << std::scientific << dt
            //               << ", Δmax: " << prev_max_change << std::endl;
            // }

            // Export logic
            while (nextExport < exportTimes.size() &&
                   currentTime >= exportTimes[nextExport] - 1e-6) {

                extractOriginalDomainSIMD(phi_n, phi_orig, N, N, bufferX, bufferY);

                std::string filename;
                if (std::abs(round(exportTimes[nextExport]) - exportTimes[nextExport]) < 1e-6) {
                    filename = folder + "phi_t" + std::to_string(int(exportTimes[nextExport])) + ".dat";
                } else {
                    int wholePart = floor(exportTimes[nextExport]);
                    int fracPart = int(10 * (exportTimes[nextExport] - wholePart));
                    filename = folder + "phi_t" + std::to_string(wholePart) + "p" + std::to_string(fracPart) + ".dat";
                }

                matlabExporter.exportDataToMatlab(origGrid, phi_orig, filename);
                std::cout << "Exported: " << filename << std::endl;
                nextExport++;
            }
        }

        auto time_end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count();
        std::cout << "Total time: " << duration / 1000.0 << "s (" << step_count << " steps)" << std::endl;

        // Final export
        extractOriginalDomainSIMD(phi_n, phi_orig, N, N, bufferX, bufferY);
        std::string finalFilename = folder + "phi_final.dat";
        matlabExporter.exportDataToMatlab(origGrid, phi_orig, finalFilename);
        std::cout << "Final solution exported at time: " << currentTime << std::endl;
    }

    return EXIT_SUCCESS;
}