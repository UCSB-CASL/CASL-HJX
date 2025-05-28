cmake_minimum_required(VERSION 3.15)
project(projectLQR2D)

# Set C++ standard to C++17 for better optimization support
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# Detect CPU architecture for optimal flags
include(CheckCXXCompilerFlag)

# Function to check and add compiler flags
function(add_compiler_flag flag)
    string(REPLACE "-" "_" flag_var ${flag})
    string(REPLACE "=" "_" flag_var ${flag_var})
    string(REPLACE "+" "_" flag_var ${flag_var})
    string(TOUPPER ${flag_var} flag_var)
    set(flag_var "HAVE_${flag_var}")

    check_cxx_compiler_flag(${flag} ${flag_var})
    if(${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        message(STATUS "Added compiler flag: ${flag}")
    else()
        message(STATUS "Flag not supported: ${flag}")
    endif()
endfunction()

# Architecture-specific optimizations
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    message(STATUS "Detected ARM64 architecture (Apple Silicon)")

    # ARM-specific optimizations
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
    add_compiler_flag("-ftree-vectorize")
    add_compiler_flag("-fomit-frame-pointer")

    # ARM NEON SIMD support
    add_compiler_flag("-mfpu=neon")

    # ARM-specific memory optimizations
    add_compiler_flag("-falign-functions=16")
    add_compiler_flag("-falign-loops=16")

elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    message(STATUS "Detected x86_64 architecture")

    # x86-specific optimizations
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")

    # SIMD support for x86
    add_compiler_flag("-mavx2")
    add_compiler_flag("-mfma")
    add_compiler_flag("-msse4.2")
    add_compiler_flag("-ftree-vectorize")
    add_compiler_flag("-fomit-frame-pointer")
    add_compiler_flag("-fprefetch-loop-arrays")

    # x86-specific alignment
    add_compiler_flag("-falign-functions=32")
    add_compiler_flag("-falign-loops=32")
else()
    message(STATUS "Unknown architecture: ${CMAKE_SYSTEM_PROCESSOR}")
    # Generic optimizations
    add_compiler_flag("-ftree-vectorize")
    add_compiler_flag("-fomit-frame-pointer")
endif()

# Universal optimizations that work on all architectures
add_compiler_flag("-fno-trapping-math")
add_compiler_flag("-fno-signed-zeros")
add_compiler_flag("-fno-math-errno")

# Enable OpenMP SIMD support (works on both ARM and x86)
if(NOT APPLE OR NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compiler_flag("-fopenmp-simd")
endif()

# Release-specific optimizations
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast -funroll-loops -ffast-math")

# Compiler-specific optimizations
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    # GCC specific optimizations
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -floop-nest-optimize")
    add_compiler_flag("-funsafe-math-optimizations")

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # Clang specific optimizations
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mllvm -vectorize-loops")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mllvm -vectorize-slp")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mllvm -unroll-threshold=1000")
endif()

# Enable Link Time Optimization for maximum performance
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
if(CMAKE_BUILD_TYPE MATCHES "Release")
    add_compiler_flag("-flto")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
endif()

# Special handling for macOS + OpenMP
if(APPLE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # Use Homebrew OpenMP on macOS
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")

        # Try different possible OpenMP library locations
        find_library(OpenMP_omp_LIBRARY
                NAMES omp gomp iomp5
                PATHS
                /opt/homebrew/opt/libomp/lib
                /usr/local/opt/libomp/lib
                /opt/local/lib
                NO_DEFAULT_PATH
        )

        if(OpenMP_omp_LIBRARY)
            message(STATUS "Found OpenMP library: ${OpenMP_omp_LIBRARY}")
        else()
            message(WARNING "OpenMP library not found. Install with: brew install libomp")
        endif()
    endif()
endif()

# Create library target for LQR-specific code
add_library(${PROJECT_NAME}_lib
        projectLQR2D_lib/CASLHamiltonianLQR2D.cpp
        projectLQR2D_lib/CASLLQRSystemDynamics.cpp
)

# Apply optimizations to library as well
target_compile_options(${PROJECT_NAME}_lib PRIVATE
        $<$<CONFIG:Release>:-Ofast -funroll-loops -ffast-math>
)

# Link Eigen3 to your library
target_link_libraries(${PROJECT_NAME}_lib
        PUBLIC
        CASLCommonLibrary
        Eigen3::Eigen
)

# Create executable with optimized build flags
add_executable(${PROJECT_NAME} main.cpp
        projectLQR2D_lib/RiccatiSolver.h)

# Add compile definitions for optimized performance
target_compile_definitions(${PROJECT_NAME}
        PRIVATE
        OPTIMIZE_GRID_320
        AGGRESSIVE_NEWTON_CONVERGENCE
        USE_PREFETCHING
        ULTRA_FAST_MODE
        $<$<CONFIG:Release>:NDEBUG>
)

# Architecture-specific definitions
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE ARM64_OPTIMIZATIONS)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE X86_64_OPTIMIZATIONS)
endif()

# Compiler-specific definitions
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    target_compile_definitions(${PROJECT_NAME} PRIVATE GCC_OPTIMIZATIONS)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_definitions(${PROJECT_NAME} PRIVATE CLANG_OPTIMIZATIONS)
endif()

# Link libraries including OpenMP
target_link_libraries(${PROJECT_NAME}
        PRIVATE
        ${PROJECT_NAME}_lib
        OpenMP::OpenMP_CXX
)

# Include directories
target_include_directories(${PROJECT_NAME}_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/projectLQR2D_lib
        ${EIGEN3_INCLUDE_DIR}
)

target_include_directories(${PROJECT_NAME}
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
)

# Test for SIMD support based on architecture
if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    # Test for ARM NEON support
    include(CheckCXXSourceRuns)
    set(NEON_TEST_CODE "
    #include <arm_neon.h>
    int main() {
        float64x2_t a = vdupq_n_f64(1.0);
        float64x2_t b = vdupq_n_f64(2.0);
        float64x2_t c = vaddq_f64(a, b);
        return 0;
    }
    ")

    check_cxx_source_runs("${NEON_TEST_CODE}" HAVE_NEON_SUPPORT)
    if(HAVE_NEON_SUPPORT)
        message(STATUS "ARM NEON support: Available")
        target_compile_definitions(${PROJECT_NAME} PRIVATE NEON_SUPPORT)
    else()
        message(WARNING "ARM NEON support: Not available")
    endif()

elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    # Test for AVX2 support
    include(CheckCXXSourceRuns)
    set(AVX2_TEST_CODE "
    #include <immintrin.h>
    int main() {
        __m256d a = _mm256_set1_pd(1.0);
        __m256d b = _mm256_set1_pd(2.0);
        __m256d c = _mm256_add_pd(a, b);
        return 0;
    }
    ")

    check_cxx_source_runs("${AVX2_TEST_CODE}" HAVE_AVX2_SUPPORT)
    if(HAVE_AVX2_SUPPORT)
        message(STATUS "AVX2 support: Available")
        target_compile_definitions(${PROJECT_NAME} PRIVATE AVX2_SUPPORT)
    else()
        message(WARNING "AVX2 support: Not available")
    endif()
endif()

# Set optimization level based on build type
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    target_compile_options(${PROJECT_NAME} PRIVATE -O0 -g)
    message(STATUS "Debug build: Optimizations disabled")
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
    target_compile_options(${PROJECT_NAME} PRIVATE -O2 -g)
    message(STATUS "RelWithDebInfo build: Moderate optimizations with debug info")
else()
    message(STATUS "Release build: Maximum optimizations enabled")
endif()

# Apple Silicon specific runtime optimization suggestions
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        message(STATUS "")
        message(STATUS "=== APPLE SILICON OPTIMIZATION TIPS ===")
        message(STATUS "For maximum performance on Apple Silicon:")
        message(STATUS "export OMP_NUM_THREADS=8")
        message(STATUS "export OMP_SCHEDULE=guided,2")
        message(STATUS "export OMP_PROC_BIND=spread")
        message(STATUS "export OMP_PLACES=cores")
        message(STATUS "=======================================")
        message(STATUS "")
    endif()
endif()

# Output comprehensive build information
message(STATUS "")
message(STATUS "=== BUILD CONFIGURATION ===")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "System: ${CMAKE_SYSTEM_NAME}")
message(STATUS "Processor: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Base CXX flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Release CXX flags: ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "OpenMP support: ${OpenMP_FOUND}")
message(STATUS "Link-time optimization: ${CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE}")
message(STATUS "============================")
message(STATUS "")

# Add a custom target for performance testing
add_custom_target(performance_test
        COMMAND echo "Running performance test on ${CMAKE_SYSTEM_PROCESSOR}..."
        COMMAND ./${PROJECT_NAME} 40 80
        DEPENDS ${PROJECT_NAME}
        COMMENT "Running performance test with 40x40 and 80x80 grids"
)

# Add a custom target to show optimization flags
add_custom_target(show_flags
        COMMAND echo "Compiler flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}"
        COMMENT "Show all compiler optimization flags"
)