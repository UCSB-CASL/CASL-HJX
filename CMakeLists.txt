# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(CASLHJB2D
        VERSION 1.0
        DESCRIPTION "CASL Hamilton-Jacobi-Bellman 2D Solver"
        LANGUAGES CXX)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Compiler flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Modern C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# LLVM configuration for macOS
if(APPLE)
    if(EXISTS "/opt/homebrew/opt/llvm/bin/clang++")
        set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
        set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    else()
        message(WARNING "LLVM not found at specified path. Using default compiler.")
    endif()

    # Check macOS version for filesystem library
    execute_process(
            COMMAND sw_vers -productVersion
            OUTPUT_VARIABLE MACOS_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(MACOS_VERSION VERSION_LESS "10.15")
        link_libraries(c++fs)
    endif()
endif()

# Create the common library
add_library(CASLCommonLibrary
        CASLCommonLibrary/CASLArray2D.cpp
        CASLCommonLibrary/CaslGrid2D.h
        CASLCommonLibrary/CaslOptions.h
        CASLCommonLibrary/CASLCppToMatlab2D.cpp
        CASLCommonLibrary/CASLHamiltonian2D.cpp
        CASLCommonLibrary/CASLHamiltonJacobi.cpp
        CASLCommonLibrary/CASLHamiltonJacobi2D.cpp
        CASLCommonLibrary/CaslHamiltonian2D.cpp
        CASLCommonLibrary/CASLRK4.cpp
        CASLCommonLibrary/CASLSecondOrderDerivative2D_1D.cpp
        CASLCommonLibrary/DPMatrix2D.cpp
)

# Set properties for the common library
set_target_properties(CASLCommonLibrary PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Add include directory for the common library
target_include_directories(CASLCommonLibrary
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CASLCommonLibrary>
        $<INSTALL_INTERFACE:include>
)

# Add the projects
add_subdirectory(CASLProjects/projectDeterministicHH2D)
add_subdirectory(CASLProjects/projectStochasticHH2D)
add_subdirectory(CASLProjects/projectLQR2D)
add_subdirectory(CASLProjects/projectLaplacian2D_1D)

add_subdirectory(CASLProjects/projectDiffusion)
add_subdirectory(CASLProjects/projectAdvectionDiffusion)
add_subdirectory(CASLProjects/projectBurgers)
add_subdirectory(CASLProjects/projectAdvection)