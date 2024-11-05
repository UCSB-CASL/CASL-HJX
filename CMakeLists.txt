# Root CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(CASLHJB2D)

# Set default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Compiler flags for different build types
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the path to LLVM explicitly
set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")

# Add filesystem support for macOS
if(APPLE)
    link_libraries(c++fs)
endif()

# Create the common library directly here
add_library(CASLCommonLibrary
        CASLCommonLibrary/CASLArray2D.cpp
        CASLCommonLibrary/CaslGrid2D.h
        CASLCommonLibrary/CaslOptions.h
        CASLCommonLibrary/CASLCppToMatlab2D.cpp
        CASLCommonLibrary/CASLHamiltonian2D.cpp
        CASLCommonLibrary/CASLHamiltonJacobi.cpp
        CASLCommonLibrary/CASLHamiltonJacobi2D.cpp
        CASLCommonLibrary/CASLRK4.cpp
        CASLCommonLibrary/CASLSecondOrderDerivative2D_1D.cpp
)

# Add include directory for the common library
target_include_directories(CASLCommonLibrary
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/CASLCommonLibrary
)

# Add the projects
#add_subdirectory(CASLProjects/projectDeterministicHH2D)
#add_subdirectory(CASLProjects/projectStochasticHH2D)
add_subdirectory(CASLProjects/projectLQR2D)
add_subdirectory(CASLProjects/projectLaplacian2D_1D)