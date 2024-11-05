cmake_minimum_required(VERSION 3.15)
project(projectLQR2D)

# Find Eigen3 package
find_package(Eigen3 REQUIRED)

# Create library target for LQR-specific code
add_library(${PROJECT_NAME}_lib
        projectLQR2D_lib/CASLHamiltonianLQR2D.cpp
        projectLQR2D_lib/CASLLQRSystemDynamics.cpp
)

# Link Eigen3 to your library
target_link_libraries(${PROJECT_NAME}_lib
        PUBLIC
        CASLCommonLibrary
        Eigen3::Eigen  # Link Eigen3
)

# Create executable
add_executable(${PROJECT_NAME} main.cpp)

# Link libraries
target_link_libraries(${PROJECT_NAME}
        PRIVATE
        ${PROJECT_NAME}_lib
)

# Include directories
target_include_directories(${PROJECT_NAME}_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/projectLQR2D_lib
        ${EIGEN3_INCLUDE_DIR}  # Include Eigen3 directory
)

# Optionally include Eigen3 in the executable's include directories
target_include_directories(${PROJECT_NAME}
        PRIVATE
        ${EIGEN3_INCLUDE_DIR}
)
