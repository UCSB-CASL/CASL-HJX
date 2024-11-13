# CASLProjects/projectLaplacian2D_1D/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(projectStochasticHH2D)

# Create library target for Laplacian-specific code
add_library(${PROJECT_NAME}_lib
        projectStochasticHH2D_lib/CaslHamiltonianHHModel.cpp
        projectStochasticHH2D_lib/CaslHHNeuronModel.cpp
)

# Create executable
add_executable(${PROJECT_NAME} main.cpp)

# Link libraries
# Our library needs CASLCommonLibrary and exposes it in its interface
target_link_libraries(${PROJECT_NAME}_lib
        PUBLIC              # Anyone using our lib will also need CASLCommonLibrary
        CASLCommonLibrary   # because our public headers might use it
)

# Our executable just uses our library internally
target_link_libraries(${PROJECT_NAME}
        PRIVATE                        # Nobody will link against our executable
        ${PROJECT_NAME}_lib   # so we keep this dependency private
)

# Include directories
target_include_directories(${PROJECT_NAME}_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/projectStochasticHH2D_lib
)