# CASLProjects/projectLaplacian2D_1D/CMakeLists.txt
cmake_minimum_required(VERSION 3.15)
project(projectReinitialization)

# Create library target for Laplacian-specific code
add_library(${PROJECT_NAME}_lib
        projectReinitialization_lib/CaslHamiltonianReinitialization2D.cpp
)

# Create executable
add_executable(${PROJECT_NAME} main.cpp)

# Link libraries
target_link_libraries(${PROJECT_NAME}_lib
        PUBLIC
        CASLCommonLibrary
)

target_link_libraries(${PROJECT_NAME}
        PRIVATE
        ${PROJECT_NAME}_lib
)

# Include directories
target_include_directories(${PROJECT_NAME}_lib
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/projectReinitialization_lib
)