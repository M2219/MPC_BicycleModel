cmake_minimum_required(VERSION 3.10)
project(myproject)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find dependencies
find_package(Eigen3 REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(Python3 COMPONENTS Development REQUIRED)

# Include directories
include_directories(
    ${EIGEN3_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}  # Ensure headers in project root are found
)

# Source files
set(SRC_FILES
    main.cpp
    mpc_controller.cpp
    mpc_controller.hpp
    mpc_utils.cpp
    mpc_utils.hpp
    cubic_planner.cpp
    cubic_planner.hpp
)

# Create executable
add_executable(mpc ${SRC_FILES})

# Link Libraries (Ensure Python and OSQP are properly linked)
target_link_libraries(mpc PRIVATE OsqpEigen::OsqpEigen ${Python3_LIBRARIES})

# Additional Compiler Flags
target_compile_options(mpc PRIVATE -Wall -Wextra -O3 -fPIC)

