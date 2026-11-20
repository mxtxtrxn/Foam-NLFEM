cmake_minimum_required(VERSION 3.16)
project(foam_nlfem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# --- Dependencies ------------------------------------------------------
# Eigen3 is the only hard external dependency for linear algebra.
# Prefer a system installation when available; otherwise use the vendored
# header-only copy in external/eigen (Eigen 3.4.0, MPL-2.0). This lets the
# project build on an offline desktop without apt-get, vcpkg or conan.
find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/external/eigen"
  )
  message(STATUS "Eigen3 not found on system: using vendored headers in external/eigen")
endif()

# --- Core library -------------------------------------------------------
# All reusable code (mesh, kinematics, materials, elements, assembly,
# solvers, io, simulation) lives in a static library that both the apps
# and the tests link against.
add_library(foam_core STATIC
  src/mesh/node.cpp
  src/mesh/mesh.cpp
  src/mesh/mesh_generators.cpp
  src/kinematics/deformation_gradient.cpp
  src/kinematics/strain_measures.cpp
  src/materials/ogden_hyperfoam.cpp
  src/materials/neo_hookean.cpp
  src/elements/quad4.cpp
  src/assembly/assembler.cpp
  src/solvers/linear_solver.cpp
  src/solvers/newton_raphson.cpp
  src/solvers/arc_length.cpp
  src/io/csv_writer.cpp
  src/io/vtk_writer.cpp
  src/simulation/problem.cpp
)

target_include_directories(foam_core
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(foam_core PUBLIC Eigen3::Eigen)

target_compile_options(foam_core PRIVATE
  -Wall -Wextra -Wpedantic
)

# --- Apps -----------------------------------------------------------
add_subdirectory(apps)

# --- Tests ------------------------------------------------------------
enable_testing()
add_subdirectory(tests)
