add_executable(foam_tests
  test_main.cpp
  test_mesh.cpp
  test_kinematics.cpp
  test_materials.cpp
  test_solvers.cpp
  test_io.cpp
)

target_link_libraries(foam_tests PRIVATE foam_core)
target_include_directories(foam_tests PRIVATE ${CMAKE_SOURCE_DIR}/external/doctest)

add_test(NAME foam_tests COMMAND foam_tests)
