find_package(GTest)
if(NOT GTEST_FOUND)
   find_package(GTestSources)
   if(NOT GTESTSOURCES_FOUND)
      find_package(GTestGit)
   endif()
endif()

# Target for the googletest.
add_executable(tu_gtest
  common.cpp
#  test_incidence.cpp
#  test_matrix.cpp
#  test_preprocessing.cpp
  test_main.cpp)

# Add tests for non-exported functions only for static library.
if(NOT SHARED)
  target_sources(tu_gtest
    PRIVATE
    test_tasks.cpp)
  target_include_directories(tu_gtest
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> # contains private headers.
)

endif()

# Configure tu_gtest target.
target_compile_features(tu_gtest PRIVATE cxx_auto_type)
target_link_libraries(tu_gtest gtest_main TU::tu)
   
include(GoogleTest)
gtest_discover_tests(tu_gtest)
