find_package(GTest)
if(NOT GTEST_FOUND)
   find_package(GTestSources)
   if(NOT GTESTSOURCES_FOUND)
      find_package(GTestGit)
   endif()
endif()

# Target for the googletest.
add_executable(cmr_gtest
  common.c
  test_camion.cpp
  test_ctu.cpp
  test_graph.cpp
  test_graphic.cpp
  test_hashtable.cpp
  test_matrix.cpp
  test_main.cpp
  test_network.cpp
  test_regular.cpp
  test_separation.cpp
  test_series_parallel.cpp)

# Add tests for non-exported functions only for static library.
if(NOT SHARED)
  target_sources(cmr_gtest
    PRIVATE
    test_one_sum.cpp
    )
  target_include_directories(cmr_gtest
    PRIVATE
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> # contains private headers.
)

endif()

# Configure cmr_gtest target.
target_compile_features(cmr_gtest PRIVATE cxx_auto_type)
target_link_libraries(cmr_gtest gtest_main CMR::cmr)
   
include(GoogleTest)
gtest_discover_tests(cmr_gtest)
