cmake_minimum_required(VERSION 3.5)

project(CMR
  VERSION 1.3.0
  LANGUAGES C CXX)

set(CMAKE_C_STANDARD 99)
option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})
option(GENERATORS "Compile matrix generators" OFF)
option(TESTS "Compile tests" ON)
message(STATUS "Build tests: " ${TESTS})

# Add cmake/ to CMAKE_MODULE_PATH.
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

if(GENERATORS)
  message(STATUS "Generators: ON")
else()
  message(STATUS "Generators: OFF")
endif()

if (MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wunused)
endif()

find_package(GMP)
set(CMR_WITH_GMP ${GMP_FOUND})

# Target for the CMR library.
add_library(cmr
  src/cmr/camion.c
  src/cmr/ctu.c
  src/cmr/dec.c
  src/cmr/densematrix.c
  src/cmr/element.c
  src/cmr/env.c
  src/cmr/hereditary_property.c
  src/cmr/matrix.c
  src/cmr/one_sum.c
  src/cmr/tu.c
  src/cmr/graph.c
  src/cmr/graphic.c
  src/cmr/hashtable.c
  src/cmr/heap.c
  src/cmr/equimodular.c
  src/cmr/linear_algebra.c
  src/cmr/listmatrix.c
  src/cmr/matroid.c
  src/cmr/network.c
  src/cmr/regular.c
  src/cmr/regularity.c
  src/cmr/regularity_partition.c
  src/cmr/regularity_graphic.c
  src/cmr/regularity_nested_minor_sequence.c
  src/cmr/regularity_onesum.c
  src/cmr/regularity_threesum.c
  src/cmr/regularity_r10.c
  src/cmr/regularity_series_parallel.c
  src/cmr/separation.c
  src/cmr/series_parallel.c
  src/cmr/sort.c
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(cmr PRIVATE -Wall -Wextra -pedantic)
endif()

# Target for the cmr-camion executable.
add_executable(cmr_camion
  src/main/camion_main.c) 

target_link_libraries(cmr_camion
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_camion PROPERTIES OUTPUT_NAME cmr-camion)

# Target for the cmr-tu.
add_executable(cmr_tu
  src/main/tu_main.c)
target_link_libraries(cmr_tu
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_tu PROPERTIES OUTPUT_NAME cmr-tu)

# Target for the cmr-ctu.
add_executable(cmr_ctu
  src/main/ctu_main.c)
target_link_libraries(cmr_ctu
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_ctu PROPERTIES OUTPUT_NAME cmr-ctu)

# Target for the cmr-graphic executable.
add_executable(cmr_graphic
  src/main/graphic_main.c)
target_link_libraries(cmr_graphic
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_graphic PROPERTIES OUTPUT_NAME cmr-graphic)

# Target for the cmr-equimodular.
add_executable(cmr_equimodular
  src/main/equimodular_main.c)
target_link_libraries(cmr_equimodular
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_equimodular PROPERTIES OUTPUT_NAME cmr-equimodular)

# Target for the cmr-network executable.
add_executable(cmr_network
  src/main/network_main.c) 
target_link_libraries(cmr_network
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_network PROPERTIES OUTPUT_NAME cmr-network)

# Target for the cmr-regular
add_executable(cmr_regular
  src/main/regular_main.c)
target_link_libraries(cmr_regular
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_regular PROPERTIES OUTPUT_NAME cmr-regular)

# Target for the cmr-series-parallel executable.
add_executable(cmr_series_parallel
  src/main/series_parallel_main.c)
target_link_libraries(cmr_series_parallel
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_series_parallel PROPERTIES OUTPUT_NAME cmr-series-parallel)

# Target for the cmr-matrix.
add_executable(cmr_matrix
  src/main/matrix_main.c)
target_link_libraries(cmr_matrix
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_matrix PROPERTIES OUTPUT_NAME cmr-matrix)

# Target for the cmr-k-ary
add_executable(cmr_k_ary
  src/main/k_ary_main.c)
target_link_libraries(cmr_k_ary
  PRIVATE
    CMR::cmr
    m
)
set_target_properties(cmr_k_ary PROPERTIES OUTPUT_NAME cmr-k-ary)

if(GENERATORS)
  # Target for cmr-generate-series-parallel
  add_executable(cmr_generate_series_parallel
    src/gen/series_parallel_gen.c)
  target_link_libraries(cmr_generate_series_parallel
    PRIVATE
      CMR::cmr
      m
  )
  set_target_properties(cmr_generate_series_parallel PROPERTIES OUTPUT_NAME cmr-generate-series-parallel)

  # Target for cmr-generate-graphic
  add_executable(cmr_generate_graphic
    src/gen/graphic_gen.c)
  target_link_libraries(cmr_generate_graphic
    PRIVATE
      CMR::cmr
      m
  )
  set_target_properties(cmr_generate_graphic PROPERTIES OUTPUT_NAME cmr-generate-graphic)

  # Target for cmr-generate-network
  add_executable(cmr_generate_network
    src/gen/network_gen.c)
  target_link_libraries(cmr_generate_network
    PRIVATE
      CMR::cmr
      m
  )
  set_target_properties(cmr_generate_network PROPERTIES OUTPUT_NAME cmr-generate-network)

  # Target for cmr-generate-random
  add_executable(cmr_generate_random
    src/gen/random_gen.c)
  target_link_libraries(cmr_generate_random
    PRIVATE
      CMR::cmr
      m
  )
  set_target_properties(cmr_generate_random PROPERTIES OUTPUT_NAME cmr-generate-random)

  # Target for cmr-perturb-random
  add_executable(cmr_perturb_random
    src/gen/random_perturb.c)
  target_link_libraries(cmr_perturb_random
    PRIVATE
      CMR::cmr
      m
  )
  set_target_properties(cmr_perturb_random PROPERTIES OUTPUT_NAME cmr-perturb-random)

  set(GENERATOR_EXECUTABLES cmr_generate_series_parallel cmr_generate_graphic cmr_generate_network cmr_generate_random
    cmr_perturb_random)

  find_package(GUROBI)
  if(GUROBI_FOUND)
    # Target for cmr-exract-gurobi
    add_executable(cmr_extract_gurobi
      src/gen/gurobi_extract.c)
    target_include_directories(cmr_extract_gurobi
      PRIVATE
        ${GUROBI_INCLUDE_DIRS}
      )
    target_link_libraries(cmr_extract_gurobi
      PRIVATE
        CMR::cmr
        ${GUROBI_LIBRARIES}
        m
    )
    set_target_properties(cmr_extract_gurobi PROPERTIES OUTPUT_NAME cmr-extract-gurobi)
    set(GENERATOR_EXECUTABLES ${GENERATOR_EXECUTABLES} cmr_extract_gurobi)
  endif()
else()
  set(GENERATOR_EXECUTABLES)
endif()

# Write compilation settings to cmr/config.h.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/cmr/config.h.in ${CMAKE_BINARY_DIR}/cmr/config.h @ONLY)

# Write export settings to cmr/export.h.
include(GenerateExportHeader)
generate_export_header(cmr EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/cmr/export.h)

# TODO: Why both things below?
# Add an alias so that library can be used inside the build tree.
add_library(CMR::cmr ALIAS cmr)

# Hide non-exported symbols in shared library.
set_target_properties(cmr PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(cmr PROPERTIES VISIBILITY_INLINES_HIDDEN 1)

# Set target properties.
target_include_directories(cmr
  PUBLIC
    $<INSTALL_INTERFACE:include> # <PREFIX>/include/ contains all installed headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # contains all regular headers.
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # contains all configured headers such as cmr/config.h and cmr/export.h.
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cmr/
)

target_compile_features(cmr PRIVATE cxx_auto_type)
target_compile_options(cmr PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall>)

if(Threads_FOUND)
  target_link_libraries(cmr
    PRIVATE
      Threads::Threads
  )
endif()

if(GMP_FOUND)
  target_link_libraries(cmr
    PRIVATE
      ${GMP_LIBRARIES}
  )
endif()

### Installation ###
include(GNUInstallDirs)

install(TARGETS
    cmr_camion
    cmr_ctu
    cmr_graphic
    cmr_equimodular
    cmr_matrix
    cmr_network
    cmr_regular
    cmr_series_parallel
    cmr_tu
    ${GENERATOR_EXECUTABLES}
  RUNTIME
    DESTINATION bin
)

install(TARGETS cmr
  EXPORT cmr-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install header files.
install(DIRECTORY include/cmr
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${CMAKE_BINARY_DIR}/cmr/config.h ${CMAKE_BINARY_DIR}/cmr/export.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cmr)

# Export the targets to a script.
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/cmr/)
install(EXPORT cmr-targets
  FILE
    CMRTargets.cmake
    NAMESPACE CMR::
    DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

# Create CMRConfig.cmake
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/CMRConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMRConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Create CMRConfigVersion.cmake.
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/CMRConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Install CMRConfig.cmake and CMRConfigVersion.cmake.
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/CMRConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/CMRConfigVersion.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)

# Write exported targets.
export(EXPORT cmr-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/CMRTargets.cmake NAMESPACE CMR::)

# Register package in user's package registry.
export(PACKAGE CMR)

add_subdirectory(doc)

if(TESTS)
  enable_testing()
  add_subdirectory(test)
endif()
