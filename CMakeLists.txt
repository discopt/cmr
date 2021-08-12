cmake_minimum_required(VERSION 3.5)

project(CMR
  VERSION 1.9.0
  LANGUAGES C CXX)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})
option(THREADS "Support for parallelization" ON)
option(GENERATORS "Compile matrix generators" OFF)

# Add cmake/ to CMAKE_MODULE_PATH.
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

# Set default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release"
    CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Dependency: THREADS
if (THREADS)
  find_package(Threads)
  set(CMR_WITH_THREADS ${Threads_FOUND})
else()
  set(CMR_WITH_THREADS)
endif()

if(CMR_WITH_THREADS)
  message(STATUS "Parallelization: ON")
else()
  message(STATUS "Parallelization: OFF")
endif()

if(GENERATORS)
  message(STATUS "Generators: ON")
else()
  message(STATUS "Generators: OFF")
endif()

# Target for the CMR library.
add_library(cmr
  src/cmr/decomposition.c
  src/cmr/element.c
  src/cmr/env.c
  src/cmr/matrix.c
  src/cmr/one_sum.c
  src/cmr/sign.c
  src/cmr/tu.c
  src/cmr/determinant.cpp
  src/cmr/ghouila_houri.cpp
  src/cmr/graph.c
  src/cmr/graphic.c
  src/cmr/hashtable.c
  src/cmr/heap.c
  src/cmr/logger.cpp
  src/cmr/matrix.cpp
  src/cmr/matroid_decomposition.cpp
  src/cmr/matroid_graph.cpp
  src/cmr/nested_minor_sequence.cpp
  src/cmr/network.c
  src/cmr/regular.c
  src/cmr/regular_onesum.c
  src/cmr/regular_simplesums.c
  src/cmr/separation.cpp
  src/cmr/series_parallel.c
  src/cmr/sort.c
  src/cmr/total_unimodularity.cpp
  src/cmr/unimodularity.cpp
  src/cmr/zero_plus_minus_one.cpp
)

# Target for the tu-test binary.
add_executable(tu_test
  src/test/main.cpp)
target_link_libraries(tu_test
  PRIVATE
    CMR::cmr
)
set_target_properties(tu_test PROPERTIES OUTPUT_NAME tu-test)

# Target for the cmr-test-new binary.
add_executable(cmr_totally_unimodular
  src/test/totally_unimodular_main.c)
target_link_libraries(cmr_totally_unimodular
  PRIVATE
    CMR::cmr
)
set_target_properties(cmr_totally_unimodular PROPERTIES OUTPUT_NAME cmr-totally-unimodular)

# Target for the cmr-graphic binary.
add_executable(cmr_graphic
   src/test/graphic_main.c)
target_link_libraries(cmr_graphic
   PRIVATE
      CMR::cmr
)
set_target_properties(cmr_graphic PROPERTIES OUTPUT_NAME cmr-graphic)

# Target for the cmr-series-parallel binary.
add_executable(cmr_series_parallel
  src/test/series_parallel_main.c)
target_link_libraries(cmr_series_parallel
  PRIVATE
    CMR::cmr
)
set_target_properties(cmr_series_parallel PROPERTIES OUTPUT_NAME cmr-series-parallel)

# Target for the tu-convert-matrix.
add_executable(cmr_convert_matrix
  src/test/convert_matrix_main.c)
target_link_libraries(cmr_convert_matrix
  PRIVATE
    CMR::cmr
)
set_target_properties(cmr_convert_matrix PROPERTIES OUTPUT_NAME cmr-convert-matrix)

if(GENERATORS)
  add_executable(cmr_generate_series_parallel
    src/gen/series_parallel_gen.c)
  target_link_libraries(cmr_generate_series_parallel
    PRIVATE
      CMR::cmr
  )
  set_target_properties(cmr_generate_series_parallel PROPERTIES OUTPUT_NAME cmr-generate-series-parallel)
endif()

# Target for the random_network
#add_executable(random_network
#   src/test/random_network.c)
#target_link_libraries(random_network
#   PRIVATE
#      CMR::cmr
#)
#set_target_properties(random_network PROPERTIES OUTPUT_NAME random-network)

# Write compilation settings to cmr/config.h.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/cmr/config.h.in ${CMAKE_BINARY_DIR}/cmr/config.h @ONLY)

# Write export settings to cmr/export.h.
include(GenerateExportHeader)
generate_export_header(cmr EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/cmr/export.h)

# TODO: Why both things below?
# Add an alias so that library can be used inside the build tree.
add_library(CMR::cmr ALIAS cmr)

# This is required so that the exported target has the name CMR and not just cmr.
set_target_properties(cmr PROPERTIES EXPORT_NAME CMR)

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

### Installation ###
include(GNUInstallDirs)

install(TARGETS tu_test
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

install(DIRECTORY include
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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

enable_testing()
add_subdirectory(test)

