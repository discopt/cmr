cmake_minimum_required(VERSION 3.5)

project(TU
   VERSION 2.0.0
   LANGUAGES C CXX)

option(SHARED "Build shared libraries" ON)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})
option(THREADS "Support for parallelization" ON)

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
  set(TU_WITH_THREADS ${Threads_FOUND})
else()
  set(TU_WITH_THREADS)
endif()

if(TU_WITH_THREADS)
  message(STATUS "Parallelization: ON")
else()
  message(STATUS "Parallelization: OFF")
endif()

# Target for the TU library.
add_library(tu
  src/tu/bixby_wagner.c
  src/tu/env.c
  src/tu/matrix.c
  src/tu/one_sum.c
  src/tu/sign.c
  src/tu/tu.c
  src/tu/determinant.cpp
  src/tu/ghouila_houri.cpp
  src/tu/graph.c
  src/tu/graphic.c
  src/tu/logger.cpp
  src/tu/matrix.cpp
  src/tu/matroid_decomposition.cpp
  src/tu/matroid_graph.cpp
  src/tu/nested_minor_sequence.cpp
  src/tu/regular.c
  src/tu/regular_dec.c
  src/tu/regular_onesum.c
  src/tu/regular_simplesums.c
  src/tu/separation.cpp
  src/tu/total_unimodularity.cpp
  src/tu/unimodularity.cpp
  src/tu/zero_plus_minus_one.cpp
)

# Target for the tu-test binary.
add_executable(tu_test
   src/test/main.cpp)
target_link_libraries(tu_test
   PRIVATE
      TU::tu
)
set_target_properties(tu_test PROPERTIES OUTPUT_NAME tu-test)

# Target for the random_network
add_executable(random_network
   src/test/random_network.c)
target_link_libraries(random_network
   PRIVATE
      TU::tu
)
set_target_properties(random_network PROPERTIES OUTPUT_NAME random-network)

# Write compilation settings to tu/config.h.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/tu/config.h.in ${CMAKE_BINARY_DIR}/tu/config.h @ONLY)

# Write export settings to tu/export.h.
include(GenerateExportHeader)
generate_export_header(tu EXPORT_FILE_NAME ${CMAKE_CURRENT_BINARY_DIR}/tu/export.h)

# TODO: Why both things below?
# Add an alias so that library can be used inside the build tree.
add_library(TU::tu ALIAS tu)

# This is required so that the exported target has the name TU and not just tu.
set_target_properties(tu PROPERTIES EXPORT_NAME TU)

# Hide non-exported symbols in shared library.
set_target_properties(tu PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(tu PROPERTIES VISIBILITY_INLINES_HIDDEN 1)

# Set target properties.
target_include_directories(tu
   PUBLIC
      $<INSTALL_INTERFACE:include> # <PREFIX>/include/ contains all installed headers.
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> # contains all regular headers.
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # contains all configured headers such as tu/config.h and tu/export.h.
   PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/src/tu/
)

target_compile_features(tu PRIVATE cxx_auto_type)
target_compile_options(tu PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall>)

if(Threads_FOUND)
  target_link_libraries(tu
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

install(TARGETS tu
   EXPORT tu-targets
   LIBRARY
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
   ARCHIVE
      DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Export the targets to a script.
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/tu/)
install(EXPORT tu-targets
   FILE
      TUTargets.cmake
      NAMESPACE TU::
      DESTINATION ${INSTALL_CONFIGDIR}
)

include(CMakePackageConfigHelpers)

# Create TUConfig.cmake
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/cmake/TUConfig.cmake.in
   ${CMAKE_CURRENT_BINARY_DIR}/TUConfig.cmake
   INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

# Create TUConfigVersion.cmake.
write_basic_package_version_file(
   ${CMAKE_CURRENT_BINARY_DIR}/TUConfigVersion.cmake
   VERSION ${PROJECT_VERSION}
   COMPATIBILITY SameMajorVersion
)

# Install TUConfig.cmake and TUConfigVersion.cmake.
install(FILES
   ${CMAKE_CURRENT_BINARY_DIR}/TUConfig.cmake
   ${CMAKE_CURRENT_BINARY_DIR}/TUConfigVersion.cmake
   DESTINATION ${INSTALL_CONFIGDIR}
)

# Write exported targets.
export(EXPORT tu-targets
   FILE ${CMAKE_CURRENT_BINARY_DIR}/TUTargets.cmake NAMESPACE TU::)

# Register package in user's package registry.
export(PACKAGE TU)

add_subdirectory(doc)

enable_testing()
add_subdirectory(test)

