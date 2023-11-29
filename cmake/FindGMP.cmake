# - Try to find the GMP libraries
# This module defines:
#  GMP_FOUND             - system has GMP lib
#  GMP_INCLUDE_DIR       - the GMP include directory
#  GMP_LIBRARIES_DIR     - directory where the GMP libraries are located
#  GMP_LIBRARIES         - Link these to use GMP

# TODO: support MacOSX

include(FindPackageHandleStandardArgs)

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  HINTS ENV GMP_INC_DIR
        ENV GMP_DIR
        $ENV{GMP_DIR}/include
  PATH_SUFFIXES include
  DOC "The directory containing the GMP header files"
)

find_library(GMP_LIBRARY_RELEASE NAMES gmp libgmp-10 gmp-10 mpir
  HINTS ENV GMP_LIB_DIR
        ENV GMP_DIR
        $ENV{GMP_DIR}/lib
  PATH_SUFFIXES lib
  DOC "Path to the Release GMP library"
)

find_library(GMP_LIBRARY_DEBUG NAMES gmpd gmp libgmp-10 gmp-10 mpir
  HINTS ENV GMP_LIB_DIR
        ENV GMP_DIR
        $ENV{GMP_DIR}/include
  PATH_SUFFIXES lib
  DOC "Path to the Debug GMP library"
)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  set(GMP_LIBRARIES ${GMP_LIBRARY_DEBUG})
else()
  set(GMP_LIBRARIES ${GMP_LIBRARY_RELEASE})
endif()

# Attempt to load a user-defined configuration for GMP if couldn't be found
if ( NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARIES)
  include( GMPConfig OPTIONAL )
endif()

find_package_handle_standard_args(GMP "DEFAULT_MSG" GMP_LIBRARIES GMP_INCLUDE_DIR)
